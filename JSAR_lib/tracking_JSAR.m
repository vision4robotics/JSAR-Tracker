
function [results] = tracking_JSAR(params)

%   Setting parameters for local use.

search_area_scale   = params.search_area_scale;
output_sigma_factor = params.output_sigma_factor;
features    = params.t_features;
video_path  = params.video_path;
s_frames    = params.s_frames;
pos         = floor(params.init_pos);
target_sz   = floor(params.wsize);
% gamma=params.gamma;
visualization  = params.visualization;
num_frames     = params.no_fram;
init_target_sz = target_sz;
featureRatio = params.t_global.cell_size;
global_feat_params = params.t_global;
% context_interp=params.context_interp;

%set the feature ratio to the feature-cell size
search_area = prod(init_target_sz * search_area_scale);
if search_area > params.max_image_sample_size
    CurrentSvFactor = sqrt(search_area /params.max_image_sample_size);
elseif search_area <params.min_image_sample_size
    CurrentSvFactor = sqrt(search_area /params.min_image_sample_size);
else
    CurrentSvFactor = 1.0;
end
CurrentArcFactor=1.0;
% target size at the initial scale
base_target_sz = target_sz / CurrentSvFactor;

% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);

% construct the label function- correlation output, 2D gaussian function,
% with a peak located upon the target
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg           = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg           = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs]     = ndgrid( rg,cg);
y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf           = fft2(y); %   FFT of y.
interp_sz = use_sz;

% construct cosine window
cos_window = single(hann(use_sz(1))*hann(use_sz(2))');
% construct consine window for scale estimation
window=single(hann(params.num_sv)*hann(params.num_arc)');
% Calculate feature dimension
try
    im = imread([video_path '/img/' s_frames{1}]);
catch
    try
        im = imread(s_frames{1});
    catch
        im = imread([video_path '/' s_frames{1}]);
    end
end
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end

% compute feature dimensionality
feature_dim = 0;
for n = 1:length(features)
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        feature_dim = feature_dim + features{n}.fparams.nDim;
    end
end


% Pre-computes the grid that is used for socre optimization
 ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
 kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
 sy=circshift(-floor((params.num_sv- 1)/2) : ceil((params.num_sv - 1)/2), [1, -floor((params.num_sv - 1)/2)]);
 sx=circshift(-floor((params.num_arc - 1)/2) : ceil((params.num_arc - 1)/2), [1, -floor((params.num_arc - 1)/2)]);
 newton_iterations = params.newton_iterations;

% parameters for sv and arc estimation
sv_sigma = sqrt(params.num_sv) *params.sv_sigma_factor;
arc_sigma=sqrt(params.num_arc)*params.arc_sigma_factor;
svg= (1:params.num_sv) - ceil(params.num_sv/2);
arcg=(1:params.num_arc) - ceil(params.num_arc/2);
[svs, arcs]     = ndgrid(svg,arcg);
ys = exp(-0.5 * ((svs.^2 / sv_sigma^2)+(arcs.^2/arc_sigma^2)));
ysf = fft2(ys);

ss = 1:params.num_sv;
as=1:params.num_arc;
svFactors = params.sv_step.^(ceil(params.num_sv/2) - ss);
arcFactors=params.arc_step.^(ceil(params.num_arc/2) - as);

if params.size_model_factor^2 * prod(target_sz) > params.size_model_max_area
    params.size_model_factor = sqrt(params.size_model_max_area/prod(target_sz));
end

if prod(target_sz) > params.size_model_max_area
    params.size_model_factor = sqrt(params.size_model_max_area/prod(target_sz));
end
size_model_sz = floor(target_sz* params.size_model_factor);

% set maximum and minimum scales
min_sv_factor = params.sv_step ^ ceil(log(max(20./sz)) / log(params.sv_step));
max_sv_factor = params.sv_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / (log(params.sv_step)));
min_arc_factor =1.5*1/search_area_scale^2*(init_target_sz(1)/init_target_sz(2));
max_arc_factor =search_area_scale^2/(1.5*init_target_sz(1)/init_target_sz(2));


% initialize the projection matrix (x,y,h,w)
rect_position = zeros(num_frames, 4);
time = 0;
reg_sz = floor(base_target_sz/featureRatio);

loop_frame = 1;

for frame = 1:numel(s_frames)
    %load image
    try
        im = imread([video_path '/img/' s_frames{frame}]);
    catch
        try
            im = imread([s_frames{frame}]);
        catch
            im = imread([video_path '/' s_frames{frame}]);
        end
    end
    if size(im,3) > 1 && colorImage == false
        im = im(:,:,1);
    end
   
    %do not estimate translation and scaling on the first frame, since we
    %just want to initialize the tracker there
    tic();
    if frame>1    
            pixel_template = get_pixels(im, pos, round(sz*CurrentSvFactor), sz);
            xtf = fft2(bsxfun(@times,get_features(pixel_template,features,global_feat_params),cos_window));
            responsef = sum(bsxfun(@times, conj(g_f), xtf), 3);
            % if we undersampled features, we want to interpolate the
            % response so it has the same size as the image patch
            responsef_padded = resizeDFT2(responsef, interp_sz);
            % response in the spatial domain
            response = ifft2(responsef_padded, 'symmetric');
            responsef_padded=fft2(response);
            % find maximum peak
            [disp_row, disp_col, ~] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, use_sz);
            % calculate translation
            translation_vec = round([disp_row, disp_col] * featureRatio * CurrentSvFactor);
            % update position
            pos = pos + translation_vec;
            xs = get_size_sample(im, pos, base_target_sz, CurrentSvFactor * svFactors, CurrentArcFactor*arcFactors,size_model_sz,window,frame);
            xsf = fft2(xs);
            scale_response = ifft2((sum(sf_num .* xsf, 3) ./ (sf_den + params.size_lambda)),'symmetric'); 
          
            % find the maximum scale response
            [recovered_sv,recovered_arc] = find(scale_response == max(scale_response(:)), 1);
            scale_vector=[recovered_sv-floor(size(scale_response,1)/2)-1,recovered_arc-floor(size(scale_response,2)/2)-1];
            % update the scale
            CurrentSvFactor = CurrentSvFactor * svFactors(recovered_sv);
            if CurrentSvFactor < min_sv_factor
                CurrentSvFactor = min_sv_factor;
            elseif CurrentSvFactor > max_sv_factor
                CurrentSvFactor = max_sv_factor;
            end            
            CurrentArcFactor=CurrentArcFactor*arcFactors(recovered_arc);
            if CurrentArcFactor < min_arc_factor
                CurrentArcFactor = min_arc_factor;
            elseif CurrentArcFactor > max_arc_factor
                CurrentArcFactor = max_arc_factor;
            end    
            
    end
    if frame==1||frame==5
    % extract training sample image region
    pixels = get_pixels(im,pos,round(sz*CurrentSvFactor),sz);
    % extract features and do windowing
    xf = fft2(bsxfun(@times,get_features(pixels,features,global_feat_params),cos_window));
    else 
        shift_sz=CurrentSvFactor * sz;
        shift_samp_pos = 2*pi * translation_vec ./shift_sz;
        xf = shift_sample(xtf, shift_samp_pos, kx', ky');
    end

    if  frame == 1
        [~,~,w]=init_regwindow(use_sz,reg_sz,params);
        g_pre= zeros(size(xf));
        mu = 0;
        beg_sz=reg_sz;
    else
        mu=params.init_mu;
        reg_sz(1)=round(beg_sz(1)*sqrt(CurrentArcFactor));
        reg_sz(2)=round(beg_sz(2)/sqrt(CurrentArcFactor));
        [~,~,w]=init_regwindow(use_sz,reg_sz,params);
    end
        g_f = single(zeros(size(xf)));
        h_f = g_f;
        l_f = h_f;
        gamma = 1;
        betha = 10;
        gamma_max = 10000;
        % ADMM solution    
        T = prod(use_sz);
        S_xx = sum(conj(xf) .* xf, 3);
        Sg_pre= sum(conj(xf) .* g_pre, 3);
        Sgx_pre= bsxfun(@times, xf, Sg_pre);
        iter = 1;
    %   ADMM
    while (iter <= params.admm_iterations)
            % subproblem g
            B = S_xx + T * (gamma + mu);
            Shx_f = sum(conj(xf) .* h_f, 3);
            Slx_f = sum(conj(xf) .* l_f, 3);
            g_f = ((1/(T*(gamma + mu)) * bsxfun(@times,  yf, xf)) - ((1/(gamma + mu)) * l_f) +(gamma/(gamma + mu)) * h_f) + (mu/(gamma + mu)) * g_pre - ...
                bsxfun(@rdivide,(1/(T*(gamma + mu)) * bsxfun(@times, xf, (S_xx .*  yf)) + (mu/(gamma + mu)) * Sgx_pre- ...
                (1/(gamma + mu))* (bsxfun(@times, xf, Slx_f)) +(gamma/(gamma + mu))* (bsxfun(@times, xf, Shx_f))), B);
            %   subproblem h
            lhd= T ./  (params.admm_lambda*w .^2 + gamma*T); 
            X=ifft2(gamma*(g_f + l_f));
            h=bsxfun(@times,lhd,X);
            h_f = fft2(h);
            %   update h
            l_f = l_f + (gamma * (g_f - h_f));
            %   update gamma
            gamma = min(betha* gamma,gamma_max );
            iter = iter+1;
    end
            % save the trained filters
            g_pre= g_f;
     %% Upadate Size
     if frame==1
        xs = get_size_sample(im, pos, base_target_sz, CurrentSvFactor * svFactors, CurrentArcFactor * arcFactors,size_model_sz,window,frame);
        xsf = fft2(xs);
     else
        if scale_vector(1)~=0||scale_vector(2)~=0
        shift_scale_pos = 2*pi * scale_vector ./[params.num_sv, params.num_arc];
        xsf = shift_sample(xsf, shift_scale_pos, sx', sy');
        end
     end   
    new_sf_num = ysf.*conj(xsf);
    new_sf_den = sum(xsf .* conj(xsf), 3);

    if frame == 1
        sf_den = new_sf_den;
        sf_num = new_sf_num;
    else
        sf_den = (1 - params.learning_rate_size) * sf_den + params.learning_rate_size * new_sf_den;
        sf_num = (1 - params.learning_rate_size) * sf_num + params.learning_rate_size * new_sf_num;
    end

    target_sz = base_target_sz * CurrentSvFactor;
    target_sz(1)=floor(target_sz(1)*sqrt(CurrentArcFactor));
    target_sz(2)=floor(target_sz(2)/sqrt(CurrentArcFactor));
    %save position and calculate FPS
    rect_position(loop_frame,:) = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
    
    time = time + toc();
     
%     
    %visualization
     if visualization == 1
        rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        figure(1);
        imshow(im);
        if frame == 1
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(12, 26, ['# Frame : ' int2str(loop_frame) ' / ' int2str(num_frames)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 12);
            hold off;
        else
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(12, 28, ['# Frame : ' int2str(loop_frame) ' / ' int2str(num_frames)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 12);
            text(12, 66, ['FPS : ' num2str(1/(time/loop_frame))], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 12);
            hold off;
         end
        drawnow
    end
    loop_frame = loop_frame + 1;
end
%   save resutls.
fps = loop_frame / time;
results.type = 'rect';
results.res = rect_position;
results.fps = fps;
end