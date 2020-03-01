% This function implements the ASRCF tracker.

function [results] =tracking_JSAR_Re(params)

num_frames     = params.no_fram;
newton_iterations=params.newton_iterations;
search_area_scale=params.search_area_scale;
global_feat_params = params.t_global;
featureRatio = params.t_global.cell_size;
search_area = prod(params.wsize * params.search_area_scale);
pos         = floor(params.init_pos);
target_sz   = floor(params.wsize);

[CurrentSvFactor,CurrentArcFactor,base_target_sz,reg_sz,sz,use_sz] = init_size(params,target_sz,search_area,featureRatio);
[y_0,cos_window] = init_gauss_win(params,base_target_sz,featureRatio,use_sz);
 y=y_0;
 yf=fft2(y);
[features,im,colorImage] = init_features(params);
[model,opts]=init_EdgeBox;
[ysf,svFactors,arcFactors,size_model_sz,min_sv_factor,max_sv_factor,min_arc_factor,max_arc_factor] = size_factor(params,target_sz,sz,base_target_sz,im);
% Pre-computes the grid that is used for score optimization
 ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
 kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
 sy=circshift(-floor((params.num_sv- 1)/2) : ceil((params.num_sv - 1)/2), [1, -floor((params.num_sv - 1)/2)]);
 sx=circshift(-floor((params.num_arc - 1)/2) : ceil((params.num_arc - 1)/2), [1, -floor((params.num_arc - 1)/2)]);
% initialize the projection matrix (x,y,h,w)
rect_position = zeros(num_frames, 4);
time = 0;
loop_frame = 1;
redetection=0;
success=1;
range=5;
rd_update=1;
in_sight=1;
for frame = 1:num_frames
    im = load_image(params,frame,colorImage);
    tic();  
%% main loop
if success==1
    if frame > 1
        % detection the object in new frame 
        [xtf,pos,translation_vec,response] = run_detection(im,pos,sz,CurrentSvFactor,features,cos_window,g_f,global_feat_params,use_sz,ky,kx,newton_iterations,featureRatio);
        % search for the scale of object
        [xs,CurrentSvFactor,CurrentArcFactor,size_vector]  = search_size(sf_num,sf_den,im,pos,base_target_sz,CurrentSvFactor,CurrentArcFactor,svFactors,arcFactors,size_model_sz,min_sv_factor,max_sv_factor,min_arc_factor,max_arc_factor,params);
         xsf = fft2(xs);
    end
    % update the target_sz via currentScaleFactor
        target_sz = base_target_sz * CurrentSvFactor;
        target_sz(1)=floor(target_sz(1)*sqrt(CurrentArcFactor));
        target_sz(2)=floor(target_sz(2)/sqrt(CurrentArcFactor));
        if (pos(1)-target_sz(1)/2)<1||(pos(1)+target_sz(1)/2)>size(im,1)||(pos(2)-target_sz(2)/2)<1||(pos(2)+target_sz(2)/2)>size(im,2)
            in_sight=0;
        else 
            in_sight=1;
        end
  end
    if frame>1
        peak=max(response(:));
        % reliable sample
        if peak>0.013
            rd_update=1;
        else
            rd_update=0;
        end
        % re-detection the object
        if peak<0.01&&frame>30&&prod(target_sz)>100&&prod(target_sz)<20000&&in_sight==1
            origin_sz=target_sz;
            if success==1
                range=5;
                thred=0.02;
            else
                range=range*1.05;
                thred=thred/1.005;
            end
            Bbs=run_EdgeBox(im,pos,target_sz,model,opts,range);
            if size(Bbs,1)~=0
                [pos,target_sz,success] = re_detection(params,Bbs,pos,target_sz,im,rg_f,features,global_feat_params,use_sz,search_area_scale,thred,featureRatio);
                if success==1
                    redetection=1;
%                     enable=0;
                if prod(target_sz)/prod(origin_sz)>3
                    target_sz=prod(origin_sz)*3/prod(target_sz)*target_sz;
                end
                if prod(target_sz)/prod(origin_sz)<1/3
                    target_sz=prod(origin_sz)*0.33/prod(target_sz)*target_sz;
                end
                search_area = prod(target_sz * search_area_scale);
                [CurrentSvFactor,CurrentArcFactor,base_target_sz,reg_sz,sz,use_sz] = init_size(params,target_sz,search_area,featureRatio);
                [y_0,cos_window] = init_gauss_win(params,base_target_sz,featureRatio,use_sz);
                y=y_0;
                yf=fft2(y);
                [ysf,svFactors,arcFactors,size_model_sz,min_sv_factor,max_sv_factor,min_arc_factor,max_arc_factor] = size_factor(params,target_sz,sz,base_target_sz,im);
                target_sz = base_target_sz * CurrentSvFactor;
                target_sz(1)=floor(target_sz(1)*sqrt(CurrentArcFactor));
                target_sz(2)=floor(target_sz(2)/sqrt(CurrentArcFactor));
                
                end
            
                % Pre-computes the grid that is used for score optimization
               ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
               kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
               sy=circshift(-floor((params.num_sv- 1)/2) : ceil((params.num_sv - 1)/2), [1, -floor((params.num_sv - 1)/2)]);
               sx=circshift(-floor((params.num_arc - 1)/2) : ceil((params.num_arc - 1)/2), [1, -floor((params.num_arc - 1)/2)]);
            else 
                success=0;
            end
            
        end
    end
    
     %save position
     rect_position(loop_frame,:) =[pos([2,1]) - (target_sz([2,1]))/2, target_sz([2,1])];
if success==1
        if frame==1||redetection==1  
            % extract training sample image region
             pixels = get_pixels(im,pos,round(sz*CurrentSvFactor),sz);
             pixels = uint8(gather(pixels));
             x=get_features(pixels,features,params.t_global);
             xf=fft2(bsxfun(@times,x,cos_window));
             % extract training sample for re-detection
             rxf=xf;
        else
            % use detection features
            shift_samp_pos = 2*pi * translation_vec ./(CurrentSvFactor* sz);
            xf = shift_sample(xtf, shift_samp_pos, kx', ky');
            if rd_update==1
                rxf=xf;
            end
        end
        if  frame == 1||redetection==1  
            [~,~,w]=init_regwindow(use_sz,reg_sz,params);
            rw=w;
            g_pre= zeros(size(xf));
            rg_pre= zeros(size(xf));
            mu = 0;
            rmu=0;
            beg_sz=reg_sz;
       else
            mu=params.init_mu;
            rmu=params.init_rmu;
            reg_sz(1)=round(beg_sz(1)*sqrt(CurrentArcFactor));
            reg_sz(2)=round(beg_sz(2)/sqrt(CurrentArcFactor));
            [~,~,w]=init_regwindow(use_sz,reg_sz,params);
            rw=w;
        end
            [g_f,g_pre] = run_training(xf,use_sz,g_pre,params,mu,yf,w);
            if rd_update==1
            [rg_f,rg_pre]=run_training(rxf,use_sz,rg_pre,params,rmu,yf,rw);
            end
     
        %% Update Scale
      if frame==1||redetection==1
        xs = get_size_sample(im, pos, base_target_sz, CurrentSvFactor * svFactors, CurrentArcFactor * arcFactors,size_model_sz);
        xsf = fft2(xs);
     else
        if size_vector(1)~=0||size_vector(2)~=0
        shift_scale_pos = 2*pi * size_vector ./[params.num_sv, params.num_arc];
        xsf = shift_sample(xsf, shift_scale_pos, sx', sy');
        end
     end   
     new_sf_num = ysf.*conj(xsf);
     new_sf_den = sum(xsf .* conj(xsf), 3);
    if frame == 1||redetection==1
            sf_den = new_sf_den;
            sf_num = new_sf_num;
    else
            sf_den = (1 - params.learning_rate_size) * sf_den + params.learning_rate_size * new_sf_den;
            sf_num = (1 - params.learning_rate_size) * sf_num + params.learning_rate_size * new_sf_num;
    end
 end
     time = time + toc();

     %%   visualization
     if params.visualization == 1
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
