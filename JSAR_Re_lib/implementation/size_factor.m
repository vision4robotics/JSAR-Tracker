function [ysf,svFactors,arcFactors,size_model_sz,min_sv_factor,max_sv_factor,min_arc_factor,max_arc_factor] = size_factor(params,target_sz,sz,base_target_sz,im)
%% SCALE ADAPTATION INITIALIZATION
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
min_arc_factor =1.5*1/params.search_area_scale^2*(target_sz(1)/target_sz(2));
max_arc_factor =params.search_area_scale^2/(1.5*target_sz(1)/target_sz(2));

end

