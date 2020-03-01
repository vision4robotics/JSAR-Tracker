function results = run_JSAR(seq)

addpath(genpath('ArvTrack_lib/'));
%   HOG feature parameters
hog_params.nDim   = 31;
params.video_path = seq.video_path;
%   Grayscale feature parameters
grayscale_params.colorspace='gray';
grayscale_params.nDim = 1;

cn_params.tablename = 'CNnorm';
cn_params.useForGray = false;
cn_params.cell_size = 4;
cn_params.nDim = 10;

%   Global feature parameters 
params.t_features = {
    struct('getFeature',@get_colorspace, 'fparams',grayscale_params),...  
    struct('getFeature',@get_fhog,'fparams',hog_params),...
    struct('getFeature',@get_table_feature, 'fparams',cn_params),...
 };
params.t_global.cell_size = 4;                  % Feature cell size

%   Search region + extended background parameters
params.search_area_shape = 'square';    % the shape of the training/detection window: 'proportional', 'square' or 'fix_padding'
params.search_area_scale = 5;           % the size of the training/detection area proportional to the target size
params.filter_max_area   = 50^2;        % the size of the training/detection area in feature grid cells

%   Learning parameters
params.output_sigma_factor =1/20;		% standard deviation of the desired correlation output (proportional to target)
%   Detection parameters
params.newton_iterations  = 5;           % number of Newton's iteration to maximize the detection scores
			
%   Scale parameters
params.num_sv =13;
params.num_arc=13;
params.learning_rate_size=0.014; 
params.sv_sigma_factor =1.5;
params.arc_sigma_factor=1.5;
params.size_model_factor = 1;
params.sv_step=1.03;
params.arc_step=1.02^2;
params.size_lambda = 0.01;
params.size_model_max_area = 32*16;        
params.hog_box_cell_size = 4;
params.min_image_sample_size = 200^2;   % Minimum area of image samples
params.max_image_sample_size = 200^2;   % Maximum area of image samples
params.init_mu=16;
%   size, position, frames initialization
params.wsize    = [seq.init_rect(1,4), seq.init_rect(1,3)];
params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);
params.s_frames = seq.s_frames;
params.no_fram  = seq.en_frame - seq.st_frame + 1;
params.seq_st_frame = seq.st_frame;
params.seq_en_frame = seq.en_frame;

%   ADMM parameters, # of iteration, and lambda- mu and betha are set in
%   the main function.
params.admm_iterations = 4;
params.admm_lambda =1e-2;
params.reg_window_max=1e5;
params.reg_window_min=1e-3;
%   Debug and visualization
params.visualization = 1;
%   Run the main function
results = tracking_JSAR(params);
