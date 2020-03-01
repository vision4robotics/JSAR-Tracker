
function results = run_JSAR_Re(seq)

%  Initialize path
addpath(genpath('ArvTrack_Re_lib/'));
addpath(genpath('EdgeBox/'));
addpath(genpath('models/'));
addpath(genpath('channels/'));
addpath(genpath('detector/'));
addpath(genpath('images/'));
addpath(genpath('matlab/'));
addpath(genpath('videos/'));
%  HOG feature parameters
hog_params.cell_size = 4;
hog_params.nDim   = 31;
%  ColorName feature parameters
cn_params.nDim  =10;
cn_params.tablename = 'CNnorm';
cn_params.useForGray = false;
cn_params.cell_size = 4;
%  Grayscale feature parameters
grayscale_params.nDim=1;
grayscale_params.colorspace='gray';
grayscale_params.cell_size = 4;

% Global feature parameters 
params.t_features = {
    struct('getFeature',@get_colorspace, 'fparams',grayscale_params),...  
    struct('getFeature',@get_fhog,'fparams',hog_params),...
    struct('getFeature',@get_table_feature, 'fparams',cn_params),...
};
%  Global feature parameters
params.t_global.cell_size = 4;                  % Feature cell size

%   Search region + extended background parameters
params.search_area_shape = 'square';    % the shape of the training/detection window: 'proportional', 'square' or 'fix_padding'
params.search_area_scale =5;           % the size of the training/detection area proportional to the target size
params.min_image_sample_size = 200^2;   % Minimum area of image samples
params.max_image_sample_size = 200^2;   % Maximum area of image samples
%   Gaussian response parameter
params.output_sigma_factor =0.054;		% standard deviation of the desired correlation output (proportional to target)
%   Detection parameters
params.newton_iterations= 5;           % number of Newton's iteration to maximize the detection scores
%  Set files and gt
params.name=seq.name;
params.video_path = seq.video_path;
params.img_files = seq.s_frames;
params.wsize    = [seq.init_rect(1,4), seq.init_rect(1,3)];
params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);
params.s_frames = seq.s_frames;
params.no_fram  = seq.en_frame - seq.st_frame + 1;
params.seq_st_frame = seq.st_frame;
params.seq_en_frame = seq.en_frame;
params.ground_truth=seq.init_rect;

%   ADMM parameters, # of iteration, and lambda- mu and betha are set in
%   the main function.
params.admm_iterations = 4;
params.admm_lambda =1;

params.init_mu=16;
params.init_rmu=16;
params.reg_window_max=1e5;
params.reg_window_min=1e-3;

%  Scale parameters
params.num_sv =13;
params.num_arc=13;
params.learning_rate_size=0.014; 
params.sv_sigma_factor =1.5;
params.arc_sigma_factor=1.7;
params.size_model_factor = 1;
params.sv_step=1.03;
params.arc_step=1.02^2;
params.size_lambda = 0.01;
params.size_model_max_area = 32*16;        
params.hog_scale_cell_size = 4;  
%   Debug and visualization
params.visualization =1;

%   Run the main function
results = tracking_JSAR_Re(params);
