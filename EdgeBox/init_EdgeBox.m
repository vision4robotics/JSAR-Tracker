function [model,opts] = init_EdgeBox()
%% load pre-trained edge detection model and set opts 
model=load('models/forest/modelBsds'); 
model=model.model;
model.opts.multiscale=0; 
model.opts.sharpen=2; 
model.opts.nThreads=4;
%% set up opts for edgeBoxes 
opts = edgeBoxes;
opts.alpha =0.625;     % step size of sliding window search
opts.beta  =0.85;     % nms threshold for object proposals
opts.minScore = 0.01;  % min score of boxes to detect
opts.maxBoxes = 30;  % max number of boxes to detect
end