function bbs = edgeBoxes( I, model, varargin )
% Generate Edge Boxes object proposals in given image(s).
%
% Compute Edge Boxes object proposals as described in:
%  C. Lawrence Zitnick and Piotr Dollár
%  "Edge Boxes: Locating Object Proposals from Edges", ECCV 2014.
% The proposal boxes are fast to compute and give state-of-the-art recall.
% Please cite the above paper if you end up using the code.
%
% The most important params are alpha and beta. The defaults are optimized
% for detecting boxes at intersection over union (IoU) of 0.7. For other
% settings of alpha/beta see the ECCV paper. In general larger alpha/beta
% improve results at higher IoU (but using large alpha can be quite slow).
% minScore/maxBoxes control the number of boxes returned and impact speed.
% Finally, a number of additional params listed below are set to reasonable
% defaults and in most cases should not need to be altered.
%
% We recently introduced an adaptive variant of nms that results in better
% Average Recall (AR) and also better subsequent detection performance.
% This variant is described in Section 5.E of the following paper:
%  Jan Hosang, Rodrigo Benenson, Piotr Dollár, and Bernt Schiele
%  "What makes for effective detection proposals?", arXiv 2015.
% TL;DR: to get top AR for 1000 boxes set alpha=.65, beta=.90, eta=.9996.
%
% For a faster version the proposal code runs at ~10 fps on average use:
%  model.opts.sharpen=0; opts.alpha=.625; opts.minScore=.02;
%
% The code uses the Structured Edge Detector to compute edge strength and
% orientation. See edgesDetect.m for details. Alternatively, the code could
% be altered to use any other edge detector such as Canny.
%
% The input 'I' can either be a single (color) image (or filename) or a
% cell array of images (or filenames). In the first case, the return is a
% set of bbs where each row has the format [x y w h score] and score is the
% confidence of detection. If the input is a cell array, the output is a
% cell array where each element is a set of bbs in the form above (in this
% case a parfor loop is used to speed execution).
%
% USAGE
%  opts = edgeBoxes()
%  bbs = edgeBoxes( I, model, opts )
%
% INPUTS
%  I          - input image(s) of filename(s) of input image(s)
%  model      - Structured Edge model trained with edgesTrain
%  opts       - parameters (struct or name/value pairs)
%   (1) main parameters, see above for details
%   .name           - [] target filename (if specified return is 1) % Note (Bingbin): return is as normal
%   .savename       - [] target filename to save bounding boxes (added by Bingbin)
%   .alpha          - [.65] step size of sliding window search
%   .beta           - [.75] nms threshold for object proposals
%   .eta            - [1.0] adaptation rate for nms threshold (see arXiv15)
%   .minScore       - [.01] min score of boxes to detect
%   .maxBoxes       - [1e4] max number of boxes to detect
%   (2) additional parameters, safe to ignore and leave at default vals
%   .edgeMinMag     - [.1] increase to trade off accuracy for speed
%   .edgeMergeThr   - [.5] increase to trade off accuracy for speed
%   .clusterMinMag  - [.5] increase to trade off accuracy for speed
%   .maxAspectRatio - [3] max aspect ratio of boxes
%   .minBoxArea     - [1000] minimum area of boxes
%   .gamma          - [2] affinity sensitivity, see equation 1 in paper
%   .kappa          - [1.5] scale sensitivity, see equation 3 in paper
%
% OUTPUTS
%  bbs        - [nx5] array containing proposal bbs [x y w h score]
%
% EXAMPLE
%
% See also edgeBoxesDemo, edgesDetect
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar and Larry Zitnick, 2014.
% Licensed under the MSR-LA Full Rights License [see license.txt]

  fid = fopen('time-edgeBoxes.txt', 'w');

  tic;

% get default parameters (unimportant parameters are undocumented)
% Modified: added 'savename'
  dfs={'name','', 'savename','', 'alpha',.65, 'beta',.75, 'eta',1, 'minScore', .01, ...
    'maxBoxes',1e4, 'edgeMinMag',.1, 'edgeMergeThr',.5,'clusterMinMag',.5,...
    'maxAspectRatio',3, 'minBoxArea',1000, 'gamma',2, 'kappa',1.5 };
  o=getPrmDflt(varargin,dfs,1); if(nargin==0), bbs=o; return; end

  fprintf(fid, sprintf('preparation time: %.3f', toc));
  
  tic;

  % run detector possibly over multiple images and optionally save results
  if(~iscell(I))
    bbs=edgeBoxesImg(I,model,o, 0); % Modified by Bingbin: pass i
  else
    n=length(I);
    bbs=cell(n,1);
    imgs=cell(n,1);
    parfor i=1:n
      bbs{i} = edgeBoxesImg(I{i},model,o, i); % Modified by Bingbin: pass i
      imgs{i} = I{i}; % image path
%       if(size(bbs{i}, 1)) % debug
%           disp(sprintf('size of bbs{%d}:', i));
%           disp(size(bbs{i}));
%       end
    end
  end

  fprintf(fid, sprintf('total edgeBoxesImg time: %.3f', toc));

  f=o.name;
  if(~isempty(f) && exist(f,'file'))
    disp(sprintf('edgeBoxes.m: file %s exists.\nWill overwrite.', f));
%    bbs=1; return;
  end
  d=fileparts(f); if(~isempty(d)&&~exist(d,'dir')), mkdir(d); end
  if(~isempty(f)), save(f,'bbs');
    % Modified: return bbs as normal (i.e. not set to 1) & save boxes to mat file
    % bbs=1;
    disp(sprintf('bbs saved at %s', f));
    thres = 0.05;
    % disp(max(bbs(:,5)));

    tic;
    highly_scored = bbs(bbs(:, 5)>thres,:);
    disp(size(highly_scored, 1));
    fprintf(fid, sprintf('get highly_scored time: %.3f', toc));

    tic;
    boxes = horzcat(highly_scored(:,1), highly_scored(:,2), highly_scored(:,1)+highly_scored(:,3), highly_scored(:,2)+highly_scored(:,4), highly_scored(:,5));
    % Modified: show bounding boxes on image
%{
    boxes_cell = cell(1,1); % Mark: unfinished
    boxes_cell{1} = [boxes];
    show_faster(I, boxes_cell);
%}
    % (Noted by Bingbin) It is an error if bounding box coordinates exceed the image boundary: 720 * 1280
    assert(size(nonzeros(boxes(:,[1,3])>1280),1)==0); % for width
    assert(size(nonzeros(boxes(:,[2,4])>720),1)==0); % for height

    top500 = sortrows(boxes, -5);

%    top500 = top500(1:500,:);
    if size(boxes,1)>500
      boxes = top500(1:500,:);
      disp('limited to 500');
    end

%{
    fsave = o.savename;
    disp(fsave);
    if (isempty(fsave))
      fsave = 'edgeBoxes.mat';
    end
%}
%    fsave = 'edgeBox-cat.mat';
    fsave = o.savename;
    save(fsave, 'boxes');

% Save proposal boxes to txt file
    fid = fopen(strcat(fsave, '.txt'), 'w');
    disp(size(boxes,1)); disp(size(boxes,2));
    for row = (1:size(boxes,1))
      fprintf(fid, sprintf('%d,%d,%d,%d,%.5f\n',boxes(row,1), boxes(row,2), boxes(row,3), boxes(row,4), boxes(row,5) ));
    end


%    save('top500-scored.mat', 'top500');
%    disp(sprintf('boxes scored over %d saved at %s', thres, fsave));
  end
end

function bbs = edgeBoxesImg( I, model, o, i)
%  fid = fopen('time-edgeBoxesImg.txt', 'a');
% Generate Edge Boxes object proposals in single image.
  fname = I;
  if(all(ischar(I)))
      % disp(strcat('reading image: ', I));
      I=imread(I);
  end
  model.opts.nms=0;
  [E,O]=edgesDetect(I,model);
  if(0), E=gradientMag(convTri(single(I),4));
    E=E/max(E(:));
  end
% Note: edgesNmsMex.mexa64: compiled from 'private/edgesNmsMex.cpp'
%  nmsTime=tic;
  E=edgesNmsMex(E,O,2,0,1,model.opts.nThreads);
%  fprintf(fid, sprintf('[edgeBoxesImg] [iter %d] time for edgesNmsMex: %.3f', i, toc(nmsTime)));
% each row in bbs = 4 coordinates + score (i.e. 5 elems per row)
% rows are sorted by the score, descending
%  eBMTime=tic;
  bbs=edgeBoxesMex(E,O,o.alpha,o.beta,o.eta,o.minScore,o.maxBoxes,...
    o.edgeMinMag,o.edgeMergeThr,o.clusterMinMag,...
    o.maxAspectRatio,o.minBoxArea,o.gamma,o.kappa);

    if (size(bbs, 1) == 0)
        disp('No boxes proposed');
    end
%  disp('bbs size:');
%  disp(size(bbs));
%  tmp = bbs(1, :);
%  disp('bbs(1} size:');
%  disp(size(tmp));
%  disp('tmp: ');
%  disp(tmp);
%  disp('tmp(1,:):');
%  disp(tmp(1,:));
end
