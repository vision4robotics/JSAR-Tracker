function Bbs = run_EdgeBox(im,pos,target_sz,model,opts,range)

%% detect Edge Box bounding box proposals 
target_sz = max(target_sz, 2);
num_pixels=prod(target_sz)*range^2;
if num_pixels>size(im,1)*size(im,2)
    ys=1:size(im,1);
    xs=1:size(im,2);
else
    xs = floor(pos(2)) + (1:round(sqrt(num_pixels))) - floor(1/2*round(sqrt(num_pixels)));
    ys = floor(pos(1)) + (1:round(sqrt(num_pixels))) - floor(1/2*round(sqrt(num_pixels)));
    % check for out-of-bounds coordinates, and set them to the values at
    % the borders
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(im,2)) = size(im,2);
    ys(ys > size(im,1)) = size(im,1);
end
    I=im(ys,xs,:);
% imshow(I);
% tic, 
    bbs=edgeBoxes(I,model,opts); 
% %toc
% gt=bbs;
% gt(:,5)=0; 
% [gtRes,dtRes]=bbGt('evalRes',double(gt),double(bbs),.1);
% figure(1);
% bbGt('showRes',I,gtRes,dtRes(dtRes(:,6)==1,:));
fclose all;
Bbs=zeros(size(bbs,1),5);
for i=1:size(Bbs,1)
    Bbs(i,1:2)=bbs(i,1:2)+[xs(1),ys(1)];
    Bbs(i,3:4)=bbs(i,3:4);
    Bbs(i,5)=bbs(i,5);
    if Bbs(i,1)<1
        Bbs(i,1)=1;
    end
     if Bbs(i,1)>size(im,2)
        Bbs(i,1)=size(im,2);
     end
     if Bbs(i,2)<1
        Bbs(i,2)=1;
     end
     if Bbs(i,2)>size(im,1)
        Bbs(i,2)=size(im,1);
     end
end
end