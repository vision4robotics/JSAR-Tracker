function out = get_size_sample(im, pos, base_target_sz, sv_factors, arc_factors, size_model_sz,window,frame)

num_sv= length(sv_factors);
num_arc=length(arc_factors);

for s = 1:num_sv
    for a=1:num_arc
        patch_sz = floor(base_target_sz * sv_factors(s));
        patch_sz(1) = floor(patch_sz(1)*sqrt(arc_factors(a)));
        patch_sz(2) =  floor(patch_sz(2)/sqrt(arc_factors(a)));
        %make sure the size is not to small
        patch_sz = max(patch_sz, 2);
 
        xs = floor(pos(2)) + (1:patch_sz(2)) - floor(patch_sz(2)/2);
        ys = floor(pos(1)) + (1:patch_sz(1)) - floor(patch_sz(1)/2);
    
        %check for out-of-bounds coordinates, and set them to the values at
        %the borders
        xs(xs < 1) = 1;
        ys(ys < 1) = 1;
        xs(xs > size(im,2)) = size(im,2);
        ys(ys > size(im,1)) = size(im,1);
    
        %extract image
        im_patch = im(ys, xs, :);
%         if frame==5
%             savedir='H:\IROS\ORIGIN\';
%             imwrite(im_patch,[savedir, num2str(s),'_',num2str(a),'.png']);
%             hold off;
%         end
        % resize image to model size
%         im_patch_resized = mexResize(im_patch, size_model_sz, 'auto');
        % extract scale features

%         im_patch_resized = mexResize(im_patch_resized,size_model_sz/hog_box_cell_size,'auto');
%         temp(:,:,1)=single(rgb2gray(im_patch_resized))/255-0.5;
%         temp(:,:,2:11)=fcn(im_patch_resized, 'cn', w2c);
        % resize image to model size
        im_patch_resized = mexResize(im_patch, size_model_sz, 'auto');
    
        % extract scale features
        temp_hog = fhog(single(im_patch_resized),4);
        temp = temp_hog(:,:,1:31);
    
        if a==1&&s==1
            out= zeros(num_sv, num_arc,numel(temp),'single');
        end
%         
%         if frame==5
%             fig_temp=sum(temp,3);
%             savedir='H:\IROS\FEATURE\';
%             colormap(parula);
%             han=figure(2);
%             set(han,'visible','off');  
%             surf(fig_temp);
%             shading interp;
%             axis ij;
%             axis off;
%             view([0,90]);
%             saveas(gcf,[savedir,num2str(s),'_',num2str(a),'.png']);
%             hold off;
%         end
%             
        
        
%         out(s,a,:) = temp(:)*window(s,a);
        out(s,a,:) = temp(:);
    end
end
