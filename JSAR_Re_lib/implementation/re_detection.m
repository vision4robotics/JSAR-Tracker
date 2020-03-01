function [pos,target_sz,success] = re_detection(params,Bbs,pos,target_sz,im,rg_f,features,global_feat_params,use_sz,search_area_scale,thred,featureRatio)
 
    response_max=zeros(size(Bbs,1),1);
    reability=zeros(size(Bbs,1),1);
    for i=1:size(Bbs,1)
        proposal_pos=[Bbs(i,2)+floor(Bbs(i,4)/2),Bbs(i,1)+floor(Bbs(i,3)/2)];
        proposal_size=[Bbs(i,4),Bbs(i,3)];
        search_area = prod(proposal_size * search_area_scale);
        [currentScaleFactor,~,~,~,sz,~] = init_size(params,proposal_size,search_area,featureRatio);
        proposal_pixel=get_pixels(im, proposal_pos, round(sz*currentScaleFactor), 4*[size(rg_f,1),size(rg_f,2)]);             
        pt=get_features(proposal_pixel,features,global_feat_params);
        % construct cosine window
        cos_window = single(hann(size(pt,1)+2)*hann(size(pt,2)+2)');
        cos_window = cos_window(2:end-1,2:end-1);
        ptf=fft2(bsxfun(@times,pt,cos_window));         
        presponsef=permute(sum(bsxfun(@times, conj(rg_f), ptf), 3), [1 2 4 3]);
        presponsef_padded = resizeDFT2(presponsef, use_sz);          
        presponse = ifft2(presponsef_padded, 'symmetric');
        response_max(i)=max(presponse(:));
        reability(i)=response_max(i);
        if response_max(i)<thred
            reability(i)=0;
        end
    end
        if max(reability)>0
            pro_index= reability == max(reability);
            redet_gt=Bbs(pro_index,1:4);
            pos=[redet_gt(2)+floor(redet_gt(4)/2),redet_gt(1)+floor(redet_gt(3)/2)];
            target_sz=[redet_gt(4),redet_gt(3)];
            success=1;
        else
            success=0;
        end
end

