function cosine_score = cosine_calculate(im,pos,target_sz,sz,features,global_feat_params,gtt)
        target_pixel=get_pixels(im, pos, round(target_sz), sz);             
        tt=get_features(target_pixel,features,global_feat_params);
        var=(gtt-tt).^2;
        s1=gtt.^2;
        s2=tt.^2;
        cosine_score=100*sum(var(:))/(sum(s1(:))*sum(s2(:)));
end

