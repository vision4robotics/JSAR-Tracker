function [xtf,pos,translation_vec,response] = run_detection(im,pos,sz,CurrentSvFactor,features,cos_window,g_f,global_feat_params,use_sz,ky,kx,newton_iterations,featureRatio)
            pixel_template = get_pixels(im, pos, round(sz*CurrentSvFactor), sz);
            xtf = fft2(bsxfun(@times,get_features(pixel_template,features,global_feat_params),cos_window));
            responsef = sum(bsxfun(@times, conj(g_f), xtf), 3);
            % if we undersampled features, we want to interpolate the
            % response so it has the same size as the image patch
            responsef_padded = resizeDFT2(responsef, use_sz);
            % response in the spatial domain
            response = ifft2(responsef_padded, 'symmetric');
            responsef_padded=fft2(response);
            % find maximum peak
            [disp_row, disp_col] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, use_sz);
            % calculate translation
            translation_vec = round([disp_row, disp_col] * featureRatio * CurrentSvFactor);
            % update position
            pos = pos + translation_vec;
end

