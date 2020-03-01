function [xs,CurrentSvFactor,CurrentArcFactor,size_vector]  = search_size(sf_num,sf_den,im,pos,base_target_sz,CurrentSvFactor,CurrentArcFactor,svFactors,arcFactors,size_model_sz,min_sv_factor,max_sv_factor,min_arc_factor,max_arc_factor,params)
            % update position
 
            xs = get_size_sample(im, pos, base_target_sz, CurrentSvFactor * svFactors, CurrentArcFactor*arcFactors,size_model_sz);
            xsf = fft2(xs);
            scale_response = ifft2((sum(sf_num .* xsf, 3) ./ (sf_den + params.size_lambda)),'symmetric'); 
       
            % find the maximum scale response
            [recovered_sv,recovered_arc] = find(scale_response == max(scale_response(:)), 1);
            size_vector=[recovered_sv-floor(size(scale_response,1)/2)-1,recovered_arc-floor(size(scale_response,2)/2)-1];
            % update the scale
            CurrentSvFactor = CurrentSvFactor * svFactors(recovered_sv);
            if CurrentSvFactor < min_sv_factor
                CurrentSvFactor = min_sv_factor;
            elseif CurrentSvFactor > max_sv_factor
                CurrentSvFactor = max_sv_factor;
            end            
            CurrentArcFactor=CurrentArcFactor*arcFactors(recovered_arc);
            if CurrentArcFactor < min_arc_factor
                CurrentArcFactor = min_arc_factor;
            elseif CurrentArcFactor > max_arc_factor
                CurrentArcFactor = max_arc_factor;
            end    

end

