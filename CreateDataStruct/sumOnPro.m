%% this function calculate the sum of signal on each promoter

function checStruct_new = sumOnPro(checStruct, promoter_size,GP)
 
    tss = GP.stein_tss13; 
    sample_names = fieldnames(checStruct.norm);

    for s = 1 : length(sample_names)

        curr_sample = char(sample_names(s));
        curr_sample_data = checStruct.norm.(curr_sample);

        for i = 1:length(tss)

            curr_coor = tss(i,2:3);
            curr_chr = tss(i,1);

            if isnan(curr_chr)
                sop(i) = NaN;

            elseif curr_coor(1) < curr_coor(2)   

                sop(i) = nansum(curr_sample_data{curr_chr}(curr_coor(1)- promoter_size: curr_coor(1)));

             elseif curr_coor(1) > curr_coor(2)  %reverse

                sop(i) = nansum(curr_sample_data{curr_chr}(curr_coor(1):curr_coor(1)+promoter_size));
              end

        end
        sop(isnan(sop)) = 0;
        sum_over_promoter.(curr_sample) = sop;

    end

    checStruct.sum_over_promoter = sum_over_promoter;
    checStruct_new = checStruct;

end





