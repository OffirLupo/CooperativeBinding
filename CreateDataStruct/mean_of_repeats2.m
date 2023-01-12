
%% This function averages repeats     
function chec_stuct_med = mean_of_repeats2(chec_data_struct,GP)
    
    %removing last two letters (_digit) to have names vector without
    %repeats number
    curr_samples = fieldnames(chec_data_struct.norm);
    myindices = ~cellfun(@isempty,regexp(curr_samples,'\w*_\d','Match'));

    for i = 1:length(myindices)
        curr_ind =   myindices(i);
        if curr_ind
        curr_name = curr_samples{i} ;
        curr_samples{i} = curr_name(1:end-2) ;
        end
    end
    old_names = fieldnames(chec_data_struct.norm);
    new_names = curr_samples;

    chr_len = GP.chr_len;

    for s = 1:length(new_names)  

        curr_samp = new_names(s);
        my_samples = find(strcmp(new_names,curr_samp));

        temp_med_chr = cell(1,16);

        for curr_chr = 1:16
            temp_array = zeros(length(my_samples),chr_len(curr_chr));
            for j = 1:length(my_samples)          
                temp_array(j,:) = chec_data_struct.norm.(old_names{my_samples(j)}){1,curr_chr};              
            end
            temp_med_chr{1,curr_chr} = nanmedian(temp_array,1);
        end

        chec_stuct_med.norm.(curr_samp{1}) = temp_med_chr ;

    end
end