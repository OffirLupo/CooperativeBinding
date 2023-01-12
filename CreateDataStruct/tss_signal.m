
function MychecStruct = tss_signal(MychecStruct,beforeTss,afterTss,GP)
% This function plots the signal around the TSS for each gene in each
% sample of chec struct

%create signal vector for each gene
samples_names = fieldnames(MychecStruct.norm);

   tss = GP.stein_tss13 ; 
   
    for j = 1:length(samples_names)
        disp(samples_names{j})
        tmp = nan(6701,beforeTss+afterTss);
        curr_sample = MychecStruct.norm.(samples_names{j});
        
       
            for i = 1:6701
                curr_tss = tss(i,:);

                    if curr_tss(2)<curr_tss(3)

                        size_corr = curr_tss(2) -beforeTss+1;

                        tmp(i,:) = curr_sample{curr_tss(1)}(size_corr:curr_tss(2)+afterTss);

                    elseif curr_tss(2)>curr_tss(3)

                       size_corr = curr_tss(2) +beforeTss-1;

                       tmp(i,:)= flip(curr_sample{curr_tss(1)}(curr_tss(2)-afterTss:size_corr)); 

                    elseif isnan(curr_tss(1))

                        tmp(i,:)= zeros(1,beforeTss+afterTss);
                        
                    end
                    
            end 
                     
      tss_struct.(samples_names{j}) = tmp;
      clear tmp curr_sample
      MychecStruct.signalTss = tss_struct;
    
    end
end

