% find the top motifs in promoters for a list of TFs

function [motif_table] = motif_table_sort(main_tf,tf_list,tf_motifs,chec_data,genes)

    GP = load('GeneralStructs\general_params_130711.mat');
    load('GeneralStructs\promoterLengthsGP.mat');
    cerGenome = fastaread('GeneralStructs\S288C_reference_sequence_R64-1-1_20110203.fsa');
    p_length = 700;
    tss = GP.stein_tss13; 
    
%     tres =0;
%     window_size = 50;

    %% Get motifs of main TF
    curr_motif = tf_motifs.(main_tf);
    m = 1;
    tss = tss(genes,:);
     for i = 1:length(tss)

        curr_coor = tss(i,2:3);
        curr_chr = tss(i,1);
        %p_length = promoterLengthsGP(i);
        if isnan(curr_chr)
           continue
        
        elseif curr_coor(1) < curr_coor(2)   
            promoter_seq = cerGenome(curr_chr).Sequence(curr_coor(1)- p_length: curr_coor(1));

            motif_loc_F = regexp(promoter_seq,curr_motif); % + strand;
            motif_loc_R = regexp(seqrcomplement(promoter_seq),curr_motif); % - strand
            motif_locs = [motif_loc_F, p_length - motif_loc_R];

            motif_strands = [ones(1,length(motif_loc_F)),zeros(1,length(motif_loc_R))];
            if isempty(motif_locs)
                continue
            else
                for j = 1:length(motif_locs)                   
                    tss_right(m) = 1;
                    motif_stand(m) = motif_strands(j);
                    gene_list(m) = genes(i);
                    tf_motif{m} = main_tf;
                    mot_loc(m) = curr_coor(1) - (p_length-motif_locs(j));
                    motif_seq{m} = curr_motif;
                    mot_chr(m) = curr_chr;
                    gene_name{m} = GP.gene_infoR64.name{genes(i)};
                    m = m+1;
                end
            end

            
         elseif curr_coor(1) > curr_coor(2)  %reverse
            promoter_seq = cerGenome(curr_chr).Sequence(curr_coor(1):curr_coor(1)+p_length);

            motif_loc_F = regexp(promoter_seq,curr_motif); % + strand;
            motif_loc_R = regexp(seqrcomplement(promoter_seq),curr_motif); % - strand
            motif_locs = [motif_loc_F, p_length - motif_loc_R];
            
            motif_strands = [ones(1,length(motif_loc_F)),zeros(1,length(motif_loc_R))];

            if isempty(motif_locs)
                continue
            else
                for j = 1:length(motif_locs)
                    tss_right(m) = 0;
                    motif_stand(m) = motif_strands(j);
                    gene_list(m) = genes(i);
                    tf_motif{m} = main_tf;
                    mot_chr(m) = curr_chr;
                    mot_loc(m) = curr_coor(1)+motif_locs(j);
                    motif_seq{m} = curr_motif;
                    gene_name{m} = GP.gene_infoR64.name{genes(i)};
                    m = m+1;
                end
            end

        end
     end

motif_table = table(tf_motif',motif_seq',gene_list',gene_name',mot_chr',mot_loc',tss_right',motif_stand','VariableNames',...
                    {'tf','seq','gene_idx','gene_name','chr','loc','tss_right','motif_plusStrand'});
    
    %% sort motifs dynamicly by treshold
%     motif_signal = zeros(height(motif_table),length(tf_list));
%     
%     for t = 1:length(tf_list)
%         curr_tf = tf_list{t};
%         curr_data = chec_data.norm.(curr_tf);
%         for i = 1:height(motif_table)
%             curr_chr = motif_table.chr(i);
%             curr_loc = motif_table.loc(i); 
%             motif_signal(i,t) = nansum(curr_data{curr_chr}(curr_loc-window_size:curr_loc+window_size-1));
%         end
%     end
%     motif_signal = log2(motif_signal);%
%     %motif_signal = zscore(motif_signal);
%     
% %     for t = 1:length(tf_list)
% %         currCol = motif_signal(:,t);
% %         curr_idx = find(currCol > 12);
% %         %if t == 1
% %             %currCol = motif_signal(curr_idx,t+1);
% %           %  [~,idx] = sort(currCol,'descend'); 
% %         %else
% %         [~,idx] = sort(currCol(curr_idx),'descend');
% %        %end
% %         tf_sorted_ind{t} = curr_idx(idx);
% %         motif_signal( tf_sorted_ind{t} ,:) = 0;
% %         
% %     end
% %    num_pass_tresh = length(vertcat(tf_sorted_ind{:}));
% %    tf_sorted_ind{length(tf_list)+1} = find(~sum(motif_signal,2) == 0);
% %    all_motifs_sorted = vertcat(tf_sorted_ind{:})    ; 
% % 
% %     [~,all_motifs_sorted] = sort(nanmean(motif_signal,2),'descend');
%     motif_table = motif_table(motif_signal>tres,:);
 
end
 

