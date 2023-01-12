%% this function plots the signal on a promoter

function [] = plot_promoterSlim(checDataStruct,MNaseSeq,gene_idx,TFs,promSize,motifs,GP)

    cerGenome = fastaread('GeneralStructs\S288C_reference_sequence_R64-1-1_20110203.fsa');    
    motifs = RenameField(motifs,fieldnames(motifs),upper(fieldnames(motifs)));
    all_motif_fields = fieldnames(motifs);
    
    %colormap for motifs
    cm  = cbrewer('qual','Set2',20);
    cm(cm>1) = 1;
    %colormap for relative signal
    cm_sumProm = cbrewer('seq','Blues',120);
    cm_sumProm(cm_sumProm<0)  = 0;
    cm_sumProm(cm_sumProm>1) = 1;
    cm_sumProm = cm_sumProm(1:100,:);
    p = 1;
    %%
    %find coordinates of current promoter
    curr_cord = GP.stein_tss13(gene_idx,:);
    gene = GP.gene_infoR64.name(gene_idx);
    if curr_cord(2) < curr_cord(3) %forward 
        curr_gene_seq.cer = cerGenome(curr_cord(1)).Sequence(curr_cord(2)-promSize:curr_cord(2)+100);     
    elseif curr_cord(2) > curr_cord(3) %reverse
        curr_gene_seq.cer = seqrcomplement(cerGenome(curr_cord(1)).Sequence(curr_cord(2)-100:curr_cord(2)+promSize-1));

    end
       
    for i= 1:length(TFs)
        tempStruct1.(TFs{i})  =  checDataStruct.norm.(TFs{i});    
    end
    
    curr_sample_names = fieldnames(tempStruct1);
    figure('Position',[        2660         249         242         333]);
    [ha, ~] = tight_subplot(length(TFs)+1,1,[.01 .03],[.1 .1],[.1 .1]);

    for i = 1:length(TFs)

        curr_TF = curr_sample_names{i};
        
        if curr_cord(2) < curr_cord(3)
            curr_data = tempStruct1.(curr_TF){curr_cord(1)}(curr_cord(2)-promSize:curr_cord(2)+100); %ChEC signal
 %           curr_nuc_dd = MNaseSeq.norm.ddmsn{curr_cord(1)}(curr_cord(2)-promSize:curr_cord(2)+100); %nucleosomes
            curr_nuc_WT = MNaseSeq.norm.WT{curr_cord(1)}(curr_cord(2)-promSize:curr_cord(2)+100); %nucleosomes
        elseif curr_cord(2) > curr_cord(3)
            curr_data = flip(tempStruct1.(curr_TF){curr_cord(1)}(curr_cord(2)-100:curr_cord(2)+promSize-1));%ChEC signal
%            curr_nuc_dd = flip(MNaseSeq.norm.ddmsn{curr_cord(1)}(curr_cord(2)-100:curr_cord(2)+promSize-1));%nucleosomes           
            curr_nuc_WT = flip(MNaseSeq.norm.WT{curr_cord(1)}(curr_cord(2)-100:curr_cord(2)+promSize-1));%nucleosomes
        end
        
        curr_sumProm = checDataStruct.sum_over_promoter.(curr_TF);
        [~,idx] = sort(curr_sumProm,'descend');
        norm_to_max = curr_sumProm./nanmedian(curr_sumProm(idx(1:5)));
        curr_col = ceil(norm_to_max(gene_idx)*100);
      
        axes(ha(i))
        hold on
        if size(curr_data,2) == 1
            curr_data = curr_data';
        end
        a = area(smooth_offir(curr_data,1),'FaceColor','k');
        ylim([0 1800])

        yl = ylim;
%         if yl(2) < 1500
%             ylim([0 1500])
%         end
        text(10,yl(2)-100,TFs{i},'FontSize',10,'Interpreter','none');
        set(gca,'ytick',[0,floor(yl(2))],'YTickLabel',[0,floor(yl(2))])
        if curr_col == 0
            
        elseif curr_col > 100
            curr_col = 100;
        end
       set(gca,'Color',cm_sumProm(curr_col,:))

        curr_tf = strsplit(TFs{i},'_');
        curr_tf = upper(curr_tf{1});
         if i ==1
            title(upper(gene) ,'Interpreter','none','fontWeight','Bold')
         end
        
         %draw motif
         if ismember(curr_tf,all_motif_fields)
             text(promSize,yl(2),motifs.(curr_tf),'fontsize',8);
            m = regexp(curr_gene_seq.cer,motifs.(curr_tf));
            if ~isempty(m)
                name_for_leg{i} = [regexprep(curr_TF,'_',' ' ),' - ' ,motifs.(curr_tf)];
               for j = m
                   scatter(j,0,20,'MarkerFaceColor',cm(p,:),'MarkerEdgeColor','k');
               end
            end
             
            m = regexp(seqrcomplement(curr_gene_seq.cer),motifs.(curr_tf));
            if ~isempty(m)
               for j = m
                   scatter(length(curr_gene_seq.cer)-j,-0,20,'MarkerFaceColor',cm(p,:),'MarkerEdgeColor','k');     
               end
            end
         end
        plot([promSize promSize],yl,'--k','LineWidth',2)
        box on
        %p = p+1;
        
    end
    %plot nucleosome occupancy
       axes(ha(i+1))
       area(smooth_offir(curr_nuc_WT,10),'FaceAlpha',0.5 ,'FaceColor',rgb('grey'));
       hold on
       %area(smooth_offir(curr_nuc_dd,10),'FaceAlpha',0.5 ,'FaceColor',rgb('white'));
       linkaxesInFigure('x')
       set(gcf,'color','w');
       set(gca,'xtick',[1:100:promSize+100],'xticklabel',-1*promSize:100:promSize+100,'FontSize',8,'XTickLabelRotation',45)
       axis tight
       xlabel('Distance from TSS (bp)')
        set(gcf,'color','w');
end
