%% Figure 1 - Load data 
% data for this analysis is organized in a matlab struct containing the
% normalized genomic tracks of all samples (out files). this struct was
% created by running the function checStruct_paired on all out files, and
% then running the function mean_of_repeats2 to average repeats

addpath(genpath(cd));
load('DataStructs\chec_struct_med.mat');
load('DataStructs\MNase_seq.mat');
load('GeneralStructs\tf_motifs.mat');
GP = load('GeneralStructs\general_params_130711.mat');

gene_names = GP.gene_infoR64.name;
for i = 1:6701
    if isempty(gene_names{i})
       gene_names{i} = 'nan';
    end
end

chec_struct_med = sumOnPro(chec_struct_med,700,GP);

%% Figure 4A - scatter plots
del_list_msn2 =  {'Msn2_d4' ,'Msn2_d7A', 'Msn2_d11','Msn2_d14'};

figure
for i = 1:length(del_list_msn2)   

    subplot(2,length(del_list_msn2)/2,i)
    xData = chec_struct_med.sum_over_promoter.Msn2;
    yData = chec_struct_med.sum_over_promoter.(del_list_msn2{i}) ;
    scatter(xData,yData,20,'filled','MarkerFaceColor',rgb('SkyBlue'),'markerEdgeColor','k')
    axis square
    xlabel('Msn2','FontSize',12)
    ylim([0 200000])
    yl = ylim;
    text(100,yl(2)-100,['R= ', num2str(round(corr(xData',yData'),2))])
    ylabel(del_list_msn2{i},'Interpreter','none','FontSize',12)
    hold on
    plot(xlim,xlim,'--k','LineWidth',1)
end
gname(gene_names)

set(gcf,'color','w')

%% Figure 4B - Plot example promoters
del_list_msn2 =  {'Msn2','Msn2_d4','Msn2_d7A','Msn2_d7B','Msn2_d11','Msn2_d14'};
plot_promoterSlim(chec_struct_med,MNase_seq,GP.gene_table.UGP1,del_list_msn2,700,tf_motifs,GP)
%% Figure 4C - correlation of Msn2 WT to Msn2 in deletion strains
all_samples = fieldnames(chec_struct_med.sum_over_promoter);

del_list_msn2 =  {'Msn2_d4' ,'Msn2_d7A' ,'Msn2_d7B', 'Msn2_d11','Msn2_d14'};
corrVecMsn2 = zeros(1,length(del_list_msn2));

for i = 1:length(del_list_msn2)   
    corrVecMsn2(i) = corr(chec_struct_med.sum_over_promoter.Msn2',chec_struct_med.sum_over_promoter.(del_list_msn2{i})','type','pearson');    
end


figure('Position',[2350         422         221         381]);
scatter(1:length(del_list_msn2),corrVecMsn2,60,'filled','MarkerFaceColor',rgb('Crimson'),'MarkerEdgeColor','k')
ylim([0.4 1])
xlim([0.5 length(del_list_msn2)+0.5])
grid on
set(gca,'xtick',1:length(del_list_msn2),'XTickLabel',del_list_msn2,'TickLabelInterpreter','none','XTickLabelRotation',45,'FontSize',10)
ylabel('Correlation with WT')
set(gcf,'color','w')


%% plot a matrix of binding around motifs of other TFs 

%%%%%%%%
tf_list = {'Msn2'};
data_to_plot = {'_d14'};
colors_to_plot = {rgb('Black'),rgb('Salmon')};
motif_tfs = {'Msn2','Msn4','Gis1','Rsf2','Tda9','Adr1','Yap1','Crz1','Skn7','Hsf1','Sko1'};

%%%%%%%%
n_prom = 50;
window = 300;
figure('Position',[2383          67         210         908])

[ha, pos] = tight_subplot(length(motif_tfs),length(tf_list),[.005 .005],[.1 .1],[.1 .1]);
meanVec = zeros(1,window*2);

p =1;
for t = 1:length(motif_tfs)   
    
    tf_forMotifs = motif_tfs{t};
    [proms, ~] = diff_proms(chec_struct_med,'Msn2',tf_forMotifs,n_prom,100);
    motif_table = motif_table_sort(tf_forMotifs,{tf_forMotifs},tf_motifs,chec_struct_med,proms);
    
    for i = 1:length(tf_list)
        curr_tf = tf_list{i};
        wt_ind = find(strcmp(all_samples,curr_tf)); 
        t_ind = find(contains(all_samples,curr_tf)); 
        X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
        Y = vertcat(X{:});
        t_ind =  [wt_ind,t_ind(Y)'];
        
        axes(ha(p))
        for d = 1:length(t_ind)
            motif_mat = zeros(size(motif_table,1),window*2);
            curr_strain = all_samples{t_ind(d)};
            for m = 1:size(motif_table,1)
                  curr_chr = motif_table.chr(m);
                  curr_loc  = motif_table.loc(m);
                  curr_tss  = motif_table.tss_right(m);
                  curr_data = chec_struct_med.norm.(curr_strain){curr_chr}(curr_loc-window:curr_loc+window-1);
                if curr_tss == 0
                    curr_data = flip(curr_data);
                end
                    motif_mat(m,:) = curr_data ;                            
            end
            
            motif_mat = smooth_data(motif_mat,3);          
            meanVec(d,:) = smooth_offir(nanmean(motif_mat),10);
            plot(meanVec(d,:),'Color',colors_to_plot{d},'Linewidth',0.25)
            hold on
            ylim([0 210])
        end
            
        yl = ylim;
        xl = xlim;
        plot([window window],yl,'--k','LineWidth',0.1);
        mean_signal = nanmean(meanVec(:,[window-60:window+60]),2);
                  
        text(xl(2)-50, yl(2)-5,num2str(size(motif_table,1)),'fontsize',8)  
        text(xl(1)+10, yl(2)-5,tf_forMotifs,'fontsize',8)  
         
         p = p+1;

          if t == 1
            title(curr_tf)
          end
                    
           if t <length(motif_tfs)
             set(gca,'xtick',[])
           else
                 set(gca,'xtick',0:100:window*2,'XTickLabel',[-1*flip(0:100:window),100:100:window],'XTickLabelRotation',45)
           end           
    end 
end
set(gcf,'color','w');



%% Figure 4SA - effect of VMA5 promoter

figure;

xData = chec_struct_med.sum_over_promoter.Msn2;
yData = chec_struct_med.sum_over_promoter.Msn2_VMA5 ;
scatter(xData,yData,20,'filled','MarkerFaceColor',rgb('lightBlue'))
axis square
xlabel('Native promoter')
ylabel('VMA5 promoter','FontSize',12)
hold on
yl = ylim;
plot(xlim,xlim,'--k','LineWidth',1)
set(gca,'fontsize',12);
set(gcf,'color','w')
text(100, yl(2),['R = ',num2str(round(corr(xData',yData'),2))]);

hold on
scatter(xData(GP.gene_table.UGP1),yData(GP.gene_table.UGP1),20,'filled','MarkerFaceColor',rgb('red'))
scatter(xData(GP.gene_table.TPS2),yData(GP.gene_table.TPS2),20,'filled','MarkerFaceColor',rgb('red'))

 

 