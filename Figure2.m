%% Figure 2 - Load data 
% data for this analysis is organized in a matlab struct containing the
% normalized genomic tracks of all samples (out files). this struct was
% created by running the function checStruct_paired on all out files, and
% then running the function mean_of_repeats2 to average repeats

addpath(genpath(cd));
load('DataStructs\chec_struct_med.mat');
load('DataStructs\labWTs.mat');
load('DataStructs\MNase_seq.mat');
load('GeneralStructs\tf_motifs.mat');
GP = load('GeneralStructs\general_params_130711.mat');

gene_names = GP.gene_infoR64.name;
for i = 1:6701
    if isempty(gene_names{i})
       gene_names{i} = 'nan';
    end
end

%calculate sum signal on promoter
chec_struct_med = sumOnPro(chec_struct_med,700,GP);
labWTs = sumOnPro(labWTs,700,GP);
%take the signal next to each TSS
MNase_seq = tss_signal(MNase_seq,700,100,GP);
%% load colormaps
cm_blues = cbrewer('seq','Blues',120);
cm_blues(cm_blues<0)  = 0;
cm_blues(cm_blues>1) = 1;
cm_blues = cm_blues(1:100,:);

cm_green = cbrewer('seq','Greens',120);
cm_green(cm_green<0)  = 0;
cm_green(cm_green>1) = 1;
cm_green = cm_green(1:100,:);

cm_reds = cbrewer('seq','Reds',100);
cm_reds(cm_reds<0)  = 0;
cm_reds(cm_reds>1) = 1;
cm_reds = cm_reds(1:90,:);

cm_rb = cbrewer('div','RdBu',100);
cm_rb(cm_rb<0)  = 0;
cm_rb(cm_rb>1) = 1;
cm_rb = cm_rb(6:95,:);


%% Figure 2A,B - change in binding on common vs unique targets
del_list = {'','_dmsn4','_dmsn2','_ddmsn'};
tf_list = {'Yap1','Hsf1','Crz1','Rsf2','Tda9','Skn7','Gis1','Asg1'};

bins_edges = 9:0.5:18;
n_prom = 50;

p = 1;
figure('Position',[2037         528         322         348])
for i = 1:length(tf_list)
    curr_tf =tf_list{i}; 
    [common_promoters, unique_promoters] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);
    promoter_list = [common_promoters ,unique_promoters];
    del_dataMat = zeros(length(del_list),length(promoter_list));
    
    Msn2_mat = log2(chec_struct_med.sum_over_promoter.Msn2(promoter_list));
    
    for d = 1:length(del_list)
        curr_del = curr_tf;
        if d == 1
            curr_data = log2(chec_struct_med.sum_over_promoter.(curr_del));
            del_dataMat(d,:) = curr_data(promoter_list)';   
        end
        if d >1
          curr_del =  [curr_tf ,del_list{d}];
          curr_data = log2(chec_struct_med.sum_over_promoter.(curr_del));
          del_dataMat(d,:) = curr_data(promoter_list)';
        end    
    end

    [~,Idx] = sort(del_dataMat(length(del_list),1:n_prom),'descend');
    del_dataMat(:,1:n_prom) = del_dataMat(:,Idx);
    
    subset_unique = del_dataMat(:,n_prom+1:n_prom*2);
    [~,Idx] = sort(subset_unique(length(del_list),:),'descend');
    del_dataMat(:,n_prom+1:n_prom*2) = subset_unique(:,Idx);
    
    %TFs
    a = subplot(6,length(tf_list)*2,[p,p+1,p+length(tf_list)*2,p+length(tf_list)*2+1]);
    imagesc(del_dataMat)
    set(gca,'ytick',1:length(del_list),'yTickLabel',{'WT','dmsn4','dmsn2','ddmsn'},'FontSize',8,'TickLabelInterpreter','none');
    xlabel('Promoters')
    title(curr_tf)
    set(gcf,'color','w')
    colormap(a,cm_blues)
    caxis([12 16])
    
    %Msn2
    a = subplot(6,length(tf_list)*2,[p+length(tf_list)*4,p+length(tf_list)*4+1]);
    imagesc(Msn2_mat)
    set(gca,'ytick',1,'yTickLabel',{'Msn2'},'FontSize',8,'TickLabelInterpreter','none');
    xlabel('Promoters')
    set(gcf,'color','w')
    colormap(a,cm_blues)
    caxis([12 16])
    
    %Common 
    b = subplot(6,length(tf_list)*2,p+length(tf_list)*6);
    proms_idx = 1:n_prom;
    data_mat = hist(del_dataMat(:,proms_idx)',bins_edges);
    data_mat = smooth_data(data_mat',3);
    joyPlot(data_mat',bins_edges,3,'FaceColor',nanmedian(del_dataMat(:,proms_idx),2),'FaceAlpha',0.95)
    set(gca,'ytick',1:3:12,'YTickLabel',{'ddmsn','dmsn4','dmsn2','WT'})
    caxis([13 14.5]) 
    colormap(b,cm_blues)
    
    %Unique
    c = subplot(6,length(tf_list)*2,p+length(tf_list)*6+1);
    proms_idx = n_prom+1:n_prom*2;
    data_mat = hist(del_dataMat(:,proms_idx)',bins_edges);
    data_mat = smooth_data(data_mat',3); 
    data_mat = smooth_data(data_mat,3);
    joyPlot(data_mat',bins_edges,3,'FaceColor',nanmedian(del_dataMat(:,proms_idx),2),'FaceAlpha',0.95)
    set(gca,'ytick',1:3:12,'YTickLabel',{'ddmsn','dmsn4','dmsn2','WT'})
    colormap(c,cm_blues)

    caxis([13 14]) 
    p = p+2;
end
set(gcf,'color','w')


%% Figure 2B - Plot promoter GSY1
tf_list = {'Msn2','Crz1','Crz1_ddmsn','Rsf2','Rsf2_ddmsn','Yap1','Yap1_ddmsn'};
plot_promoterSlim(chec_struct_med,MNase_seq,GP.gene_table.GSY1,tf_list,700,tf_motifs,GP)

%% Figure 1C top - plot a heatmap of binding around motifs of Msn2, TF in common and TF in unique 
clear all_tables
tf_list = {'Crz1'};
motif_tfs = {'Crz1'};
data_to_plot = {'_ddmsn'};
motif_labels = {'Msn2','Common','Unique'};
n_prom = 50;
window = 100;

all_samples = fieldnames(chec_struct_med.norm);
c = [1,4,2,5,3,6];
figure('Position',[2021         248         507         648])
meanVec = zeros(1,window*2);
p = 1 ;

for t = 1:length(tf_list)   
    curr_tf = tf_list{t};
    motif_tf = motif_tfs{t};
    t_ind = find(startsWith(all_samples,curr_tf)); % tf2
    X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
    Y = vertcat(X{:});
    t_ind = [t_ind(1), t_ind(Y)'];   
    
    %Msn2
    [~,Idx] = sort(chec_struct_med.sum_over_promoter.Msn2,'descend');
    proms = Idx(1:n_prom); 
    motif_table = motif_table_sort('Msn2',{'Msn2'},tf_motifs,chec_struct_med,proms);
    all_tables.M1 = motif_table;
    
    %TF - common
    [proms, ~] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);    
    motif_table = motif_table_sort(motif_tf,{motif_tf},tf_motifs,chec_struct_med,proms);
    all_tables.M2 = motif_table;
    
    %TF - Unique
    [~, proms] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);
    motif_table = motif_table_sort(motif_tf,{motif_tf},tf_motifs,chec_struct_med,proms);
    all_tables.M3 = motif_table;    

    
    all_table_names = fieldnames(all_tables);
    
    for i = 1:length(all_table_names)
        curr_table = all_tables.(all_table_names{i});
        sumMat = zeros(length(t_ind),1);
        
        for d = 1:length(t_ind)
            motif_mat = nan(size(curr_table,1),window*2);
            curr_strain = all_samples{t_ind(d)};
            for m = 1:size(curr_table,1)
                  curr_chr = curr_table.chr(m);
                  curr_loc  = curr_table.loc(m);
                  curr_tss  = curr_table.motif_plusStrand(m);
                  curr_data = chec_struct_med.norm.(curr_strain){curr_chr}(curr_loc-window:curr_loc+window-1);
                if curr_tss == 0
                    curr_data = flip(curr_data);
                end
                    motif_mat(m,:) = curr_data ;                            
            end
            
            meanVec(d,:) = mean(motif_mat);
            motif_mat = smooth_data(motif_mat,3);
            if d == 1
                 [~,Idx] = sort(nansum(motif_mat(:,window-10:window+10),2));
            end

            subplot(2,3,c(p))
            imagesc(motif_mat(Idx,:))
            if i == 1
                 title(curr_strain,'fontsize',12,'Interpreter','none');  
            end
            ylabel(motif_labels{i})
            caxis([5 100])
            colorbar
            p = p+1;
            hold on        
        end         
    end     
end

colormap(cm_blues)
set(gcf,'color','w');

%% Figure 2C bottom - plot a matrix of binding around motifs of Msn2, TF in common and TF in unique 

% general parameters%
tf_list = {'Crz1','Yap1','Tda9','Hsf1','Rsf2','Skn7','Gis1'};
data_to_plot = {'_ddmsn'}; %which background, can add also d7,d11...
motif_labels = {'Msn2','Common','Unique'};
colors_to_plot = {rgb('Black'),rgb('Crimson'),rgb('Orange')};
n_prom = 50;
window = 400;
minLim = 40;
maxLim = 200;

%%%%%%%%%%%

clear all_tables
meanVec = [];
all_samples = fieldnames(chec_struct_med.sum_over_promoter);

%colormap
cm_rb2 = cbrewer('seq','Blues',230);
cm_rb2 = cm_rb2(1:200,:);
cm_rb2(cm_rb2<0)  = 0;
cm_rb2(cm_rb2>1)  = 1;


c = 1;
figure('Position',[2028         204         594         518])
[ha, ~] = tight_subplot(length(tf_list),length(motif_labels),[.01 .03],[.1 .1],[.1 .1]);

 
for t = 1:length(tf_list)   
    curr_tf = tf_list{t};
    t_ind = find(contains(all_samples,curr_tf)); % tf2
    X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
    Y = vertcat(X{:});
    t_ind = [t_ind(1), t_ind(Y)'];   
    
    %Msn2
    [~,Idx] = sort(chec_struct_med.sum_over_promoter.Msn2,'descend');
    proms = Idx(1:n_prom); 
    motif_table = motif_table_sort('Msn2',{'Msn2'},tf_motifs,chec_struct_med,proms);
    all_tables.M1 = motif_table;
    
    %TF - common
    [proms, ~] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);    
    motif_table = motif_table_sort(curr_tf,{curr_tf},tf_motifs,chec_struct_med,proms);
    all_tables.M2 = motif_table;
    
    %TF - Unique
    [~, proms] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);
    motif_table = motif_table_sort(curr_tf,{curr_tf},tf_motifs,chec_struct_med,proms);
    all_tables.M3 = motif_table;    
    
    all_table_names = fieldnames(all_tables);

    for i = 1:length(all_table_names)
        axes(ha(c))

        curr_table = all_tables.(all_table_names{i});
        sumMat = zeros(length(t_ind),1);
        
        for d = 1:length(t_ind)
            motif_mat = zeros(size(curr_table,1),window*2);
            curr_strain = all_samples{t_ind(d)};
            for m = 1:size(curr_table,1)
                  curr_chr = curr_table.chr(m);
                  curr_loc  = curr_table.loc(m);
                  curr_tss  = curr_table.motif_plusStrand(m);
                  curr_data = chec_struct_med.norm.(curr_strain){curr_chr}(curr_loc-window:curr_loc+window-1);

                if curr_tss == 0
                    curr_data = flip(curr_data);
                end
                    motif_mat(m,:) = curr_data ;                            
            end
            
            motif_mat = smooth_data(motif_mat,3);
            meanVec(d,:) = smooth_offir(nanmean(motif_mat(:,:)),10);
            plot(meanVec(d,:),'Color',colors_to_plot{d},'Linewidth',0.2)
            hold on        
            axis tight
            
            if i == 1
              ylabel(curr_tf,'fontsize',12);  
            end
        end
        
         yl = ylim;
         text(1, yl(2)-2,num2str(size(curr_table,1)));
         maxVal = round(yl(2));
         if maxVal > maxLim
             maxVal = maxLim;
         elseif maxVal <= minLim
             maxVal =minLim+1;
         end
         maxVal = maxVal - minLim;
        set(gca,'Color',cm_rb2(maxVal,:))
        
        if t == 1
            title(motif_labels{i})
        end
         
        if t <length(tf_list)
             set(gca,'xtick',[])
        else
             set(gca,'xtick',[1,window,window*2-1],'XTickLabel',[-1*(window),0,window],'XTickLabelRotation',45)
        end
        c = c+1;      
    end  
end

set(gcf,'color','w');

%% Figure 2D - Nuclesome occupancy at full promoter
data_to_plot = {'WT','dmsn4','ddmsn'};
tf = 'Msn2';
colors_to_plot = {rgb('Navy'),rgb('SkyBlue'),rgb('Crimson'),rgb('Orange'),rgb('Green'),rgb('Gold')};
n_proms = 50;
window = 400;
[~,Idx] = sort(chec_struct_med.sum_over_promoter.(tf),'descend');
proms = Idx(1:n_proms);

figure
subplot(2,2,1)

for i = 1:length(data_to_plot)
    curr_data = MNase_seq.signalTss.(data_to_plot{i});
    plot(nanmean(curr_data(proms,:)),'Color',colors_to_plot{i},'Linewidth',2)
    hold on  
end

set(gca,'xtick',[])
axis tight
legend(data_to_plot,'Location','northeast','Interpreter','none')
ylabel('Nuc. Signal')
title(['Avg nucleosome ', tf, ' Promoters'],'FontSize',10); 
set(gcf,'color','w')

subplot('Position',[ 0.1328    0.4946    0.3313    0.0684])
imagesc(mean(MNase_seq.signalTss.Msn2_binding(proms,:)))
colormap(flipud(cm_blues))
set(gca,'xtick',0:100:800,'XTickLabel',[{-1*flip(100:100:700)},{'TSS'},{'+100'}],'XTickLabelRotation',45,'FontSize',12)
set(gca,'ytick',1,'YTickLabel','Msn2')
caxis([0 100])

% nucleosome around motifs
motif_table= motif_table_sort(tf,{tf},tf_motifs,chec_struct_med,proms);
motif_mat = zeros(size(motif_table,1),window*2,length(data_to_plot));
Msn2_mat = zeros(size(motif_table,1),window*2);


for i = 1:size(motif_table,1)
      curr_chr = motif_table.chr(i);
      curr_loc  = motif_table.loc(i);
      curr_tss  = motif_table.tss_right(i);
      for j = 1:length(data_to_plot)          
             curr_data = MNase_seq.norm.(data_to_plot{j}){curr_chr}(curr_loc-window:curr_loc+window-1);
             curr_Msn2 =  chec_struct_med.norm.Msn2{curr_chr}(curr_loc-window:curr_loc+window-1);
            if curr_tss == 0
                curr_data = flip(curr_data);
                curr_Msn2 = flip(curr_Msn2);
            end
            motif_mat(i,:,j) = curr_data ;
            Msn2_mat(i,:) = curr_Msn2 ;
      end
end

subplot(2,2,2)
for i = 1:length(data_to_plot)
    plot(mean(motif_mat(:,:,i)),'Color',colors_to_plot{i},'Linewidth',2)
    hold on
end

axis tight
ylabel('Nuc. Signal')
title(['Avg nucleosome ', tf, ' motifs'],'FontSize',10); 
set(gca,'xtick',[])

subplot('Position',[ 0.5739    0.4875    0.3297    0.0684])
imagesc(mean(Msn2_mat))
colorbar('position',[0.9094    0.4881    0.0211    0.0667]);
colormap(cm_blues)
set(gca,'xtick',0:100:window*2,'XTickLabel',[-1*flip(0:100:window),100:100:window],'XTickLabelRotation',45,'FontSize',12)
set(gca,'ytick',[])

%% Figure 2E,F - Gal11 in ddmsn2,4 scatter, Msn2 in dGa11
figure
a = subplot(1,2,1);
xData = chec_struct_med.sum_over_promoter.Msn2;
yData = chec_struct_med.sum_over_promoter.Msn2_dgal11;
cData = chec_struct_med.sum_over_promoter.Gal11;
scatter(xData,yData,60,cData,'filled','MarkerEdgeColor','k');
hold on
plot(xlim,xlim,'--k')
yl = ylim;
corrVal = round(corr(xData',yData'),2);
text(100,yl(2)-100,['R= ',num2str(corrVal)],'fontsize',12)
caxis([0 50000])
axis square;
set(gca,'FontSize',12)
cb = colorbar;
cb.Label.String = 'Gal11 binding';
xlabel('Msn2'); ylabel('Msn2 dGal11');
colormap(a,cm_blues)

b = subplot(1,2,2);
xData = chec_struct_med.sum_over_promoter.Gal11;
yData = chec_struct_med.sum_over_promoter.Gal11_ddmsn;
cData = chec_struct_med.sum_over_promoter.Msn2;
scatter(xData,yData,60,cData,'filled','MarkerEdgeColor','k');
hold on
plot(xlim,xlim,'--k')
yl = ylim;
corrVal = round(corr(xData',yData'),2);
text(100,yl(2)-100,['R= ',num2str(corrVal)],'fontsize',12)
caxis([0 50000])
axis square; 
set(gcf,'color','w')
set(gca,'FontSize',12)
cb = colorbar;
cb.Label.String = 'Msn2 binding';
xlabel('Gal11'); ylabel('Gal11 ddmsn2,4');
colormap(b,cm_reds)

%% Figure 2G - TF binding on common or unique targets
[~,I] = sort(chec_struct_med.sum_over_promoter.Msn2,'descend');
Msn2_top = I(1:50);

tf_list = {'Rsf2','Tda9','Crz1','Skn7','Yap1','Hsf1'};
del_list = {'','_dGal11','_ddmsn'};
n_prom = 50;

p = 1;
curr_data = [];

figure
med_delta_common = zeros(length(tf_list),2);
med_delta_unique = zeros(length(tf_list),2);

for i = 1:length(tf_list)
    curr_tf =tf_list{i}; 
    
    %take data of curr tf
    wt_data = chec_struct_med.sum_over_promoter.(curr_tf);
    dgal11_data = chec_struct_med.sum_over_promoter.([curr_tf,'_dGal11']);
    ddmsn_data = chec_struct_med.sum_over_promoter.([curr_tf,'_ddmsn']);
    
    %take most effected promoters in the ddmsn 
    [~,effected_in_dd] = sort(wt_data - ddmsn_data,'descend');
    effected_in_dd = effected_in_dd(1:n_prom);      
    common_promoters =  intersect(Msn2_top,effected_in_dd);       
    [~, unique_promoters] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);
   
    common_mat = zeros(length(common_promoters),length(del_list));
    unique_mat = zeros(length(unique_promoters),length(del_list));   
    
    for d = 1:length(del_list)
        curr_del = curr_tf;    
        if d == 1
              curr_data = zscore(chec_struct_med.sum_over_promoter.(curr_del));
              common_mat(:,d) = curr_data(common_promoters)';   
              unique_mat(:,d) = curr_data(unique_promoters)';              
         elseif d >1
              curr_del =  [curr_tf ,del_list{d}];
              curr_data = zscore(chec_struct_med.sum_over_promoter.(curr_del));
              common_mat(:,d) = curr_data(common_promoters)';   
              unique_mat(:,d) = curr_data(unique_promoters)';   
        end    
    end
    
    a = subplot(length(tf_list),2,p);
    [~,Idx1] = sort(common_mat(:,3),'descend');   
    imagesc(common_mat(Idx1,:)');
    caxis([0 6])
    set(gca,'yTickLabel',{'WT','dgal11','ddmsn'});
    ylabel(tf_list{i},'FontSize',12)
    if i ==1
        title('Common')
    elseif i == length(tf_list)
        colorbar('location','southoutside');
    end
    
    med_delta_common(i,1) = nanmedian(common_mat(:,2)) - nanmedian(common_mat(:,1));
    med_delta_common(i,2) = nanmedian(common_mat(:,3)) - nanmedian(common_mat(:,1));
    

    
    b = subplot(length(tf_list),2,p+1);
    med_mat = nanmedian(unique_mat,2);
    [~,Idx2] = sort(med_mat,'descend');  
    imagesc(unique_mat(Idx2,:)');
    caxis([0 10])
    set(gca,'yTickLabel',{'WT','dgal11','ddmsn'});
    ylabel(tf_list{i},'FontSize',12)
    if i ==1
        title('Unique')
    elseif i == length(tf_list)
        colorbar('location','southoutside');
    end   
    med_delta_unique(i,1) = nanmedian(unique_mat(:,2)) - nanmedian(unique_mat(:,1));
    med_delta_unique(i,2) = nanmedian(unique_mat(:,3)) - nanmedian(unique_mat(:,1));

    p = p+2;
    colormap(a,cm_blues)
    colormap(b,cm_green) 
end

set(gcf,'color','w')

figure;
scatter(med_delta_common(:,1),med_delta_common(:,2),90,'filled','markerFaceColor',rgb('LightBlue'),'markerEdgeColor','k')
hold on
gname(tf_list)
a = scatter(med_delta_unique(:,1),med_delta_unique(:,2),90,'filled','markerFaceColor',rgb('LightGreen'),'markerEdgeColor','k')
gname(tf_list,a)
xlim([-6 2]);
ylim([-6 2]);
plot(xlim,xlim,'--k');
legend({'Common','Unique'},'location','northwest')
axis square
hold on
set(gcf,'color','w')
xlabel('dGal11')
ylabel('ddmsn2,4')
set(gca,'fontsize',12)

