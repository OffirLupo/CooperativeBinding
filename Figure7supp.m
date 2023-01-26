%% Figure 7 - Load data 
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

% calculate sum signal on promoter
chec_struct_med = sumOnPro(chec_struct_med,700,GP);

% load colormaps
cm_YlGnBu = cbrewer('seq','YlGnBu',120);
cm_YlGnBu(cm_YlGnBu<0)  = 0;
cm_YlGnBu(cm_YlGnBu>1) = 1;
cm_YlGnBu = cm_YlGnBu(1:100,:);

cm_rb = cbrewer('div','RdBu',100);
cm_rb(cm_rb<0)  = 0;
cm_rb(cm_rb>1) = 1;
cm_rb = cm_rb(6:95,:);
cm_rb(37:54,:) = repmat(cm_rb(46,:),18,1);

cm_green = cbrewer('seq','Greens',120);
cm_green(cm_green<0)  = 0;
cm_green(cm_green>1) = 1;
cm_green = cm_green(1:100,:);

cm_blues = cbrewer('seq','Blues',120);
cm_blues(cm_blues<0)  = 0;
cm_blues(cm_blues>1) = 1;
cm_blues = cm_blues(1:100,:);

% Create a sumProm matrix of all samples
all_samples = fieldnames(chec_struct_med.sum_over_promoter);
all_combind_mat = zeros(6701,length(all_samples));
all_samples(strcmp(all_samples,'Msn2DBDonly_d11')) = [];

for i =1:length(all_samples)   
    all_combind_mat(:,i) = chec_struct_med.sum_over_promoter.(all_samples{i})';
end

%% Figure S7A - Msn2DBD in d11
tf_list1 = {'Msn2','Msn2_d11','Msn2DBDonly','Msn2DBDonly_d11'};
    
all_combined_mat = zeros(6701,length(tf_list1));

for i =1:length(tf_list1)
    all_combined_mat(:,i) = chec_struct_med.sum_over_promoter.(tf_list1{i});
end

figure
subplot(1,5,1)
imagesc(corr(all_combined_mat,'rows','pairwise'))
axis square
set(gcf,'color','w')
cb = colorbar;
cb.Label.String = 'Pearson r';
title('Promoter preference');
set(gca,'xtick',1:length(tf_list1),'XTickLabel',tf_list1,'XTickLabelRotation',45,'FontSize',12,'TickLabelInterpreter','none')
set(gca,'ytick',1:length(tf_list1),'yTickLabel',tf_list1)

xTF_toPlot = {'Msn2','Msn2_d11'};
yTF_toPlot = {'Msn2DBDonly','Msn2DBDonly_d11'};

for i = 1:length(xTF_toPlot)
    
    a = subplot(1,5,i+1);
    xData = chec_struct_med.sum_over_promoter.(xTF_toPlot{i});
    yData = chec_struct_med.sum_over_promoter.(yTF_toPlot{i});
    scatter(xData,yData,30,'filled')
    corrVal = round(corr(xData',yData'),2);
    yl = ylim;
    xl = xlim;
    text(xl(2)-10000,yl(2)-100,['R= ',num2str(corrVal)],'fontsize',12)
    axis square
    xlabel(xTF_toPlot{i})
    ylabel(yTF_toPlot{i})
    set(gca,'FontSize',12)
   caxis([0 50000])
   set(gca,'fontsize',12)

end

set(gcf,'color','w')
colormap(cm_YlGnBu)
%% Figure S7B - scatter plots effect of nonDBD removal
tf_list = {'Msn2','Yap1','Crz1','Skn7','Gis1','Sko1','Sok2','Phd1','Ixr1','Dot6','Mig1'};

figure('position',[1941         316        1764         517])
for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    subplot(2,6,i)
    xData = chec_struct_med.sum_over_promoter.(curr_tf);
    yData = chec_struct_med.sum_over_promoter.([curr_tf,'DBDonly']);
    scatter(xData,yData,50,'filled','MarkerFaceColor',rgb('lightBlue'))
    corrVal = round(corr(xData',yData'),2);
    yl = ylim;
    xl = xlim;
    text(xl(2)-10000,yl(2)-100,['R= ',num2str(corrVal)],'fontsize',12)
    axis square
    title(curr_tf)
    xlabel('WT')
    ylabel('DBDonly')
    set(gca,'FontSize',12)

end

set(gcf,'color','w')


%% Figure 7C - binding on different promoter groups of the DBD
bins_edges = -2:0.5:12;
prom_groups = {'DBD Bound','DBD lost'};
cLimits= {[0,6] , [0, 8]};
numProms = 100;

%choose which group to plot
%cyc8 module
tf_list = {'Ixr1','Phd1','Sok2','Dot6','Skn7','Sko1','Mig1'};
del_order = {'DBDonly','d12','d16'};
%med15 module
% tf_list = {'Msn2','Skn7','Crz1','Gis1','Sko1','Yap1'};
% del_order = {'DBDonly','ddmsn','d11'};

all_combind_zscore = zscore(all_combind_mat);
med_sum = nan(length(tf_list),length(del_order)+1,3);
p = 1;
figure('Position',[2442         583         445         358])

for i = 1:length(tf_list)
    
    %get ordered indices of data
    curr_tf = tf_list{i};
    tf_ind = find(startsWith(all_samples,curr_tf));
    wt_ind = find(strcmp(all_samples,curr_tf));
    X = cellfun(@(c) find(endsWith(all_samples(tf_ind),c)),del_order,'UniformOutput',false);
    Y = vertcat(X{:});
    sorted_ind = tf_ind(Y)';  
    tf_ind = [sorted_ind(1),wt_ind,sorted_ind(2:end)];
    
    %Promoter groups
     WT_data = zscore(chec_struct_med.sum_over_promoter.(curr_tf));
     DBD_data =  zscore(chec_struct_med.sum_over_promoter.([curr_tf,'DBDonly']));     
    [~,WT_sorted]  = sort(WT_data,'descend');
    [~,DBD_sorted]  = sort(DBD_data,'descend');    
    DBD_promoters = setdiff(DBD_sorted(1:numProms),WT_sorted(1:numProms));
    DBD_lost = setdiff(WT_sorted(1:numProms),DBD_sorted(1:numProms));
    promoter_list = {DBD_promoters, DBD_lost};

    for j = 1:length(promoter_list)
        curr_proms = promoter_list{j};
        if isempty(curr_proms)
            p = p+1;
            continue
        end
        mat_subset = all_combind_zscore(curr_proms,tf_ind)';
        med_sum(i,1:length(tf_ind),j) = nanmedian(mat_subset,2)';
        
        med_mat = nanmedian(mat_subset);
        [~,Idx] = sort(med_mat,'descend');  
        mat_subset(:,:) = mat_subset(:,Idx);      

        %plot as histograms
        a(j) = subplot(length(tf_list),length(promoter_list),p);        
        mat_subset(mat_subset>max(bins_edges)) = max(bins_edges);
        mat_subset(mat_subset<min(bins_edges)) = min(bins_edges);        
        data_mat = hist(mat_subset',bins_edges);
        data_mat = smooth_data(data_mat',3);
        joyPlot(data_mat',bins_edges,6,'FaceColor',nanmedian(mat_subset,2),'FaceAlpha',0.95)
        set(gca,'ytick',1:6:length(tf_ind)*6,'YTickLabel',flipud(all_samples(tf_ind)))

        p = p+1;
        caxis(cLimits{j})
        if i == 1
            title(prom_groups{j});
        end
        
        if i == length(tf_list)
            colorbar('location','southoutside')
        end
    end
   colormap(a(1),cm_green)
   colormap(a(2),cm_blues)

end

set(gcf,'color','w')

