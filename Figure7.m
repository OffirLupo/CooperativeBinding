
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


%% Figure 7B - Example Sko1
xData = chec_struct_med.sum_over_promoter.Sko1';
yData = chec_struct_med.sum_over_promoter.Sko1DBDonly';
cData = chec_struct_med.sum_over_promoter.Sko1_d16;

figure;
scatter(xData,yData,30,cData,'filled','MarkerEdgeColor','k');
axis square; hold on; plot(xlim,xlim,'--k')
set(gcf,'color','w')
xlabel('Sko1');ylabel('Sko1DBDonly'); cb =colorbar; cb.Label.String = 'Sko1 d16';
colormap(cm_YlGnBu)
caxis([ 0 20000])
corrVal = round(corr(xData,yData),2);
yl = ylim;
xl = xlim;
text(100,yl(2)-100,['r= ',num2str(corrVal)],'fontsize',12)

%% Figure 7C top -  heatmap of binding on different promoter groups of the DBD
tf_list = {'Msn2','Skn7','Crz1','Gis1','Sko1','Yap1'};
del_order = {'DBDonly','ddmsn','d11'};
prom_groups = {'DBD Bound','DBD lost'};
cLimits= {[0,6] , [0, 8]};

numProms = 100;
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

        %plot
        a(j) = subplot(length(tf_list),length(promoter_list),p);        
        imagesc(mat_subset)
        set(gca,'yTick',1:length(tf_ind),'yTickLabel',all_samples(tf_ind),'TickLabelInterpreter','none');
        caxis(cLimits{j})
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

%% Figure 7C bottom heatmap of binding on different promoter groups of the DBD
tf_list = {'Ixr1','Phd1','Sok2','Dot6','Skn7','Sko1','Mig1'};
del_order = {'DBDonly','d12','d16'};
prom_groups = {'DBD targets','DBD lost'};
numProms = 100;
all_combind_zscore = zscore(all_combind_mat);
med_sum = nan(length(tf_list),length(del_order)+1,3);
cLimits= {[0,6] , [0, 8]};

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
        
        %plot heatmap
        a(j) = subplot(length(tf_list),length(promoter_list),p);        
        imagesc(mat_subset)
        set(gca,'yTick',1:length(tf_ind),'yTickLabel',all_samples(tf_ind),'TickLabelInterpreter','none');
         
        caxis(cLimits{j})
        p = p+1;
        
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
%% Figure 7D - Correlation to WT vs DBD
figure;
tf_list = {'Msn2','Yap1','Crz1','Skn7','Gis1','Sko1'};
corr_Vecs = nan(length(tf_list),3);

for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    corr_Vecs(i,1) = corr(chec_struct_med.sum_over_promoter.(curr_tf)',chec_struct_med.sum_over_promoter.([curr_tf,'_d11'])');
    corr_Vecs(i,2) = corr(chec_struct_med.sum_over_promoter.([curr_tf,'DBDonly'])',chec_struct_med.sum_over_promoter.([curr_tf,'_d11'])');
    corr_Vecs(i,3) = corr(chec_struct_med.sum_over_promoter.(curr_tf)',chec_struct_med.sum_over_promoter.([curr_tf,'DBDonly'])');  
end

subplot(2,1,1)
scatter(corr_Vecs(:,2), corr_Vecs(:,1),120, corr_Vecs(:,3) - corr_Vecs(:,2),'filled','MarkerEdgeColor','k')

axis square
xlim([0 1]);
ylim([0 1]);
caxis([-0.2 0.2])
hold on
plot(xlim,xlim,'--k')
set(gcf,'color','w')
ylabel('Corr d11 to WT');
xlabel('Corr d11 to DBDonly');
set(gca,'FontSize',12)
cb = colorbar;
cb.Label.String = 'Corr (WT to DBD) - corr(d11 to DBD)';
colormap(cm_rb)
gname(tf_list)

tf_list = {'Sok2','Sko1','Phd1','Skn7','Ixr1','Dot6','Mig1'};
corr_Vecs = nan(length(tf_list),3);
subplot(2,1,2)
for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    if strcmp(curr_tf,'Mig1')
        curr_del = '_d12';
    else
        curr_del = '_d16';
    end
    corr_Vecs(i,1) = corr(chec_struct_med.sum_over_promoter.(curr_tf)',chec_struct_med.sum_over_promoter.([curr_tf,curr_del])');
    corr_Vecs(i,2) = corr(chec_struct_med.sum_over_promoter.([curr_tf,'DBDonly'])',chec_struct_med.sum_over_promoter.([curr_tf,curr_del])');
    corr_Vecs(i,3) = corr(chec_struct_med.sum_over_promoter.(curr_tf)',chec_struct_med.sum_over_promoter.([curr_tf,'DBDonly'])');  
end

scatter(corr_Vecs(:,2), corr_Vecs(:,1),120, corr_Vecs(:,3) - corr_Vecs(:,2),'filled','MarkerEdgeColor','k')
axis square
xlim([0 1]);
ylim([0 1]);
hold on
plot(xlim,xlim,'--k')
set(gcf,'color','w')
xlabel('Corr WT-DBD');
ylabel('Corr d16-DBD');
set(gca,'FontSize',12)
cb = colorbar;
cb.Label.String = 'Corr WT-d16';
colormap(cm_rb)
caxis([-0.2 0.2])
gname(tf_list)


