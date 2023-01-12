%% load truncations data

load('DataStructs\MNase_seq.mat');
load('GeneralStructs\tf_motifs.mat');
GP = load('GeneralStructs\general_params_130711.mat');
load('DataStructs\chec_trun.mat');
chec_trun = sumOnPro(chec_trun,700,GP);
%% load colormaps

cm_ylbu = cbrewer('seq','YlGnBu',120);
cm_ylbu(cm_ylbu<0)  = 0;
cm_ylbu(cm_ylbu>1) = 1;
cm_ylbu = cm_ylbu(1:100,:);


%% Figure 3B,C 200 AA removal - correlation heatmap
all_samples = fieldnames(chec_trun.sum_over_promoter);
TFs = {'Crz1','Hsf1','Skn7'};
dels = {'dMsn4','ddmsn','d200','d150_350','d250_450','d350_550','d450_642'};
n_proms = 40;

[~,I] = sort(chec_trun.sum_over_promoter.Msn2,'descend');
Msn2_top = I(1:n_proms);

figure('Position',[ 2229         440        1505         266]);
for i = 1:length(TFs)
    all_ind = [];
    curr_TF = TFs{i};
    
    %get samples index
    allTFs_ind = find(contains(all_samples,curr_TF));   
    X = cellfun(@(c) find(endsWith(all_samples(allTFs_ind),c)),dels,'UniformOutput',false);
    Y = vertcat(X{:});
    tf_ind = allTFs_ind([1,Y']);
    all_ind = [all_ind,tf_ind'];
    
    all_combined = zeros(6701,length(all_ind));
    for j =1:length(all_ind)
        curr_sample_data = chec_trun.sum_over_promoter.(all_samples{all_ind(j)});
        all_combined(:,j) = curr_sample_data';
    end
    
    [~,effected_in_dd] = sort(chec_trun.sum_over_promoter.(curr_TF) - chec_trun.sum_over_promoter.([curr_TF,'_ddmsn']),'descend');
    effected_in_dd = effected_in_dd(1:n_proms);
    proms = intersect(Msn2_top,effected_in_dd);
  
    corr_mat = corr(all_combined(proms,:),'rows','pairwise');
    
    subplot(1,3,i)
    h = heatmap(round(corr_mat(:,:),2));
    h.FontSize = 10;
    h.YDisplayLabels = strrep(all_samples(all_ind),'_d*','');
    h.NodeChildren(3).YAxis.TickLabelInterpreter = 'none';
    h.Colormap = cm_ylbu;
    set(gcf,'color','w')
end
%% Figure 3C,F - Plot example promoter - UGP1

prom = 'UGP1';
curr_tf = 'Crz1';
dels = {'dMsn4','ddmsn','d200','d150_350','d250_450','d350_550','d450_642'};

allTFs_ind = find(contains(all_samples,curr_tf));   
X = cellfun(@(c) find(endsWith(all_samples(allTFs_ind),c)),dels,'UniformOutput',false);
Y = vertcat(X{:});
tf_ind = allTFs_ind([1,Y']);
plot_promoterSlim(chec_trun,MNase_seq,GP.gene_table.(prom),all_samples(tf_ind),700,tf_motifs,GP)

curr_tf = 'Hsf1';
trun_order = {'dMotA','d135','d200','dTAD','d320','d370','d420','d470','d520','d570','d604'};
allTFs_ind = find(contains(all_samples,curr_tf));   
X = cellfun(@(c) find(endsWith(all_samples(allTFs_ind),c)),trun_order,'UniformOutput',false);
Y = vertcat(X{:});
tf_ind = allTFs_ind([1,Y']);
plot_promoterSlim(chec_trun,MNase_seq,GP.gene_table.(prom),all_samples(tf_ind),700,tf_motifs,GP)
  
%% Figure 3E
numProms = 50;
bins_edges = 8:0.5:16.5;
tf = 'Hsf1';
trun_order = {'dMotA','d135','d200','dTAD','d320','d370','d420','d470','d520','d570','d604'};

all_samples = fieldnames(chec_trun.sum_over_promoter);
all_combined = zeros(6701,length(all_samples));

clear all_combined
for i =1:length(all_samples)
    curr_sample_data = chec_trun.sum_over_promoter.(all_samples{i});
    all_combined(:,i) = curr_sample_data';
end

[~,effected_in_dd] = sort(chec_trun.sum_over_promoter.(tf) - chec_trun.sum_over_promoter.([tf,'_ddmsn']),'descend');
effected_in_dd = effected_in_dd(1:numProms);
[~,Msn2_top] = sort(chec_trun.sum_over_promoter.Msn2, 'descend');
Msn2_top = Msn2_top(1:numProms);
common_promoters = intersect(Msn2_top,effected_in_dd);

tf2_ind = find(contains(all_samples,tf)); 
X = cellfun(@(c) find(contains(all_samples(tf2_ind),c)),trun_order,'UniformOutput',false);
Y = vertcat(X{:});
tf2_ind = [tf2_ind(1),tf2_ind(Y)'];
temp =all_combined(common_promoters,tf2_ind);

figure('Position',[1948         434        1912         399]);
subplot(1,3,1)
imagesc(corr(temp))
set(gca,'xtick',1:length(tf2_ind),'XTickLabel',['WT',trun_order],'XTickLabelRotation',45,'FontSize',12,'TickLabelInterpreter','none')
set(gca,'ytick',1:length(tf2_ind),'yTickLabel',['WT',trun_order])
colorbar

subplot(1,3,2)
temp = log2(temp);
[~,Idx2] = sort(nanmedian(temp,2),'descend');
imagesc(temp(Idx2,:)')
colorbar
caxis([11 16])
set(gca,'ytick',1:length(tf2_ind),'yTickLabel',['WT',trun_order],'FontSize',12,'TickLabelInterpreter','none')
title(tf)
ylabel('promoters')

subplot(1,3,3)
h1 = hist(temp,bins_edges);
h1 = smooth_data(h1,3);
joyPlot(h1,bins_edges,5,'FaceColor',mean(temp,1),'FaceAlpha',0.95)
hold on
set(gca,'ytick',1:5:length(tf2_ind)*5,'yTickLabel',flip(['WT',trun_order]),'FontSize',12,'TickLabelInterpreter','none')
colorbar
colormap(cm_ylbu)
set(gcf,'color','w')

%% Figure 3G - median truncation effect all TFs
tf_1 = 'Msn2';
tf_2_all = {'Hsf1','Crz1','Skn7'};
trun_order = {'dMotA','d135','d200','dTAD','d320','d370','d420','d470','d520','d570','d604'};
OL200_order = {'d150_350','d250_450','d350_550','d450_642'};
numProms = 50;
all_samples = fieldnames(chec_trun.sum_over_promoter);

%colormap for figure
cm_ylGn = cbrewer('seq','YlGnBu',length(trun_order));
cm_set = cm_ylGn([4,6,8,10],:);

%Create a mat with all sumProm data
all_samples = fieldnames(chec_trun.sum_over_promoter);
all_combined_log = zeros(6701,length(all_samples));
for i =1:length(all_samples)
    curr_sample_data = chec_trun.sum_over_promoter.(all_samples{i});
    all_combined_log(:,i) = log2(curr_sample_data');
end

%Msn2 top targets
[~,Msn2_top] = sort(chec_trun.sum_over_promoter.Msn2, 'descend');
Msn2_top = Msn2_top(1:100);
    
figure
for i  = 1:length(tf_2_all)
    
    tf_2 = tf_2_all{i};   
    if strcmp(tf_2,'Skn7')
       trun_order = setdiff(trun_order,'d570','stable'); % noisy sample, remove for proper comparison
    end
    all_order = [trun_order,OL200_order];

    %take Msn2-dependent promoter
    [~,effected_in_dd] = sort(chec_trun.sum_over_promoter.(tf_2) - chec_trun.sum_over_promoter.([tf_2,'_ddmsn']),'descend');
    effected_in_dd = effected_in_dd(1:numProms);
    common_promoters = intersect(Msn2_top,effected_in_dd);
          
    tf1_ind = find(contains(all_samples,tf_1)); % Msn2
    X = cellfun(@(c) find(endsWith(all_samples(tf1_ind),c)),all_order,'UniformOutput',false);
    Y = vertcat(X{:});
    tf1_ind = tf1_ind(Y)';  
    xData = nanmedian(all_combined_log(common_promoters,tf1_ind));

    tf2_ind = find(contains(all_samples,tf_2)); % tf2
    X = cellfun(@(c) find(endsWith(all_samples(tf2_ind),c)),all_order,'UniformOutput',false);
    Y = vertcat(X{:});
    tf2_ind = tf2_ind(Y)';
    yData = nanmedian(all_combined_log(common_promoters,tf2_ind));
    
    subplot(1,3,i)
    scatter(xData(1:length(trun_order)),yData(1:length(trun_order)),60,1:length(trun_order),'filled','MarkerEdgeColor','k');
    hold on
    scatter(xData(length(trun_order)+1:length(all_order)),yData(length(trun_order)+1:length(all_order)),120,(length(trun_order)+1:length(all_order)),'filled','MarkerEdgeColor','k','Marker','s','MarkerFaceAlpha',0.5);

    s = polyfit(xData,yData,1); 
    y = polyval(s,xData); 
    plot(xData,y,'--k')
   
    cb = colorbar;
    cb.Ticks = 1:length(all_order);
    cb.TickLabels = all_order;
    cb.TickLabelInterpreter = 'none';
       
    axis square
    xlabel(tf_1)
    ylabel(tf_2)
    set(gca,'fontsize',12)
    hold on
    xl = xlim;
    yl = ylim;
    text(xl(1)+0.1,yl(2)-0.1,['R = ' ,num2str(round(corr(xData',yData'),2))],'FontSize',12)
    legend({'Truncations','200AA'},'Location','southeast')
    gname(all_order)
    colormap(vertcat(cm_ylGn,cm_set))

end
set(gcf,'color','w')
colormap(vertcat(cm_ylGn,cm_set))

