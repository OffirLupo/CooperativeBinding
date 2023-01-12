%% Figure 5 - Load data 
% data for this analysis is organized in a matlab struct containing the
% normalized genomic tracks of all samples (out files). this struct was
% created by running the function checStruct_paired on all out files, and
% then running the function mean_of_repeats2 to average repeats

addpath(genpath(cd));
load('DataStructs\chec_struct_med.mat');
load('Datastructs\MNase_seq.mat');
load('GeneralStructs\tf_motifs.mat');
GP = load('GeneralStructs\general_params_130711.mat');

% calculate sum signal on promoter
chec_struct_med = sumOnPro(chec_struct_med,700,GP);
MNase_seq = tss_signal(MNase_seq,700,100,GP);

% load colormap
cm_YlGnBu = cbrewer('seq','YlGnBu',120);
cm_YlGnBu(cm_YlGnBu<0)  = 0;
cm_YlGnBu(cm_YlGnBu>1) = 1;
cm_YlGnBu = cm_YlGnBu(1:100,:);

cm_blues = cbrewer('seq','Blues',120);
cm_blues(cm_blues<0)  = 0;
cm_blues(cm_blues>1) = 1;
cm_blues = cm_blues(1:100,:);

% Create a sumProm matrix of all samples
all_samples = fieldnames(chec_struct_med.sum_over_promoter);
all_combind_mat = zeros(6701,length(all_samples));

for i =1:length(all_samples)   
    all_combind_mat(:,i) = chec_struct_med.sum_over_promoter.(all_samples{i})';
end


[~,Idx] = sort(chec_struct_med.sum_over_promoter.Msn2,'descend');
Msn2_top100 = Idx(1:100);


%% Figure 5A - Plot heatmap top promoters of list of TFs - all together
clear sumProm_list tf_added sorted_proms sumprommat medBindMat

%the order of TFs and and backgrounds
tf_list1 = {'Msn2','Crz1','Skn7','Gis1','Rsf2','Tda9','Yap1'};
tf_list2 = {'Msn2','Crz1_ddmsn','Skn7_ddmsn','Gis1_ddmsn','Rsf2_ddmsn','Tda9_ddmsn','Yap1_ddmsn'};
tf_list3 = {'Msn2_d11','Crz1_d11','Skn7_d11','Gis1_d11','Rsf2_d11','Tda9_d11','Yap1_d11',};

tf_listProms = {'Msn2','Crz1','Skn7','Gis1','Rsf2','Yap1'}; % A list of the TFs to take the top promoters from
tf_list_all  = [{tf_list1}, {tf_list2}, {tf_list3}];

numClusters = 2;
numProms = 100;         
           
% build list based on WTs
for i = 1:length(tf_listProms)
    curr_tf = tf_listProms{i};
    sumProm = chec_struct_med.sum_over_promoter.(curr_tf);
    [~,idx] = sort(sumProm,'descend');
    if i == 1
        sumProm_list = idx(1:numProms);
         tf_added{i} = idx(1:numProms);
         count = numProms;
    else
    sumProm_list = [sumProm_list , idx(1:numProms)];
    sumProm_list = unique(sumProm_list,'stable');
    tf_added{i} =sumProm_list(count+1:length(sumProm_list));
    count = length(sumProm_list);
    end
end
 

%%cluster list based on deletion strain
% build sumProm mat
sumProm_mat_del = zeros(length(tf_list3),6701);    
for i = 1:length(tf_list3)
    sumProm_mat_del(i,:) = zscore(chec_struct_med.sum_over_promoter.(tf_list3{i}));
end
 
%sort promoters
count = 1;
for i  = 1:length(tf_added)
    curr_proms = tf_added{i};
    if length(curr_proms)< 10
        sorted_proms{i} = curr_proms;
        clus_proms{count} = [];
        count = count+1;
    else
    subset = sumProm_mat_del(:,curr_proms);
    Idx = kmeans(subset',numClusters,'Distance','sqeuclidean');
    C = hist(Idx,1:numClusters);
    [~,sortedClusters] = sort(C,'descend');
    X = cellfun(@(c) find(ismember(Idx,c)),num2cell(sortedClusters),'UniformOutput',false);
    for cl = 1:length(X)
        clus_proms{count} = curr_proms(X{cl});
        count = count+1;
    end
    X = vertcat(X{:});
    sorted_proms{i} = curr_proms(X);
    end
end
sorted_proms_all = [sorted_proms{:}];


%plot data on sorted promoters
figure('position',[2408          56         508         914])
for l = 1:length(tf_list_all)
    tf_list = tf_list_all{l};
    
    % build sumProm mat
    sumProm_mat = zeros(length(tf_list),length(sorted_proms_all));            
    for i = 1:length(tf_list)
        curr_data = chec_struct_med.norm.(tf_list{i});
        curr_sumProm = zscore(chec_struct_med.sum_over_promoter.(tf_list{i}));
        sumProm_mat(i,:) = curr_sumProm(sorted_proms_all);
    end
    
    sumprommat{l} = sumProm_mat;
    a = subplot(length(tf_list_all)+2,1,l);
    imagesc(sumProm_mat)
    hold on
    
    % draw lines for clusters
    curr_line = 1;
     for cl = 1:length(tf_added)
         curr_line = curr_line+ length(tf_added{cl});
         plot([curr_line,curr_line],ylim,'-k','LineWidth',1);
     end
     
    cb = colorbar; cb.Label.String = 'Z-score'; cb.Location = 'eastoutside';
    set(gca,'ytick',1:length(tf_list),'yTickLabel',tf_list,'FontSize',10,'TickLabelInterpreter','none');
    caxis([0 5])
    set(gcf,'color','w')

end
colormap(cm_YlGnBu)


%% Figure 5B - Rescue
tf_list = {'Skn7','Yap1','Crz1','Tda9'};
del_order = {'ddmsn','d11','d10'};
numProms = 60;
all_combind_zscore = zscore(all_combind_mat);

figure
for i = 1:length(tf_list)

    curr_tf = tf_list{i};
    tf_1 = 'Msn2';

    tf_ind = find(contains(all_samples,curr_tf));
    X = cellfun(@(c) find(contains(all_samples(tf_ind),c)),del_order,'UniformOutput',false);
    Y = vertcat(X{:});
    tf_ind = [tf_ind(1),tf_ind(Y)'];

    WT_data = zscore(chec_struct_med.sum_over_promoter.(curr_tf));
    dd_data =  zscore(chec_struct_med.sum_over_promoter.([curr_tf,'_ddmsn']));
    [~,Idx] = sort(WT_data - dd_data,'descend');
    proms = intersect(Idx(1:numProms),Msn2_top100);

    mat_subset = all_combind_zscore(proms,tf_ind)';
        
    %sort promoters
    med_mat = nanmedian(mat_subset);
    [~,Idx] = sort(med_mat,'descend');  
    mat_subset(:,:) = mat_subset(:,Idx);      

    subplot(length(tf_list),1,i);        
    imagesc(mat_subset)
    set(gca,'yTick',1:length(tf_ind),'yTickLabel',all_samples(tf_ind),'TickLabelInterpreter','none');
    caxis([0 6])
end
colormap(cm_blues)
set(gcf,'color','w')

%% Figure 5C - Nuclesome occupancy no rescue
data_to_plot = {'WT','ddmsn','Msn2_d11_SM','d12_SM'};
colors_to_plot = {rgb('Navy'),rgb('SkyBlue'),rgb('Crimson'),rgb('Orange'),rgb('Green'),rgb('Gold')};
n_proms = 50; %how many Msn2 targets
window = 400;  %window size around motif

[~,Idx] = sort(chec_struct_med.sum_over_promoter.Msn2,'descend');
proms = Idx(1:n_proms);

%Nucleosomes at promoters
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
title('Mean nucleosome Msn2 targets','FontSize',10); 
set(gcf,'color','w')

subplot('Position',[ 0.1328    0.4946    0.3313    0.0684])
imagesc(mean(MNase_seq.signalTss.Msn2_binding(proms,:)))
colormap(flipud(cm_blues))
set(gca,'xtick',0:100:800,'XTickLabel',[{-1*flip(100:100:700)},{'TSS'},{'+100'}],'XTickLabelRotation',45,'FontSize',12)
set(gca,'ytick',1,'YTickLabel','Msn2')
caxis([0 100])

%Nucleosomes around motifs
motif_table= motif_table_sort('Msn2',{'Msn2_binding'},tf_motifs,MNase_seq,proms);

motif_mat = zeros(size(motif_table,1),window*2,length(data_to_plot));
Msn2_mat = zeros(size(motif_table,1),window*2);

for i = 1:size(motif_table,1)
      curr_chr = motif_table.chr(i);
      curr_loc  = motif_table.loc(i);
      curr_tss  = motif_table.tss_right(i);
      for j = 1:length(data_to_plot)
            curr_data = MNase_seq.norm.(data_to_plot{j}){curr_chr}(curr_loc-window:curr_loc+window-1);
            curr_Msn2 =  MNase_seq.norm.Msn2_binding{curr_chr}(curr_loc-window:curr_loc+window-1);
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
title('Mean nucleosome Msn2 motifs','FontSize',10); 
set(gca,'xtick',[])

subplot('Position',[ 0.5739    0.4875    0.3297    0.0684])
imagesc(mean(Msn2_mat))
cb = colorbar('position',[0.9094    0.4881    0.0211    0.0667]);
colormap(cm_blues)
set(gca,'xtick',0:100:window*2,'XTickLabel',[-1*flip(0:100:window),100:100:window],'XTickLabelRotation',45,'FontSize',12)
set(gca,'ytick',[])



%% Figure 5D - corr to WT in deletion strains
%TFs and deletions
tf_list = {'Skn7','Crz1','Gis1','Rsf2','Yap1','Tda9'};
del_order = {'ddmsn','d4','d7A','d11'};
cm_YlGnSmall = flipud(cbrewer('seq','YlGnBu',length(tf_list)+1));

d_vec = nan(length(del_order),length(tf_list));
figure('Position',[2095         299         332         476])
for i = 1:length(tf_list)
    
    curr_tf = tf_list{i};
    for d = 1:length(del_order)      
    curr_del = [curr_tf,'_',del_order{d}];
        if ismember(curr_del,all_samples)
            d_vec(d,i) = corr(chec_struct_med.sum_over_promoter.(curr_del)',...
                                                chec_struct_med.sum_over_promoter.([curr_tf,'_ddmsn'])');                                           
        else
            continue
        end
    end
    notNanIdx = find(~isnan(d_vec(1:length(del_order),i)));
    plot(notNanIdx+0.5,d_vec(notNanIdx,i)','-o','Color',cm_YlGnSmall(i,:),'LineWidth',0.2,'MarkerFaceColor',cm_YlGnSmall(i,:))
    hold on
end

grid on
ylim([0.4 1])
xlim([0.5 length(del_order)+1.5])
legend(tf_list,'Location','southwest','FontSize',8);
ylabel('Correlation to ddmsn2,4');
set(gca,'XTick',1.5:length(del_order)+0.5,'XTickLabel',del_order,'FontSize',12,'XTickLabelRotation',45,'TickLabelInterpreter','none');
set(gcf,'color','w')

%% Each TF only on its own motif
tf_list = {'Tda9','Yap1','Gis1','Skn7','Rsf2','Crz1'};
data_to_plot = {'_dd','_d11'};

colors_to_plot = {rgb('Black'),rgb('Red')};
cm_motifs = cbrewer('seq','Blues',230);
cm_motifs = cm_motifs(1:200,:);
minLim = 80;
maxLim = 200;

%%%%%%%%
n_prom = 50;
window = 400;
figure('Position',[2284         563         578         199])
[ha, ~] = tight_subplot(1,length(tf_list),[.02 .02],[.1 .1],[.1 .1]);

p =1;
for t = 1:length(tf_list)   
    curr_tf = tf_list{t};

    [~,Idx] = sort(chec_struct_med.sum_over_promoter.([tf_list{t},'_ddmsn']),'descend');
    proms = Idx(1:n_prom);
   
    motif_table = motif_table_sort(curr_tf,{curr_tf},tf_motifs,chec_struct_med,proms);

    wt_ind = find(strcmp(all_samples,curr_tf)); 
    t_ind = find(contains(all_samples,curr_tf)); 
    X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
    Y = vertcat(X{:});
    t_ind =  t_ind(Y)';
    axes(ha(p))
    
    meanVec = zeros(length(t_ind),window*2);
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
    end
    axis square
    title(curr_tf) 
    
    %ste the background color as the maxSignal, limit according the values
    %assigned on top
    yl = ylim;
    plot([window window],yl,'--k','LineWidth',0.1)
    xl = xlim;    
    maxVal = round(yl(2));
    if maxVal > maxLim
       maxVal = maxLim;
    elseif maxVal <= minLim
       maxVal =minLim+1;
    end
    maxVal = maxVal - minLim;
    set(gca,'Color',cm_motifs(maxVal,:))
    
    
    text(xl(2)-200, yl(2)-5,num2str(size(motif_table,1)),'fontsize',8)
    p = p+1;
    set(gca,'xtick',[1,window,window*2-1],'XTickLabel',[-1*(window),0,window],'XTickLabelRotation',45)
    set(gca,'ytick',[])  
       
end
set(gcf,'color','w');

%% plot a matrix of binding around motifs of other TFs 

%colormap
cm_motifs = cbrewer('seq','Blues',230);
cm_motifs = cm_motifs(1:200,:);
minLim = 80;
maxLim = 200;

%%%%%%%%
tf_list = {'Tda9','Yap1','Gis1','Skn7','Rsf2','Crz1'};
motif_tfs ={'Tda9','Yap1','Gis1','Skn7','Rsf2','Crz1','Sko1'};
data_to_plot = {'_ddmsn','_d11'};
colors_to_plot = {rgb('Black'),rgb('Crimson')};
%%%%%%%%
n_prom = 50;
window = 400;
figure('Position',[2287          93         577         447])
[ha, pos] = tight_subplot(length(motif_tfs),length(tf_list),[.005 .005],[.1 .1],[.1 .1]);
meanVec = zeros(1,window*2);

p =1;
for t = 1:length(motif_tfs)   
    
    %find motifs in promoter bound by the current TF
    tf_forMotifs = motif_tfs{t};
   if  strcmp(tf_forMotifs,'Sko1') %dont have data from the ddmsn, so take 
       [~,Idx] = sort(chec_struct_med.sum_over_promoter.(tf_forMotifs),'descend');
        proms = Idx(1:n_prom); 
   else
        [~,Idx] = sort(chec_struct_med.sum_over_promoter.([tf_forMotifs,'_ddmsn']),'descend');
        proms = Idx(1:n_prom);
   end    
   
    motif_table = motif_table_sort(tf_forMotifs,{tf_forMotifs},tf_motifs,chec_struct_med,proms);
    
    %plot the signal of all other TFs on this motif
    for i = 1:length(tf_list)
        curr_tf = tf_list{i};
        t_ind = find(contains(all_samples,curr_tf)); 
        X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
        Y = vertcat(X{:});
        t_ind =  t_ind(Y);
        
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
            set(gca,'ytick',[])
        end
            
        yl = ylim;
        xl = xlim;
        if i == length(tf_list)
             text(xl(2)+100, yl(2)-5,num2str(size(motif_table,1)),'fontsize',8)
        end
        
        %color background
        maxVal = round(yl(2));
        if maxVal > maxLim
           maxVal = maxLim;
        elseif maxVal <= minLim
           maxVal =minLim+1;
        end
        maxVal = maxVal - minLim;
        set(gca,'Color',cm_motifs(maxVal,:)) 
        
        %indicate titles and motif
        p = p+1;
         if i == 1
            ylabel(tf_forMotifs)
         end
         if t == 1
            title(curr_tf)
         end
        %put xticks only at bottom row
         if t <length(motif_tfs)
             set(gca,'xtick',[])
         else
              set(gca,'xtick',[1,window,window*2-1],'XTickLabel',[-1*(window),0,window],'XTickLabelRotation',45)
         end          
    end 
end
set(gcf,'color','w');

%% Figure 5F top Scatter ddmsn vs d11 for Yap1 and Tda9

TF_toPlot = {'Tda9'};

figure
for i = 1:length(TF_toPlot)
    subplot(1,length(TF_toPlot),i)
    xData = chec_struct_med.sum_over_promoter.([TF_toPlot{i},'_ddmsn']);
    yData = chec_struct_med.sum_over_promoter.([TF_toPlot{i},'_d11']);
    cData = chec_struct_med.sum_over_promoter.Crz1;
    scatter(xData,yData,50,cData,'filled')
    axis square
    xlabel([TF_toPlot{i},' ddmsn'])
    ylabel([TF_toPlot{i},' d11'])
    hold on
    plot(xlim,xlim,'--k')
    set(gca,'FontSize',14)
    cb = colorbar;
    cb.Label.String = 'Crz1';
    caxis([0 60000])

end
ylim([0 160000])
xlim([0 160000])

colormap(cm_YlGnBu)
set(gcf,'color','w')

%% Figure 5F bottom

plot_promoterSlim(chec_struct_med,MNase_seq,GP.gene_table.HOR7,{'Crz1','Tda9','Tda9_ddmsn','Tda9_d11'}, ...
    700,5,tf_motifs,GP)



