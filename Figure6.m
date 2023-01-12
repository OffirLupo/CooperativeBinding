%% Figure 6 - Load data 
% data for this analysis is organized in a matlab struct containing the
% normalized genomic tracks of all samples (out files). this struct was
% created by running the function checStruct_paired on all out files, and
% then running the function mean_of_repeats2 to average repeats

addpath(genpath(cd));
load('DataStructs\chec_struct_med.mat');
load('Datastructs\labWTs.mat');
load('GeneralStructs\tf_motifs.mat');
GP = load('GeneralStructs\general_params_130711.mat');

% calculate sum signal on promoter
chec_struct_med = sumOnPro(chec_struct_med,700,GP);
labWTs = sumOnPro(labWTs,700,GP);

% load colormap
cm_YlGnBu = cbrewer('seq','YlGnBu',120);
cm_YlGnBu(cm_YlGnBu<0)  = 0;
cm_YlGnBu(cm_YlGnBu>1) = 1;
cm_YlGnBu = cm_YlGnBu(1:100,:);

cm_YlOrBr = cbrewer('seq','YlOrBr',120);
cm_YlOrBr(cm_YlOrBr<0)  = 0;
cm_YlOrBr(cm_YlOrBr>1) = 1;
cm_YlOrBr = cm_YlOrBr(1:100,:);

% Create a sumProm matrix of all samples
all_samples = fieldnames(chec_struct_med.sum_over_promoter);
all_combind_mat = zeros(6701,length(all_samples));

for i =1:length(all_samples)   
    all_combind_mat(:,i) = chec_struct_med.sum_over_promoter.(all_samples{i})';
end



%% Figure 6A - corr to Cyc8

Cyc8_sumProm = labWTs.sum_over_promoter.Cyc8';
T = readtable('Tables/tableS1.xlsx');
tf_list = T.GeneName(T.ChEC_Profile == 1);
tf_list(strcmp(tf_list,'Cyc8')) = [];
tf_list(strcmpi(tf_list,'CAD1')) = {'YAP2'};
tf_list(strcmpi(tf_list,'CIN5')) = {'YAP4'};

corrToCyc8 = zeros(1,length(tf_list));

for i = 1:length(tf_list)
    curr_tf = lower(tf_list{i});
    curr_tf(1) = upper(curr_tf(1));
    curr_data = labWTs.sum_over_promoter.(curr_tf)';
    corrToCyc8(i) = corr(Cyc8_sumProm,curr_data,'type','pearson');
end

figure
[~,I] = sort(corrToCyc8,'descend');
scatter(length(tf_list)+1 - (1:length(I)),corrToCyc8(I),10,'MarkerFaceColor',rgb('Grey'),'MarkerEdgeColor','none');
ylim([-0.1 1]);
xlabel('TF')
ylabel('Correlation to Cyc8')
set(gcf,'color','w')
axis square
set(gca,'fontsize',12)
gname(tf_list(I))
hold on

tf_subset = {'Yap4','Ixr1','Mig1','Mig2','Mig3','Mot3','Nrg1','Nrg2','Sko1','Sok2','Skn7','Phd1','Dot6','Sut2','Rgt1','Sfl1','Rox1'};
r = find(ismember(tf_list(I), upper(tf_subset)));
b = scatter(length(tf_list)+1 - r,corrToCyc8(I(r)),30,'MarkerFaceColor',rgb('lightGreen'),'MarkerEdgeColor','k');
gname(tf_list(I(r)),b)


%% Figure 6B - corr to WT in deletion strains
tf_list = {'Skn7','Ixr1','Sko1','Dot6','Sok2','Mot3','Cyc8','Mig3','Phd1','Mig2','Nrg1','Nrg2','Mig1','Rox1'};
del_order = {'d6','d12','d16'};

%create a colormap
cm_YlGnSmall = flipud(cbrewer('seq','YlGnBu',length(tf_list)));
cm_YlGnSmall(cm_YlGnSmall<0)  = 0;
cm_YlGnSmall(cm_YlGnSmall>1) = 1;

d_vec = nan(length(del_order),length(tf_list));

figure('Position',[2095         299         332         476])
for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    for d = 1:length(del_order)      
        curr_del = [curr_tf,'_',del_order{d}];
        if ismember(curr_del,all_samples)
            d_vec(d,i) = corr(chec_struct_med.sum_over_promoter.(curr_del)',...
                                                chec_struct_med.sum_over_promoter.(curr_tf)');                                            
        else
            continue
        end
    end
    notNanIdx = find(~isnan(d_vec(:,i)));
    plot(notNanIdx',d_vec(notNanIdx,i)','-o','Color',cm_YlGnSmall(i,:),'LineWidth',0.2,'MarkerFaceColor',cm_YlGnSmall(i,:))
    hold on
    
end
grid on
ylim([0.4 1])
xlim([0.7 3.5])
legend(tf_list,'Location','southwest','FontSize',8);
ylabel('Correlation to WT');
set(gca,'XTick',1:length(del_order),'XTickLabel',del_order,'FontSize',12,'XTickLabelRotation',45,'TickLabelInterpreter','none');
set(gcf,'color','w')


%% Figure 6C - example scatter plots: Sok2 and Nrg1
tf_list = {'Nrg1','Sok2'};
figure
for i = 1:length(tf_list)
    xData = chec_struct_med.sum_over_promoter.(tf_list{i});
    yData = chec_struct_med.sum_over_promoter.([tf_list{i},'_d12']);
    cData = chec_struct_med.sum_over_promoter.Cyc8;
    subplot(1,length(tf_list),i)
    scatter(xData,yData,30,cData,'filled','MarkerEdgeColor','k')
    axis square
    cb = colorbar;
    cb.Label.String = 'Cyc8 binding';
    corrVal = round(corr(xData',yData'),2);
    caxis([0 30000]);
    xlabel(tf_list{i});
    ylabel([tf_list{i},'_d12'])
    set(gca,'FontSize',12)
    hold on
    plot(xlim,xlim,'--k')  
    yl = ylim;
    xl = xlim;
    text(xl(1)+100,yl(2)-5000,['r= ',num2str(corrVal)],'fontsize',12)
end
colormap(cm_YlGnBu);
set(gcf,'color','w')

%% Figure 6D,E Plot heatmap top promoters of list of TFs
clear sumProm_list tf_added sorted_proms sumprommat medBindMat

tf_list1 = {'Cyc8','Nrg1','Nrg2','Rox1','Mig1','Mig2','Mig3','Ixr1','Sok2','Phd1','Sko1','Dot6','Skn7','Mot3',};
tf_list2 = {'Cyc8_d16','Nrg1_d12','Nrg2_d12','Rox1_d16','Mig1_d12','Mig2_d12','Mig3_d12','Ixr1_d16','Sok2_d16','Phd1_d16',...
                'Sko1_d16','Dot6_d16','Skn7_d16','Mot3_d16',};
           
tf_listProms = {'Cyc8','Nrg1','Mig3','Rox1','Mig1','Ixr1','Sok2','Phd1','Sko1','Dot6','Skn7','Mot3'};
tf_list_all  = [{tf_list1}, {tf_list2}];

z_tresh = 3
numClusters = 3;
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
 

%cluster list based on deletion strain
% build sumProm mat
sumProm_mat_del = zeros(length(tf_list2),6701);    
sumProm_mat_wt = zeros(length(tf_list2),6701);            

for i = 1:length(tf_list2)
    sumProm_mat_del(i,:) = zscore(chec_struct_med.sum_over_promoter.(tf_list2{i}));
    sumProm_mat_wt(i,:) = zscore(chec_struct_med.sum_over_promoter.(tf_list1{i}));
end
 
%sort and cluster promoters
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

%plot the binding on the clustered promoter group
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
    
    % draw line for clusters
    curr_line = 1;
     for cl = 1:length(tf_added)
         curr_line = curr_line+ length(tf_added{cl});
         plot([curr_line,curr_line],ylim,'-k','LineWidth',1);
     end
    cb = colorbar; cb.Label.String = 'Z-score'; cb.Location = 'eastoutside';
    set(gca,'ytick',1:length(tf_list),'yTickLabel',tf_list,'FontSize',10,'TickLabelInterpreter','none');
    caxis([1 6])
    set(gcf,'color','w')
    colormap(a,cm_YlGnBu)
    hold on

end

c =subplot(length(tf_list_all)+2,1,l+1);
imagesc(sum(sumprommat{1}>z_tresh));
cb = colorbar; 
cb.Location = 'eastoutside';
set(gca,'ytick',1,'yTickLabel','# bound TFs');
colormap(c,cm_YlOrBr)


%% Figure 6F top - matrix of binding around own motif

tf_list = {'Nrg1','Nrg2','Mig1','Sok2','Skn7','Dot6','Sko1','Mig2','Mig3','Rox1','Phd1','Mot3'};
data_to_plot = {'_d12','_d16'};
%colormap
cm_motifs = cbrewer('seq','Blues',180);
cm_motifs = cm_motifs(1:150,:);
minLim = 30;
maxLim = 150;

%%%%%%%%
n_prom = 50;
window = 400;
figure('Position',[2474          42         363         954])
[ha, ~] = tight_subplot(2,length(tf_list)/2,[.01 .01],[.1 .1],[.1 .1]);

for t = 1:length(tf_list)   
    curr_tf = tf_list{t};

    [~,Idx] = sort(chec_struct_med.sum_over_promoter.(curr_tf),'descend');
    proms = Idx(1:n_prom);
   
    motif_table = motif_table_sort(curr_tf,{curr_tf},tf_motifs,chec_struct_med,proms);

    wt_ind = find(strcmp(all_samples,curr_tf)); 
    t_ind = find(contains(all_samples,curr_tf)); 
    X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
    Y = vertcat(X{:});
    t_ind =  [wt_ind, t_ind(Y)'];
    axes(ha(t))
    
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
            %color of line
            if contains(curr_strain,'d12')
                curr_color = rgb('Salmon');
            elseif contains(curr_strain,'d16')
                curr_color = rgb('red');
            else
                curr_color = 'k';
            end
            plot(meanVec(d,:),'Color',curr_color,'Linewidth',0.25)
            hold on   
            set(gca,'ytick',[])            
    end
    
    axis square
    title(curr_tf)  
    yl = ylim;
    max_signal = round(yl(2));  
    if max_signal < 60
        ylim([0 60])
    end      
    
    if max_signal > maxLim
       max_signal = maxLim;
    elseif max_signal <= minLim
       max_signal =minLim+1;
    end        
    max_signal = max_signal - minLim;
    set(gca,'Color',cm_motifs(max_signal,:))       
    set(gca,'xtick',0:window:window*2,'XTickLabel',[-1*window,0,window],'XTickLabelRotation',45)
end
set(gcf,'color','w');
%% Figure 6F bottom - matrix of binding around motifs of other TFs 

%%%%%%%%
tf_list = {'Nrg1','Nrg2','Mig1','Mig2','Mig3','Rox1','Sok2','Phd1','Skn7','Dot6','Sko1','Mot3'};
motif_tfs = {'Nrg1','Nrg2','Mig1','Mig2','Mig3','Rox1','Sok2','Phd1','Skn7','Dot6','Sko1','Mot3'};
data_to_plot = {'_d12','_d16'};
%%%%%%%%
n_prom = 50;
window = 400;
%%%%%%%%


%colormap
cm_motifs = cbrewer('seq','Blues',180);
cm_motifs = cm_motifs(1:150,:);
minLim = 30;
maxLim = 150;

figure('Position',[ 2373         109        1007         733])
[ha, ~] = tight_subplot(length(motif_tfs),length(tf_list),[.01 .01],[.1 .1],[.1 .1]);

p =1;
for t = 1:length(motif_tfs)   
    tf_forMotifs = motif_tfs{t};
    [~,Idx] = sort(chec_struct_med.sum_over_promoter.(tf_forMotifs),'descend');
    proms = Idx(1:n_prom);
   
    motif_table = motif_table_sort(tf_forMotifs,{tf_forMotifs},tf_motifs,chec_struct_med,proms);
    
    for i = 1:length(tf_list)
        curr_tf = tf_list{i};
        wt_ind = find(strcmp(all_samples,curr_tf)); 
        t_ind = find(contains(all_samples,curr_tf)); 
        X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
        Y = vertcat(X{:});
        t_ind =  [wt_ind, t_ind(Y)'];
        
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
            
            %color of line
            if contains(curr_strain,'d12')
                curr_color = rgb('Salmon');
            elseif contains(curr_strain,'d16')
                curr_color = rgb('red');
            else
                curr_color = 'k';
            end
            plot(meanVec(d,:),'Color',curr_color,'Linewidth',0.25)
            hold on   
            set(gca,'ytick',[])
        end
         yl = ylim;
         max_signal = round(max(max(meanVec,[],2)));
         max_signal = round(yl(2));
        if max_signal < 60
            ylim([0 60])
        end
        
        %set color
        xl = xlim;    
        if max_signal > maxLim
           max_signal = maxLim;
        elseif max_signal <= minLim
           max_signal =minLim+1;
        end
        max_signal = max_signal - minLim;
        set(gca,'Color',cm_motifs(max_signal,:))
        
        %number of sites
        if i == length(tf_list)
             text(xl(2)+100, yl(2)-5,num2str(size(motif_table,1)),'fontsize',8)
        end
        
        %motif TF
         if i == 1
             ylabel(tf_forMotifs)
         end
         
        %binding TF
         if t == 1
            title(curr_tf)
         end
         
         % xticks
         if t <length(motif_tfs)
             set(gca,'xtick',[])
         else
             set(gca,'xtick',0:window:window*2,'XTickLabel',[-1*window,0,window],'XTickLabelRotation',45)
         end
           
           p = p+1;         
    end 
end
set(gcf,'color','w');





