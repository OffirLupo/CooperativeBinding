%% Figure S2 supplementary
% data for this analysis is organized in a matlab struct containing the
% normalized genomic tracks of all samples (out files). this struct was
% created by running the function checStruct_paired on all out files, and
% then running the function mean_of_repeats2 to average repeats

addpath(genpath(cd));
load('DataStructs\chec_struct_med.mat');
load('DataStructs\labWTs.mat');
load('GeneralStructs\tf_motifs.mat');
GP = load('GeneralStructs\general_params_130711.mat');
load('GeneralStructs\opnScore.mat');
gene_names = GP.gene_infoR64.name;
for i = 1:6701
    if isempty(gene_names{i})
       gene_names{i} = 'nan';
    end
end


%calculate sum signal on promoter
chec_struct_med = sumOnPro(chec_struct_med,700,GP);
labWTs = sumOnPro(labWTs,700,GP);

%colormaps for figure
cm_blues = cbrewer('seq','Blues',100);
cm_blues(cm_blues<0)  = 0;
cm_blues(cm_blues>1) = 1;
cm_blues = cm_blues(1:90,:);

cm_YlGnBu = cbrewer('seq','YlGnBu',120);
cm_YlGnBu(cm_YlGnBu<0)  = 0;
cm_YlGnBu(cm_YlGnBu>1) = 1;
cm_YlGnBu = cm_YlGnBu(1:100,:);

%% S2A - Effect of single deletions Msn2/4 in common and unique promoters

tf_list = {'Yap1','Hsf1','Crz1','Rsf2','Tda9','Skn7','Gis1','Asg1'};
del_list = {'','','_dmsn4','_dmsn2','_ddmsn'};

bins_edges = 9:0.5:18;
n_prom = 50;

p = 1;
figure('Position',[2051         473        1541         298])
for i = 1:length(tf_list)
    
    curr_tf =tf_list{i}; 
    [common_promoters, unique_promoters] = diff_proms(chec_struct_med,'Msn2',curr_tf,n_prom,100);
    promoter_list = [common_promoters ,unique_promoters];

    del_dataMat = zeros(length(del_list),length(promoter_list));
    
    del_dataMat(1,:) = log2(chec_struct_med.sum_over_promoter.Msn2(promoter_list));
    
    for d = 2:length(del_list)
        curr_del = curr_tf;
        if d == 2
        curr_data = log2(chec_struct_med.sum_over_promoter.(curr_del));
        del_dataMat(d,:) = curr_data(promoter_list)';   
        end
        if d >2
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
    
    %open new figure after 4 TFs
    if i == 5
        figure('Position',[2052         141        1541         298])
        p = 1;
    end
    a = subplot(2,length(tf_list),[p:p+1]);
    imagesc(del_dataMat)
    set(gca,'ytick',1:5,'yTickLabel',{'Msn2','WT','dmsn4','dmsn2','ddmsn'},'FontSize',8,'TickLabelInterpreter','none');
    xlabel('Promoters')
    title(curr_tf)
    set(gcf,'color','w')
    colormap(a,cm_blues)
    caxis([12 16])
    
    %Common 
    b = subplot(2,length(tf_list),p+length(tf_list));
    proms_idx = 1:n_prom;
    data_mat = hist(del_dataMat(2:length(del_list),proms_idx)',bins_edges);
    data_mat = smooth_data(data_mat',3);
    joyPlot(data_mat',bins_edges,3,'FaceColor',nanmedian(del_dataMat(2:length(del_list),proms_idx),2),'FaceAlpha',0.95)
    set(gca,'ytick',1:3:12,'YTickLabel',{'ddmsn','dmsn4','dmsn2','WT'})
    caxis([13 14.5]) 
    colormap(b,cm_blues)
    ylabel(curr_tf)
    
    %Unique
    c = subplot(2,length(tf_list),p+length(tf_list)+1);
    proms_idx = n_prom+1:n_prom*2;
    data_mat = hist(del_dataMat(2:length(del_list),proms_idx)',bins_edges);
    data_mat = smooth_data(data_mat',3); 
    data_mat = smooth_data(data_mat,3);
    joyPlot(data_mat',bins_edges,3,'FaceColor',nanmedian(del_dataMat(2:length(del_list),proms_idx),2),'FaceAlpha',0.95)
    set(gca,'ytick',1:3:12,'YTickLabel',{'ddmsn','dmsn4','dmsn2','WT'})
    colormap(c,cm_blues)

    caxis([13 15]) 
    p = p+2;
end
set(gcf,'color','w')


%% Figure S2B,C - Msn2 with Rpn4 DBD
tf_list = {'Msn2','Crz1','Hsf1','Rsf2','Tda9','Skn7','Yap1'};
figure('position',[1918         430        1470         439]);

axis_ind = [{[1,2,6,7]},num2cell(3:5),num2cell(8:10)];
for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    if i == 1
        xData = chec_struct_med.sum_over_promoter.(curr_tf);
        yData = chec_struct_med.sum_over_promoter.([curr_tf,'_Rpn4ZF']);
    else
        xData = chec_struct_med.sum_over_promoter.([curr_tf,'_dmsn2']);
        yData = chec_struct_med.sum_over_promoter.([curr_tf,'_Msn2zfRpn4']);
    end
     cData = labWTs.sum_over_promoter.Rpn4;
 
    subplot(2,5,axis_ind{i})
    scatter(log2(xData),log2(yData),10,log2(cData),'filled')
    yl = ylim;
    text(yl(1), yl(2) ,['R =  ' , num2str(round(corr(xData',yData'),2))])
    cb = colorbar;
    cb.Label.String = 'Rpn4';
    set(gcf,'color','w')
    set(gca,'fontsize',12)
    caxis([8 15])
    xlabel('WT');
    ylabel('Msn2 wRpn4ZF');
    title(curr_tf)  
end

colormap(cm_YlGnBu)


%% Figure 2D -  Msn2 fails ro recruit itself to Rpn4 targets

figure
xData = chec_struct_med.sum_over_promoter.Msn2;
yData = chec_struct_med.sum_over_promoter.Msn2_Msn2wRpn4DBD_A;
cData = labWTs.sum_over_promoter.Rpn4;
scatter(log2(xData),log2(yData),10,log2(cData),'filled')
axis square
yl = ylim;
text(yl(1), yl(2) ,['R =  ' , num2str(round(corr(xData',yData'),2))])
cb = colorbar;
cb.Label.String = 'Rpn4';
set(gcf,'color','w')
set(gca,'fontsize',12)
caxis([8 15])
xlabel('WT');
ylabel('Msn2 wRpn4ZF');
title('Msn2')  

colormap(cm_YlGnBu)

%% Figure 2E - Nuclesome occupancy at full promoter
data_to_plot = {'WT','dmsn4','ddmsn'};
tf = 'Crz1';
colors_to_plot = {rgb('Navy'),rgb('SkyBlue'),rgb('Crimson'),rgb('Orange'),rgb('Green'),rgb('Gold')};
n_proms = 100;
window = 400;

[~, proms] = diff_proms(chec_struct_med,'Msn2',tf,n_proms,100);
% [~,Idx] = sort(chec_struct_med.sum_over_promoter.(tf),'descend');
% proms = Idx(1:n_proms);

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





%% Figure 2F - correlation matrix
tf_list1 = {'Gal11','Msn2','Gis1','Skn7','Rsf2','Tda9','Crz1','Yap1','Hsf1'};
tf_list2 = {'Gal11','Msn2_dgal11','Gis1','Skn7_dGal11','Rsf2_dGal11','Tda9_dGal11','Crz1_dGal11','Yap1_dGal11','Hsf1_dGal11'};
tf_list3 = {'Gal11_ddmsn','Msn2','Gis1_ddmsn','Skn7_ddmsn','Rsf2_ddmsn','Tda9_ddmsn','Crz1_ddmsn','Yap1_ddmsn','Hsf1_ddmsn'};
tf_list4 = {'Gal11_ddmsn','Msn2_d11','Gis1_d11','Skn7_d11','Rsf2_d11','Tda9_d11','Crz1_d11','Yap1_d11','Hsf1_ddmsn'};

titles = {'WT','dmed15','ddmsn2,4','d11'};
tf_list_all = {tf_list1, tf_list2 ,tf_list3,tf_list4};

all_samples = fieldnames(chec_struct_med.sum_over_promoter);
all_combind_mat = zeros(6701,length(all_samples));

all_samples(strcmp(all_samples,'Msn2nonDBD_d11')) = [];
all_samples(strcmp(all_samples,'Msn2DBDonly_d11')) = [];

for i =1:length(all_samples)   
    all_combind_mat(:,i) = chec_struct_med.sum_over_promoter.(all_samples{i})';
end

figure('position',[ 2034         382        1460         420])
for i  = 1:length(tf_list_all)
    curr_list = tf_list_all{i};
    X = cellfun(@(c) find(strcmp(all_samples,c)),curr_list,'UniformOutput',false);
    X = vertcat(X{:});   
    
    subplot(1,4,i)
    imagesc(corr(all_combind_mat(:,X),'rows','pairwise'));
    plotgrid(corr(all_combind_mat(:,X)),1,1)
    set(gca,'xtick',1:length(X),'xticklabel',all_samples(X),'xticklabelrotation',45,'TickLabelInterpreter','none')
    set(gca,'ytick',1:length(X),'yticklabel',all_samples(X),'TickLabelInterpreter','none')
    title(titles{i})
    caxis([ 0.2 1])
    axis square
    
end
cb = colorbar('position',[    0.9215    0.3476    0.0097    0.3381]);
cb.Label.String = 'Corr. promoter preference';
set(gcf,'color','w')
colormap(cm_YlGnBu)


%% plot a matrix of binding around motifs of Msn2, TF in common and TF in unique 
% 
% % general parameters%
% tf_list = {'Skn7','Crz1','Tda9','Rsf2','Yap1','Hsf1'};
% data_to_plot = {'_ddmsn','_dGal11'}; %which background, can add also d7,d11...
% motif_labels = {'Msn2'};
% colors_to_plot = {rgb('Black'),rgb('Crimson'),rgb('Orange')};
% n_prom = 100;
% window = 400;
% %%%%%%%%%%%
% 
% all_tables = [];
% meanVec = [];
% all_samples = fieldnames(chec_struct_med.sum_over_promoter);
% 
% c = 1;
% figure('Position',[2028         569        1263         153])
% [ha, ~] = tight_subplot(1,length(tf_list),[.01 .03],[.1 .1],[.1 .1]);
% 
%  
% for t = 1:length(tf_list)   
%     curr_tf = tf_list{t};
%     t_ind = find(contains(all_samples,curr_tf)); % tf2
%     X = cellfun(@(c) find(contains(all_samples(t_ind),c)),data_to_plot,'UniformOutput',false);
%     Y = vertcat(X{:});
%     t_ind = [t_ind(1), t_ind(Y)'];   
%     
%     %Msn2 motifs
%     [~,Idx] = sort(chec_struct_med.sum_over_promoter.Msn2,'descend');
%     proms = Idx(1:n_prom); 
%     motif_table = motif_table_sort('Msn2',{'Msn2'},tf_motifs,chec_struct_med,proms);
%     all_tables.M1 = motif_table;
%       
%     all_table_names = fieldnames(all_tables);
% 
%     for i = 1:length(all_table_names)
%         axes(ha(c))
% 
%         curr_table = all_tables.(all_table_names{i});
%         sumMat = zeros(length(t_ind),1);
%         
%         for d = 1:length(t_ind)
%             motif_mat = zeros(size(curr_table,1),window*2);
%             curr_strain = all_samples{t_ind(d)};
%             for m = 1:size(curr_table,1)
%                   curr_chr = curr_table.chr(m);
%                   curr_loc  = curr_table.loc(m);
%                   curr_tss  = curr_table.tss_right(m);
% 
%                   curr_data = zscore(chec_struct_med.norm.(curr_strain){curr_chr}(curr_loc-window:curr_loc+window-1));
%                    if curr_tss == 0
%                      curr_data = flip(curr_data);
%                     end
%                   motif_mat(m,:) = curr_data ;                            
%             end
%             
%             motif_mat = smooth_data(motif_mat,3);            
%             meanVec(d,:) = smooth_offir(nanmean(motif_mat(:,:)),10);
%             plot(meanVec(d,350:450),'Color',colors_to_plot{d},'Linewidth',0.2)
%             hold on        
%             axis tight          
%         end
%               
%          title(curr_tf)
%         
%         set(gca,'xtick',1:50:101,'XTickLabel',{'-50','0','50'},'XTickLabelRotation',45)
% 
%         c = c+1;
%         
%     end 
%     
% end
% 
% set(gcf,'color','w');