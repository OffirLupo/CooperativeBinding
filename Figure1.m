%% Figure 1 - Load data 
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
% create colormaps
cm_green = cbrewer('seq','Greens',120);
cm_green(cm_green<0)  = 0;
cm_green(cm_green>1) = 1;
cm_green = cm_green(1:100,:);

cm_YlGn = cbrewer('seq','YlGn',120);
cm_YlGn(cm_YlGn<0)  = 0;
cm_YlGn(cm_YlGn>1) = 1;
cm_YlGn = cm_YlGn(1:100,:);

cm_YlOrBr = cbrewer('seq','YlOrBr',120);
cm_YlOrBr(cm_YlOrBr<0)  = 0;
cm_YlOrBr(cm_YlOrBr>1) = 1;
cm_YlOrBr = cm_YlOrBr(1:100,:);


%% Figure 1B - TF -selection  
T = readtable('Tables/tableS1.xlsx');
all_TFs = T.GeneName(T.ChEC_Profile == 1);
all_TFs(strcmpi(all_TFs,'CAD1')) = {'YAP2'}; %rename the TF
all_TFs(strcmpi(all_TFs,'CIN5')) = {'YAP4'};
all_TFs(strcmpi(all_TFs,'MSN2')) = []; 

Msn2_sumProm = labWTs.sum_over_promoter.Msn2';
Gal11_sumProm = labWTs.sum_over_promoter.Gal11';

corrToMsn2 = zeros(1,length(all_TFs));
corrToGal11 = zeros(1,length(all_TFs));


for i = 1:length(all_TFs)
    curr_tf = lower(all_TFs{i});
    curr_tf(1) = upper(curr_tf(1));
    curr_data = labWTs.sum_over_promoter.(curr_tf)';
    corrToMsn2(i) = corr(Msn2_sumProm,curr_data);
    corrToGal11(i) = corr(Gal11_sumProm,curr_data);
end

figure;
scatter(corrToMsn2,corrToGal11,80,'MarkerFaceColor',rgb('Grey'),'MarkerEdgeColor','none');
axis square
xlim([-0.1 1]); ylim([-0.1 1]);
set(gcf,'color','w')
xlabel('Correlation to Msn2')
ylabel('Correlation to Gal11');
set(gca,'fontsize',12)
hold on
gname(all_TFs)

tf_subset = {'Yap1', 'Skn7', 'Crz1','Msn4','Mss11','Sko1','Hot1','Rsf2','Tda9','Asg1','Gis1','Adr1','Ixr1','Nfi1','Hsf1'};
r = find(ismember(all_TFs, upper(tf_subset)));
b = scatter(corrToMsn2(r),corrToGal11(r),80,'MarkerFaceColor',rgb('lightGreen'),'MarkerEdgeColor','k');
gname(all_TFs(r),b)


%% Figure 1C - plot example scatter plots  
TF_toPlot = {'Tda9','Crz1','Hsf1'};

figure
for i = 1:length(TF_toPlot)
    
    subplot(1,length(TF_toPlot),i)
    xData = chec_struct_med.sum_over_promoter.Msn2;
    yData = chec_struct_med.sum_over_promoter.(TF_toPlot{i});
    cData = chec_struct_med.sum_over_promoter.Gal11;  
    scatter(xData,yData,20,cData,'filled','MarkerEdgeColor','k','LineWidth',0.1)
    corrVal = round(corr(xData',yData','type','pearson'),2);
    yl = ylim;
    xl = xlim;
    text(xl(2)-10000,yl(2)-100,['R= ',num2str(corrVal)],'fontsize',12)
    axis square
    xlabel('Msn2')
    ylabel(TF_toPlot{i})
    set(gca,'FontSize',12)
   caxis([0 40000])
   set(gca,'fontsize',12)

end

cb = colorbar;
cb.Label.String = 'Med15';
colormap(cm_YlGn)
set(gcf,'color','w')


%% Figure 1D - plot the length of the IDR for all TFs
T = readtable('Tables/tableS1.xlsx');
tf_list = {'Msn2','Msn4','Nfi1','Adr1','Asg1','Mss11','Rsf2','Tda9','Ixr1','Yap1','Crz1','Skn7','Gis1','Hsf1','Sko1','Hot1'};

all_families = unique(T.DBD_YETFASCO);
X = cellfun(@(c) find(contains(T.DBD_YETFASCO,c)),all_families,'UniformOutput',false);

%sort by num of TF in group
[~,Idx] = sort(cellfun(@length, X),'descend');
%move Zinc cluster (largest group) to the end
Idx(length(Idx)+1) = Idx(1);
Idx = Idx(2:end);
X = X(Idx);
groups_length = cellfun(@length, X);
all_families = all_families(Idx);

%sort family by length
for i  = 1:length(X)
    curr_ind = X{i};
    [~,a] = sort(T.IDR_length(curr_ind),'descend');
    X{i} = curr_ind(a);
end

all_idx_sorted = vertcat(X{:});
T = T(all_idx_sorted,:);

figure('position',[ 2201         452         818         205])
%plot all TFs
bar(T.IDR_length,'FaceColor',rgb('lightGrey'),'EdgeColor',rgb('lightGrey'),'FaceAlpha',0.9)
axis tight
hold on
bar(find(T.ChEC_Profile),T.IDR_length(find(T.ChEC_Profile)),'FaceColor',rgb('Grey'),'EdgeColor',rgb('Grey'))
yl = ylim;

%plot selected TFs
X = cellfun(@(c) find(contains(T.GeneName,c)),upper(tf_list),'UniformOutput',false);
Y = vertcat(X{:});
bar(Y,T.IDR_length(Y),'FaceColor',rgb('lightGreen'),'EdgeColor','k')
set(gcf,'color','w')
text(Y,T.IDR_length(Y)+0.05,tf_list,'Color','k')
ylabel('IDR length (AA)')
xlabel('TFs')
set(gca,'fontsize',12)

c = 0;
for i = 1:length(groups_length)
    c = c + groups_length(i);
    b= plot([c,c],yl,'--k')   ; 
    text(c,yl(2),all_families(i),'Rotation',90)
end 
plot(xlim,[median(T.IDR_length),median(T.IDR_length)],'--k')

legend({'All','profiled','Selected'},'location','northeast','FontSize',8)

%% Figure 1E - Plot heatmap top promoters of list of TFs
tf_list = {'Msn2','Msn4','Gal11','Nfi1','Adr1','Asg1','Mss11','Rsf2','Tda9','Ixr1','Yap1','Crz1','Skn7','Gis1','Hsf1','Sko1','Hot1'};
z_tresh = 3;
clear sumProm_list

% take the top promoters of each selected TF
sumProm_all = zeros(6701,length(tf_list));
numProms = 100;

for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    sumProm = chec_struct_med.sum_over_promoter.(curr_tf);
    sumProm_all(:,i) = sumProm';
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

%build a sumProm mat on this list
sumProm_mat = zeros(length(sumProm_list),length(tf_list));            
for i = 1:length(tf_list)
    curr_data = chec_struct_med.norm.(tf_list{i});
    curr_sumProm = zscore(chec_struct_med.sum_over_promoter.(tf_list{i}));
    sumProm_mat(:,i) = curr_sumProm(sumProm_list);
end
sumProm_mat = sumProm_mat';

[~,Idx] = sort(nanmedian(sumProm_mat(:,1:numProms)),'descend');
sumProm_mat(:,1:numProms) = sumProm_mat(:,Idx);


figure
a = subplot(13,1,1:6);
imagesc(sumProm_mat)
cb = colorbar; cb.Label.String = 'Z-score'; cb.Location = 'eastoutside';
set(gca,'ytick',1:length(tf_list),'yTickLabel',tf_list,'FontSize',10,'TickLabelInterpreter','none');
caxis([1 6])
set(gcf,'color','w')
colormap(a,cm_green)
hold on
set(gca,'xtick',[]);

c =subplot(13,1,7);
imagesc(sum(sumProm_mat>z_tresh));
cb = colorbar; 
cb.Location = 'eastoutside';
set(gca,'ytick',1,'yTickLabel','# bound TFs');
colormap(c,cm_YlOrBr)

%% Figure 1F - Plot example promoters
tf_list = {'Msn2','Msn4','Gal11','Nfi1','Adr1','Asg1','Mss11','Rsf2','Tda9','Ixr1','Yap1','Crz1','Skn7','Gis1','Hsf1','Sko1','Hot1'};

plot_promoterSlim(chec_struct_med,MNase_seq,GP.gene_table.YPL230W,tf_list,700,tf_motifs,GP)
plot_promoterSlim(chec_struct_med,MNase_seq,GP.gene_table.SED1,tf_list,700,tf_motifs,GP)


