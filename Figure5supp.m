%% Figure S5 - Load data 
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


%% Effect of VMA5 promoter
tf_list = {'Tda9','Asg1','Gis1','Skn7'};
%tf_list = all_samples(find(contains(all_samples,'VMA5')));
figure;

for i = 1:length(tf_list)

    subplot(2,length(tf_list)/2,i)
   % curr_tfWT = strsplit(tf_list{i},'_VMA5');
    xData = chec_struct_med.sum_over_promoter.(tf_list{i});
    yData = chec_struct_med.sum_over_promoter.([tf_list{i},'_VMA5']) ;
    scatter(xData,yData,20,'filled','MarkerFaceColor',rgb('lightBlue'))
    title(tf_list{i});
    axis square
    xlabel('Native promoter')
    ylabel('VMA5 promoter','FontSize',12)
    hold on
    yl = ylim;
    plot(xlim,xlim,'--k','LineWidth',1)
    set(gca,'fontsize',12);
    set(gcf,'color','w')
    text(100, yl(2),['R = ',num2str(round(corr(xData',yData'),2))]);

% hold on
% scatter(xData(GP.gene_table.UGP1),yData(GP.gene_table.UGP1),20,'filled','MarkerFaceColor',rgb('red'))
% scatter(xData(GP.gene_table.TPS2),yData(GP.gene_table.TPS2),20,'filled','MarkerFaceColor',rgb('red'))
end

%% Figure S5B - corr to WT in deletion strains
%TFs and deletions

tf_list = {'Msn2','Skn7','Msn4','Sko1','Crz1','Gis1','Rsf2','Tda9','Yap1'};
del_order = {'dmsn4','dmsn2','ddmsn','d4','d7A','d11','d10_wMsn2'};
cm_YlGnSmall = flipud(cbrewer('seq','YlGnBu',length(tf_list)+1));

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
    notNanIdx = find(~isnan(d_vec(1:length(del_order),i)));
    plot(notNanIdx,d_vec(notNanIdx,i)','-o','Color',cm_YlGnSmall(i,:),'LineWidth',0.2,'MarkerFaceColor',cm_YlGnSmall(i,:))
    hold on
end
axis tight
grid on
ylim([0.4 1])
xlim([0.5 length(del_order)+1.5])
legend(tf_list,'Location','southwest','FontSize',8);
ylabel('Correlation to ddmsn2,4');
set(gca,'XTick',1:length(del_order),'XTickLabel',del_order,'FontSize',12,'XTickLabelRotation',45,'TickLabelInterpreter','none');
set(gcf,'color','w')

