%% Figure S1 - Load data 
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
% create colormaps
cm_green = cbrewer('seq','Greens',120);
cm_green(cm_green<0)  = 0;
cm_green(cm_green>1) = 1;
cm_green = cm_green(1:100,:);

cm_RdGy = cbrewer('div','RdGy',120);
cm_RdGy(cm_RdGy<0)  = 0;
cm_RdGy(cm_RdGy>1) = 1;
cm_RdGy = cm_RdGy(1:100,:);

cm_YlGnBu = cbrewer('seq','YlGnBu',120);
cm_YlGnBu(cm_YlGnBu<0)  = 0;
cm_YlGnBu(cm_YlGnBu>1) = 1;
cm_YlGnBu = cm_YlGnBu(1:100,:);

%% Figure S1A
T = readtable('Tables/tableS1.xlsx');
all_TFs = T.GeneName(T.ChEC_Profile == 1);
all_TFs(strcmpi(all_TFs,'CAD1')) = {'YAP2'};
all_TFs(strcmpi(all_TFs,'CIN5')) = {'YAP4'};
all_TFs(strcmpi(all_TFs,'NFI1')) = [];

Msn2_sumProm = labWTs.sum_over_promoter.Msn2';
Gal11_sumProm = labWTs.sum_over_promoter.Gal11';

corrToMsn2 = zeros(1,length(all_TFs));
corrToGal11 = zeros(1,length(all_TFs));
sumProm_mat = zeros(6701,length(all_TFs)+2);
sumProm_mat(:,1) = Msn2_sumProm;
sumProm_mat(:,2) = Gal11_sumProm;

for i = 1:length(all_TFs)
    curr_tf = lower(all_TFs{i});
    curr_tf(1) = upper(curr_tf(1));
    curr_data = labWTs.sum_over_promoter.(curr_tf)';
    sumProm_mat(:,i+2) = curr_data;
    corrToMsn2(i) = corr(Msn2_sumProm,curr_data);
    corrToGal11(i) = corr(Gal11_sumProm,curr_data);
end

figure
[~,I] = sort(corrToGal11,'descend');
scatter(length(I)+1 - (1:length(I)),corrToGal11(I),10,'MarkerFaceColor',rgb('Grey'),'MarkerEdgeColor','K');
ylim([-0.1 1]);
xlabel('TF rank')
ylabel('Correlation to Med15')
set(gcf,'color','w')
axis square
set(gca,'fontsize',12)
gname(all_TFs(I))

hold on
tf_subset = {'Msn2','Msn4'};
r = find(ismember(all_TFs(I), upper(tf_subset)));
b = scatter(length(I)+1 - r,corrToGal11(I(r)),30,'MarkerFaceColor',rgb('lightGreen'),'MarkerEdgeColor','k');
gname(all_TFs(I(r)),b)


%% Figure S1B compare motifs
%need to load all cisBP motifs for comparisson
%PART I: build a temporal struct of specific TFs and do motif-enrichment analysis (7mer)
clear tempstruct 
tf_list = {'Msn2','Msn4','Adr1','Asg1','Rsf2','Tda9','Crz1','Skn7','Gis1','Sko1','Ixr1','Yap1'};

for i =1:length(tf_list)
    tempstruct.norm.(tf_list{i}) = chec_struct_med.norm.(tf_list{i});
end

tempstruct = motifFind(tempstruct,7,GP,1:6701);

%PART II: build PWMs and plot 7Mer vs. cisBP seqlogos

figure('Position',[2273         204         758         549])
[ha, pos] = tight_subplot(length(tf_list),2,[.01 .08],[.1 .1],[.1 .1]);
T = importdata('GeneralStructs\cisBP\TF_Information.txt');
T = cell2table(T.textdata(2:end,1:27),'VariableNames',T.textdata(1,1:27));
p = 1;

for i  = 1:length(tf_list)*2
    %plot 7mer seqlogo of the top 10 motifs for each TF
    axes(ha(p)); 
    seqlogoPlot(tempstruct,tf_list{i},'mer7',10);
    if p == 1 
        title('Data - 7mer','FontSize',6)
    end
    
    if i < length(tf_list)
     set(gca,'xtick',[])
    end
    p = p+1;
 
    %plot the seqlogo of the PWM from cisBP
    axes(ha(p));   
    curr_tf_ind = find(strcmpi(T.TF_Name,tf_list{i}));
    if isempty(curr_tf_ind)
        p = p+1;
        continue
    end
    
    JAS_ind = strcmp(T.MSource_Identifier(curr_tf_ind),'JASPAR');
    PBM_ID = strcmp(T.Motif_Type(curr_tf_ind),'PBM');
    DeBoer11_ID = strcmp(T.MSource_Identifier(curr_tf_ind),'DeBoer11');
    
    if strcmp(tf_list{i},'Crz1') 
        PBM_ID = 0;
    elseif strcmp(tf_list{i},'Yap1')
        PBM_ID = 0;
        JAS_ind = 0;
    end

    if  sum(PBM_ID)>0
        curr_tf_ind = curr_tf_ind(PBM_ID);
         curr_ID =  T.Motif_ID(curr_tf_ind);
         pwm = importdata(['GeneralStructs\cisBP\pwms_all_motifs\',curr_ID{1},'.txt']);
         pwm = pwm.data(:,2:end)';
         pwm = pwm(:,median(pwm)<0.2);
         SeqLogoFig(pwm,'CUTOFF',0)
         ylim([0 2.5])
         text(1,2.2,'PBM')
         
    elseif sum(JAS_ind)>0
         curr_tf_ind = curr_tf_ind(JAS_ind);
         curr_ID =  T.Motif_ID(curr_tf_ind);
         pwm = importdata(['GeneralStructs\cisBP\pwms_all_motifs\',curr_ID{1},'.txt']);
         SeqLogoFig(pwm.data(:,2:end)','CUTOFF',0)
         ylim([0 2.5])
         text(1,2.2,'JASPAR')
         
    elseif sum(DeBoer11_ID)>0
         curr_tf_ind = curr_tf_ind(DeBoer11_ID);
         curr_ID =  T.Motif_ID(curr_tf_ind(2));
         pwm = importdata(['GeneralStructs\cisBP\pwms_all_motifs\',curr_ID{1},'.txt']);
         pwm =pwm.data(:,2:end)';
         % do reverse complement
         pwm = rot90(pwm,2) ;
         SeqLogoFig(pwm,'CUTOFF',0)
         ylim([0 2.5])
         text(1,2.2,'YETFASCO')
    end

    ylabel('')
    if i < length(tf_list)-1
     set(gca,'xtick',[])
    end

    p = p+1;
    
end
%% Figure S1C find # bound motif each tf

window = 100;
%z_tres = 1;
tf_list = {'Msn2','Msn4','Adr1','Rsf2','Tda9','Crz1','Skn7','Gis1','Sko1','Yap1','Hsf1'};
figure('position',[ 1966         357        1871         275])

for i = 1:length(tf_list)
    
    curr_tf = tf_list{i};
    motif_table = motif_table_sort(curr_tf,{curr_tf},tf_motifs,labWTs,1:6701);
    motif_mat = nan(size(motif_table,1),window*2);
   
    for m = 1:size(motif_table,1)
          curr_chr = motif_table.chr(m);
          curr_loc  = motif_table.loc(m);
          curr_strand  = motif_table.motif_plusStrand(m);
          curr_data = labWTs.norm.(curr_tf){curr_chr}(curr_loc-window:curr_loc+window-1);
        if curr_strand == 0
            curr_data = flip(curr_data);
        end
            motif_mat(m,:) = curr_data ;                            
    end

    [~,Idx] = sort(nansum(motif_mat,2));    
    subplot(1,length(tf_list),i)
    imagesc(flipud(log2(motif_mat(Idx,:))))
    title(curr_tf)
    xlabel(tf_motifs.(curr_tf));  
    caxis([2 8])
    set(gca,'xtick',[1,window,window*2-1],'XTickLabel',[-1*(window),0,window],'XTickLabelRotation',45)

end

cb = colorbar('position',[0.9180    0.1091    0.0029    0.8145]);
cb.Label.String  = 'Signal (log2)';
colormap(cm_YlGnBu)

%% Figure S1D  - plot a correlation matrix of sum ver promoter for all TFs

all_list = fieldnames(labWTs.norm);
tf_list = {'Gal11','Msn2','Msn4','Nfi1','Adr1','Mss11','Rsf2','Tda9','Asg1','Yap1', 'Skn7', 'Crz1','Gis1','Ixr1','Hsf1','Hot1',...
             'Cyc8','Nrg1','Nrg2','Sok2','Phd1','Mig1','Mig2','Mig3','Sko1','Dot6','Sut2','Sfl1','Rox1','Mot3','Yap4','Rgt1'};
         
other_TFs = all_list(~ismember(all_list,tf_list));
final_list = vertcat(tf_list',other_TFs);

% create a matrix of all sum signal data
all_sumProm = zeros(6701,length(final_list));
for i = 1:length(final_list)
    all_sumProm(:,i) = labWTs.sum_over_promoter.(final_list{i})';
end

figure;
subplot(1,2,1)
imagesc(corr(all_sumProm,'rows','pairwise'))
axis square
caxis([0 1])
colorbar
set(gcf,'color','w')
set(gca,'ytick',1:length(final_list),'yticklabel',final_list,'fontsize',6)
set(gca,'xtick',1:length(final_list),'xticklabel',final_list,'xticklabelrotation',45)

%create another matrix only for the selected subset
all_sumProm = zeros(6701,length(tf_list));
for i = 1:length(tf_list)
    all_sumProm(:,i) = labWTs.sum_over_promoter.(tf_list{i})';
end

subplot(1,2,2)
imagesc(corr(all_sumProm,'rows','pairwise'))
axis square
caxis([0 1])
colorbar
set(gcf,'color','w')
set(gca,'ytick',1:length(tf_list),'yticklabel',tf_list)
set(gca,'xtick',1:length(tf_list),'xticklabel',tf_list,'xticklabelrotation',45)

colormap(cm_YlGnBu)
%% Figure S1E - calculate for each tf it num of motifs in all promoters

tf_list = {'Msn2','Msn4','Gal11','Nfi1','Adr1','Asg1','Mss11','Rsf2','Tda9','Ixr1','Yap1','Crz1','Skn7','Gis1','Hsf1','Sko1','Hot1'};
motifs_tfs = {'Msn2','Rsf2','Yap1','Crz1','Skn7','Hsf1','Sko1'};
motif_climits = [7,7,3,7,12,6,3];


sumProm_all = zeros(6701,length(tf_list));
numProms = 100;

sumProm_list = [];
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

sumProm_mat = zeros(length(sumProm_list),length(tf_list));            
for i = 1:length(tf_list)
    curr_data = chec_struct_med.norm.(tf_list{i});
    curr_sumProm = zscore(chec_struct_med.sum_over_promoter.(tf_list{i}));
    sumProm_mat(:,i) = curr_sumProm(sumProm_list);
end
sumProm_mat = sumProm_mat';

[~,Idx] = sort(nanmean(sumProm_mat(:,1:numProms)),'descend');
 sumProm_mat(:,1:numProms) = sumProm_mat(:,Idx);

figure
a = subplot(14,1,1:6);
imagesc(sumProm_mat)
cb = colorbar; cb.Label.String = 'Z-score'; cb.Location = 'eastoutside';
set(gca,'ytick',1:length(tf_list),'yTickLabel',tf_list,'FontSize',10,'TickLabelInterpreter','none');
caxis([1 6])
set(gcf,'color','w')
colormap(a,cm_green)
hold on
set(gca,'xtick',[]);

motif_mat = zeros(length(motifs_tfs),length(sumProm_list));
for i = 1:length(motifs_tfs)
    curr_motif = tf_motifs.(motifs_tfs{i});
    
    % call for function to count number of motifs occurances in promoters
    num_motifs = numOfMotifs(sumProm_list,curr_motif,GP);
    motif_mat(i,:) = num_motifs';
    
    % plot
    e =subplot(14,1,i+7);
    imagesc(num_motifs');
    caxis([0 motif_climits(i)])
    colormap(e,cbrewer('seq','Greys',5))
    cb = colorbar ; cb.Location = 'eastoutside'; 
    cb.Label.String = curr_motif;
    if i <length(motifs_tfs)
         set(gca,'xtick',[])
    end
    set(gca,'ytick',1,'yTickLabel',[motifs_tfs{i}],'FontSize',8);  
    box off
end
xlabel('Promoters')

b =subplot(14,1,14);
imagesc(opnScore(sumProm_list)');
cb = colorbar; 
cb.Location = 'eastoutside';
set(gca,'ytick',1,'yTickLabel','OPN score');
set(gca,'xtick',[]);
colormap(b,flipud(cm_RdGy))





