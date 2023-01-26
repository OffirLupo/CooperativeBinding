%% Figure6 supplementary
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


%% Figure S6A - calculate the length of the IDR for all TFs
T = readtable('Tables/tableS1.xlsx');
tf_list = {'Skn7','Ixr1','Sko1','Dot6','Sok2','Mot3','Mig3','Phd1','Mig2','Nrg1','Nrg2','Mig1','Rox1','Sfl1','Cin5','Rgt1','Sut2'};

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
%% Figure S6B - motif analysis
clear tempstruct 
tf_list = {'Nrg1','Nrg2','Mig1','Mig2','Mig3','Sok2','Phd1','Sko1','Sut2','Mot3','Rox1','Dot6'};

for i =1:length(tf_list)
    tempstruct.norm.(tf_list{i}) = labWTs.norm.(tf_list{i});
end

tempstruct = motifFind(tempstruct,7,GP,1:6701);

figure('Position',[2273         204         758         549])
[ha, pos] = tight_subplot(length(tf_list),2,[.01 .08],[.1 .1],[.1 .1]);
T = importdata('GeneralStructs\cisBP\TF_Information.txt');
T = cell2table(T.textdata(2:end,1:27),'VariableNames',T.textdata(1,1:27));
p = 1;

for i  = 1:length(tf_list)*2
    %plot 7mer seqlogo of the top 3 motifs for each TF
    axes(ha(p)); 
    seqlogoPlot(tempstruct,tf_list{i},'mer7',3);
    if p == 1 
        title('Data - 7mer','FontSize',6)
    end
    
    if i < length(tf_list)
     set(gca,'xtick',[])
    end
    p = p+1;
    
    axes(ha(p));
    
    %plot the seqlogo of the PWM from cisBP
    curr_tf_ind = find(strcmpi(T.TF_Name,tf_list{i}));
    if isempty(curr_tf_ind)
        p = p+1;
        continue
    end
    
    JAS_ind = strcmp(T.MSource_Identifier(curr_tf_ind),'JASPAR');
    PBM_ID = strcmp(T.Motif_Type(curr_tf_ind),'PBM');
    DeBoer11_ID = strcmp(T.MSource_Identifier(curr_tf_ind),'DeBoer11');
    
    if strcmp(tf_list{i},'Rgt1') 
        PBM_ID = 0;
    elseif strcmp(tf_list{i},'Yap1')
        PBM_ID = 0;
        JAS_ind = 0;
    end

    if  sum(PBM_ID)>0
        curr_tf_ind = curr_tf_ind(PBM_ID);
         curr_ID =  T.Motif_ID(curr_tf_ind);
         pwm = importdata(['GeneralStructs\\cisBP\pwms_all_motifs\',curr_ID{1},'.txt']);
         pwm = pwm.data(:,2:end)';
         pwm = pwm(:,find(median(pwm)<0.2));
         SeqLogoFig(pwm,'CUTOFF',0)
         ylim([0 2.5])
        % title('PBM','FontSize',6)
         text(1,2.2,'PBM')
         
    elseif sum(JAS_ind)>0
         curr_tf_ind = curr_tf_ind(JAS_ind);
         curr_ID =  T.Motif_ID(curr_tf_ind);
         pwm = importdata(['GeneralStructs\\cisBP\pwms_all_motifs\',curr_ID{1},'.txt']);
         SeqLogoFig(pwm.data(:,2:end)','CUTOFF',0)
          ylim([0 2.5])
         %title('JASPAR','FontSize',6)
         text(1,2.2,'JASPAR')
         
    elseif sum(DeBoer11_ID)>0
         curr_tf_ind = curr_tf_ind(DeBoer11_ID);
         curr_ID =  T.Motif_ID(curr_tf_ind(2));
         pwm = importdata(['GeneralStructs\\cisBP\pwms_all_motifs\',curr_ID{1},'.txt']);
         pwm =pwm.data(:,2:end)';
         % do reverse complement
         pwm = fliplr(pwm) ;
         pwm = flipud(pwm) ;

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

%% Degron - Figure S6C,D
% Plot heatmap top promoters of list of TFs
% Degron correlation matrix
tf_list1 = {'Cyc8','Ixr1','Ixr1_cyc8_deg_0','Ixr1_cyc8_deg_30','Ixr1_d16',...
            'Sok2','Sok2_cyc8_deg_0','Sok2_cyc8_deg_30','Sok2_d16',...
            'Sko1','Sko1_cyc8_deg_0','Sko1_cyc8_deg_30','Sko1_d16',...
            'Skn7','Skn7_cyc8_deg_0','Skn7_cyc8_deg_30','Skn7_d16',...
            'Mot3','Mot3_cyc8_deg_0','Mot3_cyc8_deg_30','Mot3_d16'};
all_combined_mat = zeros(6701,length(tf_list1));

for i =1:length(tf_list1)
    all_combined_mat(:,i) = chec_struct_med.sum_over_promoter.(tf_list1{i});
end

figure
imagesc(corr(all_combined_mat,'rows','pairwise'))
axis square
set(gcf,'color','w')
cb = colorbar;
cb.Label.String = 'Pearson r';
title('Promoter preference');
set(gca,'xtick',1:length(tf_list1),'XTickLabel',tf_list1,'XTickLabelRotation',45,'FontSize',12,'TickLabelInterpreter','none')
set(gca,'ytick',1:length(tf_list1),'yTickLabel',tf_list1)
colormap(cm_YlGnBu)
%% heatmap
clear sumProm_list tf_added sorted_proms sumprommat medBindMat

tf_list1 = {'Cyc8','Ixr1','Sok2','Sko1','Skn7','Mot3'};
tf_list2 = {'Cyc8','Ixr1_cyc8_deg_0','Sok2_cyc8_deg_0','Sko1_cyc8_deg_0','Skn7_cyc8_deg_0','Mot3_cyc8_deg_0'};
tf_list3 = {'Cyc8','Ixr1_cyc8_deg_30','Sok2_cyc8_deg_30','Sko1_cyc8_deg_30','Skn7_cyc8_deg_30','Mot3_cyc8_deg_30'};
tf_list4 = {'Cyc8','Ixr1_d16','Sok2_d16','Sko1_d16','Skn7_d16','Mot3_d16'};

     
tf_listProms = {'Cyc8','Ixr1','Sok2','Sko1','Skn7','Mot3'};

tf_list_all  = [{tf_list1}, {tf_list2}, {tf_list3},{tf_list4}];

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
sumProm_mat_del = zeros(length(tf_list3),6701);    
sumProm_mat_wt = zeros(length(tf_list3),6701);            

for i = 1:length(tf_list3)
    sumProm_mat_del(i,:) = zscore(chec_struct_med.sum_over_promoter.(tf_list2{i}));
    sumProm_mat_wt(i,:) = zscore(chec_struct_med.sum_over_promoter.(tf_list1{i}));
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
    subset1 = sumProm_mat_del(:,curr_proms);
    subset2 = sumProm_mat_wt(:,curr_proms);
    subset = subset1;% - subset2;
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
    %set(gca,'xtick',[]);
    medBindperTF = zeros(size(sumProm_mat,1),length(tf_added));
    count = 1;
    for p = 1:length(tf_added)
        curr_proms = count:count+length(tf_added{p}) -1;
        medBindperTF(:,p) = nanmedian(sumProm_mat(:,curr_proms),2);
        count = count+length(tf_added{p}) -1;
    end
    medBindMat{l} = medBindperTF;
end


%% Figure s6E,F - Msn2 in dcyc8

figure
xData = chec_struct_med.sum_over_promoter.Msn2;
yData = chec_struct_med.sum_over_promoter.Msn2_d12cycmod;
cData = chec_struct_med.sum_over_promoter.Cyc8;
scatter(xData,yData,30,cData,'filled','MarkerEdgeColor','k')
axis square
yl = ylim;
text(yl(1), yl(2) ,['R =  ' , num2str(round(corr(xData',yData'),2))])
cb = colorbar;
cb.Label.String = 'Cyc8';
set(gcf,'color','w')
set(gca,'fontsize',12)
xlabel('WT');
ylabel('Msn2 d12_cycModule');
title('Msn2')  
hold on
plot(xlim,xlim,'--k')
colormap(cm_YlOrBr)

tf_list = {'Nrg1','Msn2','Msn2_d12cycmod'};
plot_promoterSlim(chec_struct_med,MNase_seq,GP.gene_table.SUC2,tf_list,700,tf_motifs,GP)


%% Figure S6G - Mig1
tf_list = {'Mig1','Mig2','Mig3'};
size_list = {'Mig2','Mig3','Mig2'};
color_list = {'Mig3_d12','Mig1_d12','Mig1_d12'};
figure
for i = 1:length(tf_list)
    curr_tf = tf_list{i};
    xData = chec_struct_med.sum_over_promoter.(curr_tf);
    yData = chec_struct_med.sum_over_promoter.([curr_tf,'_d12']);

    cData = chec_struct_med.sum_over_promoter.(color_list{i});
    sData = chec_struct_med.sum_over_promoter.(size_list{i});
    
    [~,idx] = sort(sData,'descend');
    norm_to_max = sData./nanmedian(sData(idx(1:5)));
    sData = rescale(norm_to_max,10,100);
    subplot(1,length(tf_list),i)

    scatter(xData,yData,sData,cData,'filled')
    axis square
    yl = ylim;
    text(yl(1), yl(2) ,['R =  ' , num2str(round(corr(xData',yData'),2))])
    text(yl(1), yl(2)-100 ,['size =  ' , size_list{i}])
    cb = colorbar;
    cb.Label.String = color_list{i};
    set(gca,'fontsize',12)
    xlabel(curr_tf);
    ylabel([curr_tf,'_d12'],'Interpreter','none');
    gname(gene_names)
    hold on
    plot(ylim,ylim,'--k')

end
set(gcf,'color','w')
colormap(cm_YlOrBr)


