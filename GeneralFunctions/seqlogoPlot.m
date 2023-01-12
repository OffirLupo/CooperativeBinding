function seqlogoPlot(checStruct,tf,mer,num_of_motifs)

sample = checStruct.motifs.(mer).meanSignalPerMotif.(tf);

[meanSignal_sorted, meanSignal_sorted_idx] =sort(sample, 'descend');
meanSignal_sorted=meanSignal_sorted(1:num_of_motifs);
meanSignal_sorted_idx=meanSignal_sorted_idx(1:num_of_motifs);

MotSet = checStruct.motifs.(mer).motifs_seq(meanSignal_sorted_idx);

flip_tfs = {'Msn4','Rsf2','Tda9','Yap1','Gis1','Asg1','Crz1','Ixr1','Mig1','Mot3','Nrg1','Nrg2','Sut2','Rgt1'};
if  ismember(tf,flip_tfs)
    MotSet{1} = seqrcomplement(MotSet{1});
end

meanSignal_sorted  = round(meanSignal_sorted);
motRank  = [];
for i = 1:length(MotSet)
    motRank (i) = round((meanSignal_sorted(i) / sum(meanSignal_sorted)),2) * 100;
end

motRank = int64(motRank);

for i = 2:length(MotSet)
    [~,forAlign] = nwalign(MotSet{1},MotSet{i},'glocal',true);
    forAlignMatch = sum(forAlign(2,:)=='|');
    [~,revAlign] = nwalign(MotSet{1},seqrcomplement(MotSet{i}),'glocal',true);
    revAlignMatch = sum(revAlign(2,:)=='|');
    if revAlignMatch > forAlignMatch
        MotSet{i} = seqrcomplement(MotSet{i});
    end
end

finalMatforAlign = cell(sum(motRank)*2,1);
j=0;
for i = 1:2:sum(motRank )*2-1
    j=j+1;
    finalMatforAlign{i} = ['>',num2str(j)];
end
e = 0;
for f = 1:length(motRank )
    currSize = motRank (f);
    for i = 1:currSize
        finalMatforAlign{int64((e+i)*2)} = MotSet{f};
    end
    e = e + motRank (f);
end

formatAlign = fastaread(char(finalMatforAlign));
finalAlign = multialign(formatAlign,'terminalGapAdjust',true);

sumA = [];
sumG = [];
sumC = [];
sumT = [];
totAGCT = [];

for h= 1:length(finalAlign(1).Sequence)
    currPos = char();
    for i = 1:length(finalAlign)
        currPos = [currPos,char(finalAlign(i).Sequence(h))];
        sumA(h) = sum(currPos=='A');
        sumG(h) = sum(currPos=='G');
        sumC(h) = sum(currPos=='C');
        sumT(h) = sum(currPos=='T');
        totAGCT(h) = length(regexp(currPos,'[AGCT]'));
    end
end

[~,iMax] = max(smooth(totAGCT));
if  contains(mer,'mer6')
 Positions2Plot = [iMax-2:iMax+2];
elseif contains(mer,'mer5')
 Positions2Plot = [iMax-2:iMax+2];
elseif contains(mer,'mer7')
     Positions2Plot = [iMax-3:iMax+3];
end

sumAfreq = sumA ./ totAGCT;
sumGfreq = sumG ./ totAGCT;
sumCfreq = sumC ./ totAGCT;
sumTfreq = sumT ./ totAGCT;
finalFreq = [sumAfreq;sumCfreq;sumGfreq;sumTfreq];

SeqLogoFig(finalFreq(:,Positions2Plot(1):Positions2Plot(end)),'CUTOFF',0)
ylim([0 2])

set(gcf,'color','w')
ylabel(tf,'Interpreter','none')

end



