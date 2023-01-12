function checStruct_new = motifFind(checStruct,xMer,GP,genes)
% motif analysis

%%
minReads = 20; % under this number signal will be 0; 
slidingWindowSize = 30; % for smoothing the signal
promoterLength = 500; % the length of upstream bp from TSS to consider as promoters

%%

    refGenome = fastaread('Genral structs\S288C_reference_sequence_R64-1-1_20110203.fsa');
    tss = GP.stein_tss13; % 6701 x 3. rows: genes, columns: chromosome,start,end
    
%% 
% convert genome into 0-3
% A=0 T=1 G=2 C=3
for i = 1:16
    numberGenome{i} = zeros(1,length(refGenome(i).Sequence));
end

for i = 1:16
    numberGenome{i}(strfind(refGenome(i).Sequence,'A')) = 0;
    numberGenome{i}(strfind(refGenome(i).Sequence,'T')) = 1;
    numberGenome{i}(strfind(refGenome(i).Sequence,'G')) = 2;
    numberGenome{i}(strfind(refGenome(i).Sequence,'C')) = 3;
end


%% maping the genome sequence to x-mers  

genomeXmers = [];

xmerBits = 4.^[0:(xMer-1)];

for i = 1:16
    chrXmers =conv(numberGenome{i}, xmerBits);
    chrXmers(1:xMer-1) = [];
    genomeXmers =[genomeXmers,chrXmers];
end
% can't handele zeros xmers
genomeXmers(genomeXmers==0) = 4.^xMer;


checStruct.motifs.(['mer',num2str(xMer)]).genomeXmers = genomeXmers;

%% generating a matrix of all checSeq TFs, rows: TFs, columns: genome positions (chromosomes concatenated)
my_samples = fieldnames(checStruct.norm);

checSeqSignal = zeros(length(my_samples),length(genomeXmers));

for tf = 1:length(my_samples)
    curr_tf = char(my_samples(tf))
    tmp = [];
    for i = 1:16
        a = size(checStruct.norm.(curr_tf){i});
        if a(1) > 1
           tmp = [tmp ,checStruct.norm.(curr_tf){i}'];
        else
           tmp = [tmp ,checStruct.norm.(curr_tf){i}]; 
        end       
    end
    checSeqSignal(tf,:) = tmp; 
end

%checStruct.checSeqSignal = checSeqSignal;


%% smoothing the signal
curr_window = ones(slidingWindowSize,1);
edge = slidingWindowSize-floor(slidingWindowSize/2)-1-round(xMer/2);
checSeq_smooth = nan(1,length(checSeqSignal)+edge);


for tf = 1:length(my_samples)
    tmp = checSeqSignal(tf,:);
    tmp(tmp < minReads) = 0;
    tf_smooth = conv(tmp ,curr_window);
    tf_smooth(1:slidingWindowSize-1-edge) = [];
    checSeq_smooth(tf,:) = tf_smooth;
    
end

%checStruct.checSeq_smooth = checSeq_smooth;
    
%% Annotating promoters (logical)

promotersRegion = [];
tmp = []; %temporary promoters regions

for i = 1:16
    tmp{i} = zeros(1,length(refGenome(i).Sequence));
end

genes_with_annotated_stanscript = find(~isnan(tss(:,1)));
genes_with_annotated_stanscript = intersect(genes_with_annotated_stanscript,genes,'stable');
    
for i = 1:length(genes_with_annotated_stanscript)
    gene = genes_with_annotated_stanscript(i);
    locations = tss(gene,:);
    strand = sign(locations(3)-locations(2));
    tmp{locations(1)}((locations(2)-strand*promoterLength):strand:locations(2)) = 1;
end

for i = 1:16
    promotersRegion = [promotersRegion,tmp{i}];
end

promotersRegion = logical(promotersRegion);

checStruct.motifs.promotersRegion = promotersRegion;

%% Saving all possible x-mers in xmerSeq 
xmerSeq = [];
num2nuc{1} = 'A';
num2nuc{2} = 'T';
num2nuc{3} = 'G';
num2nuc{4} = 'C';

for i = 1:(4.^xMer-1)
    curr_motif_binary = fliplr(de2bi(i,xMer,4));
    curr_motif_seq = [];
    for n = 1:xMer
        curr_nuc_binary = curr_motif_binary(n);
        curr_motif_seq = [curr_motif_seq,num2nuc{curr_nuc_binary+1}];
    end
    xmerSeq{i} = curr_motif_seq;
end

xmerSeq{4.^xMer} = repmat('A',1,xMer);


%% For each motif find it's reverse complement motif, also combine forward & reverse into one representative motif
reverseMotif = [];
combineFR = [];
palindroms = [];
for i = 1:4^xMer
    if mod(i,100) ==0
        fprintf('in %d...\n',i);
    end
    curr_motif_seq = xmerSeq{i};
    curr_RevComp = seqrcomplement(curr_motif_seq);
    reverseMotif(i) = find(strcmp(xmerSeq,curr_RevComp));
    if i == reverseMotif(i)
        palindroms = [palindroms ,i];
    end
end


%% save to chec struct
checStruct.motifs.(['mer',num2str(xMer)]).palindroms = palindroms; %without half sites
checStruct.motifs.(['mer',num2str(xMer)]).motifs_seq = xmerSeq;
checStruct.motifs.(['mer',num2str(xMer)]).reverse_motif = reverseMotif;


%% calculate mean signal for each motif in the promoters regions

my_samples = fieldnames(checStruct.norm);

for i = 1:length(my_samples)
    curr_tf = char(my_samples(i));
    fprintf('working on %s...\n',curr_tf);
    tmp_meanSignalPerMotif = accumarray(genomeXmers(promotersRegion)', checSeq_smooth(i,promotersRegion)',[],@mean);
    tmp_meanSignalPerMotif(end+1:4^xMer) = 0;
 
    meanSignalPerMotif = zeros(size(tmp_meanSignalPerMotif));
    for j = 1:4.^xMer
        rev = reverseMotif(j);
        if j ~= rev
            meanSignalPerMotif(min(rev,j)) = tmp_meanSignalPerMotif(j)+ tmp_meanSignalPerMotif(rev);
        else
            meanSignalPerMotif(j) = tmp_meanSignalPerMotif(j);
        end
    end
    checStruct.motifs.(['mer',num2str(xMer)]).meanSignalPerMotif.(curr_tf) = meanSignalPerMotif;
end

checStruct_new = checStruct;

disp('Finished!');
end