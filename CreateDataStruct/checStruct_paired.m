%%%% This function takes path of out files, normalize for read number, and
%%%% create a mat struct for further processing

function [samplesStruct] = checStruct_paired(path,checOrMNase)

GP = load('CreateDataStruct\general_params_130711.mat');

chr_length = GP.chr_len;
%chr_length(17)  = 85779; %This is mitochondria
folders = dir(path);

for i = 1:length(folders)
    if  ~isempty(strfind(folders(i).name,'out')) == 1
        fileFullPath = ([path,'/',folders(i).name]);
        disp(['now working on ' fileFullPath]);
        currName = strrep(folders(i).name,'.out','');
        disp(currName)
        currName = currName(1:end-1);
        currName = regexprep(currName,'-','_');
        
        curr_data = importdata(fileFullPath);       
        readsByChrom = zeros(1,16); 
        
        % normalizing the data to 10000000 reads        
        c=1;
        for j=1:16            
            profile{j} = curr_data(c:c+chr_length(j)-1)';
            c = c+chr_length(j);
            readsByChrom(j) = sum(profile{j});
        end
            %samplesStruct.mitoReads.(currName) = sum(readsByChrom(17))/2;
            samplesStruct.totalReads.(currName) =  sum(readsByChrom(1:16))/2;
            
        for j=1:16
            if strcmp(checOrMNase,'chec')
                 samplesStruct.norm.(currName){1,j} =  (profile{j} ./ samplesStruct.totalReads.(currName)) .* 5000000;
            elseif strcmp(checOrMNase,'MNase')
                samplesStruct.norm.(currName){1,j} =  (profile{j} ./ samplesStruct.totalReads.(currName)) .* 5000000 *100;
            end
        end
        
    end
end
end