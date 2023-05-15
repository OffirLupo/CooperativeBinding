%this function load mat files of normalized read counts and merges them to
%one struct.

function outStruct = loadMatFiles(file_dir)
    file_list = dir(file_dir);
    for i = 1:length(file_list)
        if ~contains(file_list(i).name,'.mat')
            continue
       else
            curr_file = load([file_dir,'\',file_list(i).name]);
            curr_name = file_list(i).name(1:end-4);
            outStruct.norm.(curr_name) = curr_file.norm;   
       end
    end
end