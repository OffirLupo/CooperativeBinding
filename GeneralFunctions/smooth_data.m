function [smoothed_data] = smooth_data(data_mat, box_size)

    smoothed_data = zeros(size(data_mat));
    for i= 1: size(data_mat,1)
        smoothed_data(i,:) = smooth_offir(data_mat(i,:),box_size);
    end
end