function [common_promoters, unique_promoters] = diff_proms(checStruct,tf_1,tf_2,numProms,rank_treshold)

    [~,tf_1Idx] = sort(checStruct.sum_over_promoter.(tf_1),'descend');
    [~,tf_2Idx] = sort(checStruct.sum_over_promoter.(tf_2),'descend');

    n = 1;
    while length(intersect(tf_1Idx(1:n),tf_2Idx(1:n))) < numProms
            n = n+1;
    end
    common_promoters = intersect(tf_1Idx(1:n),tf_2Idx(1:n));

    n = 1;
    while length(setdiff(tf_2Idx(1:n),tf_1Idx(1:rank_treshold))) < numProms
        n = n+1;
    end
    unique_promoters = setdiff(tf_2Idx(1:n),tf_1Idx(1:rank_treshold),'stable');

end