function y_new = smooth_offir(y,n)
    norm = ones(1,length(y));
    norm(isnan(y)) = 0;
    y(isnan(y)) = 0;
    y_new = filter2(ones(1,n),y) ./ filter2(ones(1,n),norm);
end