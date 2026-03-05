function chiN = side(gr,gfov1,gfov2,n)
    ss = zeros(1,n);
    for i = 1:n
        ss(i) = sig(gr(:,i)'*(gfov2-gfov1));
    end
    cs0 = countj(0,ss);
    cs1 = countj(1,ss);
    cs2 = countj(2,ss);
    chiN = -1;
    if 2 <= max(cs1,cs2) && max(cs1,cs2) == n-cs0 && n-cs0 <= n
        chiN = 1;
    elseif n-1 <= cs0 && cs0 <= n
        chiN = 0;
    end
end