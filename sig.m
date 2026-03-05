function value = sig(ss)
    value = 0;
    if ss < 0
        value = 1;
    elseif ss > 0
        value = 2;
    end
end