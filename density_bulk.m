function y=density_bulk(alpha, gamma)

if alpha<0.5
    if gamma<=1-alpha
        y = alpha;
    else
        y = gamma;
    end
else
    if gamma<=0.5
        y = 0.5;
    else
        y = gamma;
    end
end

end