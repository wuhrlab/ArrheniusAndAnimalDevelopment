function k = calcKT(Ea,lnA)

    R = 8.3144598*(10^(-3));
    A = exp(lnA);
    rangeT = [285.15 305.15];
    T = rangeT(1) + diff(rangeT)/2;
    k = 1/(sum(exp(Ea./(R*T))./A));

end