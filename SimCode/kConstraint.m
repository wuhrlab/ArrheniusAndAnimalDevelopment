function [c, ceq] = kConstraint(Ea,lnA)

    %c(x) < 0

    %limiting individuak ks
    c(:,1) = calcKT(Ea(1),lnA(1)) - 1; %less than or eq to 1, use minus
    c(:,2) = -calcKT(Ea(1),lnA(1)) + 1/240000; %more than or eq to 1/240000, use +
    c(:,3) = calcKT(Ea(2),lnA(2)) - 1;
    c(:,4) = -calcKT(Ea(2),lnA(2)) + 1/240000;
    
    %limiting overal k
    c(:,5) = calcKT(Ea,lnA) - 1/60;
    c(:,6) = -calcKT(Ea,lnA) + 1/300000;
    
    ceq = [];
    
end