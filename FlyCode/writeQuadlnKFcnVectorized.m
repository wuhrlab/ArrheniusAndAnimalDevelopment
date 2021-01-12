function f = writeQuadlnKFcnVectorized(Ea,lnA,rou)
%Constants, take out of ln space
R = 8.3144598*(10^(-3));
A = exp(lnA);
%write vectored function
syms f(T)
f(T) = -log(sum(exp((Ea./(R*T))+(rou./T^2))./A));
end