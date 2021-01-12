function curvature = calcCurvature(Ea, lnA)
%set constants 
R = 8.3144598*(10^(-3));
rangeT = [285.15 305.15];
T = rangeT(1) + diff(rangeT)/2;
%convert lnA -> A
A = exp(lnA);
%calculate curvature = abs(f''(x))/(1+(f'(x))^2)^[3/2]
%f''(x) = num1/denom1 - (num2.1+num2.2)/denom2
    %num1 = (sum((Ea.*exp(Ea./(R*T)))./((R*T^2).*A)))^2
        %scalar - runs
    %denom1 = (sum((exp(Ea./(R*T)))./(A)))^2
        %scalar - runs
    %num2.1 = sum((2.*Ea.*exp(Ea./(R*T)))./((R*T^3).*A))
        %scalar - runs
    %num2.2 = sum((Ea.^2.*exp(Ea./(R*T)))./((R^2*T^4).*A))
        %scalar -> this had the problem: Ea.^2 was Ea^2 - now fixed, runs
    %denom2 = sum((exp(Ea./(R*T)))./(A))
        %scalar - runs
%f'(x) = num3/denom3
    %num3 = sum((Ea.*exp(Ea./(R*T)))./((R*T^2).*A))
        %scalar - runs
    %denom3 = sum((exp(Ea./(R*T)))./(A))
        %scalar - runs
%Scalars = T, R
%Vectors = Ea, A
curvature = 0 - (abs((sum((Ea.*exp(Ea./(R*T)))./((R*T^2).*A)))^2/...
    (sum((exp(Ea./(R*T)))./(A)))^2 - (sum((2.*Ea.*exp(Ea./(R*T)))./...
    ((R*T^3).*A)) + sum((Ea.^2.*exp(Ea./(R*T)))./((R^2*T^4).*A)))/...
    sum((exp(Ea./(R*T)))./(A)))/(1+(sum((Ea.*exp(Ea./(R*T)))./...
    ((R*T^2).*A))/sum((exp(Ea./(R*T)))./(A)))^2)^(3/2));

end