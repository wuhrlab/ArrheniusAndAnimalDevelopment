function BICcomp = BICExpCompCalc(data,mdl,mdl2,stage)

    %load and parse data
        %constant
        n = data.numPoints(stage);
        yi = data.yi{stage,1};
        xi = data.xi{stage,1};
        %will change if poly degree is increased
        k = length(mdl);
        %will change if poly degree is increased
        k2 = length(mdl2);
        
    %AIC analysis on data
        esti = mdl.a*exp(mdl.b*xi);
        ResidualSS = (yi - esti).^2;
        SSR = sum(ResidualSS); %normalized for units somehow?
        esti2 = polyval(mdl2,xi);
        ResidualSS2 = (yi - esti2).^2;
        SSR2 = sum(ResidualSS2);
    
    %logL = n*(log(2*pi) + log(1/n) + log(SSR) + 1);% + (2*k^2 + 2*k)/(n-k-1);
    %logL2 = n*(log(2*pi) + log(1/n) + log(SSR2) + 1);
    %AIC = n*log(SSR/n) + k*log(n); %check if there is a way to correct for units so that ssr 
            %term does not overpower the correction term 
    %BIC 2 -> log(n)
    BICcomp = exp((log(n)*(k-k2) + n*log(SSR/SSR2))/2);

end