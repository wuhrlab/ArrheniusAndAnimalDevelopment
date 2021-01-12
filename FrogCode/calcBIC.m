function BICval = calcBIC(data,mdl)

  %load and parse data
        %constant
        n = data.numPoints;
        yi = data.yi{1};
        xi = data.xi{1};
        %will change if poly degree is increased
        k = length(mdl);
        %will change if poly degree is increased
        
    %AIC analysis on data
        esti = polyval(mdl,xi);
        ResidualSS = (yi - esti).^2;
        SSR = sum(ResidualSS); %normalized for units somehow?

    BICval = k*log(n)+n*log(SSR/n); 
        %BIC = k*ln(n)+n*ln(RSS/n), 

end