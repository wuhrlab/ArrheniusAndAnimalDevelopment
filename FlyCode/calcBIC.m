function BICval = calcBIC(data,mdl,S,mu)

  %load and parse data
        %constant
        n = data.numPoints;
        yi = data.yi{1};
        xi = data.xi{1};
        %will change if poly degree is increased
        k = length(mdl);
        %will change if poly degree is increased
        
    %AIC analysis on data
        [esti,delta] = polyval(mdl,xi,S,mu);
        ResidualSS = (yi - esti).^2;
        SSR = sum(ResidualSS); %normalized for units somehow?

    BICval = k*log(n)+n*log(SSR/n); 
    BICval = k*log(n)+n*log(SSR/n);%k*2+n*log(SSR/n)+(2*k^2+2*k)/(n-k-1);
        %BIC = k*ln(n)+n*ln(RSS/n), 

end