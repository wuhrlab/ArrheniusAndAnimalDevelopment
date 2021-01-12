function BICcomp = fitBICCompCalc(data,mdl,mdl2,stage)

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
        esti = polyval(mdl,xi);
        ResidualSS = (yi - esti).^2;
        SSR = sum(ResidualSS); %normalized for units somehow?
        esti2 = polyval(mdl2,xi);
        ResidualSS2 = (yi - esti2).^2;
        SSR2 = sum(ResidualSS2); %normalized for units somehow?
        %correction for 'small' n, improve a tiny bit
        c1 = (2*k^2 + 2*k)/(n - k - 1);
        c2 = (2*k2^2 + 2*k2)/(n - k2 - 1);
    %logL = n*(log(2*pi) + log(1/n) + log(SSR) + 1);% + (2*k^2 + 2*k)/(n-k-1);
    %logL2 = n*(log(2*pi) + log(1/n) + log(SSR2) + 1);
    %AIC = n*log(SSR/n) + k*log(n); %check if there is a way to correct for units so that ssr 
            %term does not overpower the correction term 
    %BIC 2 -> log(n)
    %AICcomp = exp((2*(k-k2) - n*log(SSR2/SSR) + (c1-c2))/2);%AIC with
        %correction
    BICcomp = exp((log(n)*(k-k2) - n*log(SSR2/SSR))/2); 
        %BIC drastically improve, 
        %but still seems very off, cuts every roughly in hald
    %deltaBIC = log(n)*(k-k2) - n*log(SSR2/SSR);
end
%AICcomp = exp((minAIC-maxAIC)/2)
%AIC = 2k - 2ln(Lhat)
        %k = cum estimated parameters
        %Lhat = maximum value of likelyhood function
            %= -(n/2)ln(2pi) - (n/2)ln(1/n) -
            %for linear models (in the parameters) as follows
            %(n/2)ln*sum(yi-b0+b1xi+bnxi^n)^2-(n/2)
                %sum(yi-b0+b1xi+bnxi^n)^2 = sum of residual squares
            %simplified = -(n/2)(log(2pi) + log(1/n) + log(SSR) + 1)
        %AIC = 2k - 2(-(n/2)(log(2pi) + log(1/n) + log(SSR) + 1))
            %= 2k + n(log(2pi) + log(1/n) + log(SSR) + 1)
        %minimum = best fit
        
%Small data set <=40
    %assuming univariate and linear in parameters
    %AICc = AIC + (2k^2 + 2k)/(n - k -1)
    
%Comparing Fits (relative likelyhood)
    %minimum = best fit
    %prob min loss = exp((AIC_min - AIC_i)/2)
        %compares the models with the minimum (best model) to determine how
        %probable (as times the min model) to minimize info loss. Ex. 0.368
        %means model 2 is 0.368 times as probable to minimize loss as the
        %minimum/best model (thus less likely to be the "true" model)
