function pVal = ancova(regAn)

    %load and parse data
    numPoints = regAn.numPoints;
    SX = regAn.SX;
    SY = regAn.SY;
    SXX = regAn.SXX;
    SYY = regAn.SYY;
    SXY = regAn.SXY;
    AveX = regAn.AveX;
    AveY = regAn.AveY;
    
    %Ancova analysis on data
    df = numPoints - 2;
    SdX2 = SXX-SX.^2./numPoints;
    SdXdY = SXY-(SX.*SY)./numPoints;
    SdY2 = SYY-SY.^2./numPoints;
    b = SdX2./SdXdY;
    ResidualSS = SdY2-((SdXdY.^2)./SdX2);
    MS = ResidualSS./df;
    pDf = sum(df);
    pResidualSS = sum(ResidualSS);
    pMS = pResidualSS/pDf;
    regDf = length(numPoints); %number of regressions
    cDf = pDf + regDf;
    cSdX2 = sum(SdX2);
    cSdXdY = sum(SdXdY);
    cSdY2 = sum(SdY2);
    cResidualSS = cSdY2-((cSdXdY^2)/cSdX2);
    cMS = cResidualSS/cDf;
    regResidualSS = cResidualSS - pResidualSS;
    regMS = regResidualSS/regDf;
    F = regMS/pMS;
    numDf = regDf;
    denomDf = pDf;
    %slopes = p>0.05, slopes ~= p<0.05
    pVal = 1 - fcdf(F,numDf,denomDf);

end