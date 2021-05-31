%%Key
%Objective of the following code sub-chunk
for example = 1:9
%Objective of following code sub-chunk
    code here
        %any specific notes about function or specific coded values
end

%Each time course (with multiple embryos) is a single batch, take mean and
%treat as a single replicate
%%
%Data Import

%Destination
    filename = 'FlyScoresFinalCombined.xlsx';
    folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Fly Scores/';

%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{1});

%find number of events scored ("scores") in file
    isSearch = cellfun(@(x)isequal(x,'Score Number'),data);
    [row,col] = find(isSearch);
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+2;
        %+1 to include T0 as a score
    
%Initialize arrays to store different data types        
    accumulated = cell(sizeData,1);
    perStage = cell(sizeData,1);
    Temps = zeros(sizeData,1);
    stageNames = strings(1,numScores);
    stageAbbs = strings(1,numScores);

%iterate over each sheet to parse data types into their respective arrays        
for i = 1:sizeData
%read in specific data of sheet i
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{i});
    
%parse out "accumulated" data rows for each embryo of sheet i and store in 
%cell array
    isSearch = cellfun(@(x)isequal(x,'Accumulated'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize cell i to contain each accum data row from a single sheet
        accumulated{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
        for embryo = 1:length(row)
        %iterate over all of the columns (or stages of interest) from sheet
            for column = 1:numScores
             %check if data is not missing, if so, copy over into matrix
                if ismissing(data{row(embryo),(col+1+column)})==0
                    accumulated{i}(embryo,column) = cell2mat(data(row...
                        (embryo),(col+1+column)))*60; %convert to seconds
                else
                end
            end
        end
        
%parse out per stage data for each sheet in same manner
    isSearch = cellfun(@(x)isequal(x,'Per stage'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize a cell to contain each perStage data row from a single sheet
        perStage{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
    for embryo = 1:length(row)
    %iterate over all of the columns (or stages of interest) from sheet
        for column = 1:numScores
        %check is data is not missing, if so, copy over into matrix
            if ismissing(data{row(embryo),(col+1+column)})==0
                perStage{i}(embryo,column) = cell2mat(data(row...
                    (embryo),(col+1+column)))*60; %convert to seconds
            else
            end
        end
    end

%parse out temperature data for each sheet and store in array
    isSearch = cellfun(@(x)isequal(x,'Temp'),data);
    [row,col] = find(isSearch);
    Temps(i) = cell2mat(data(row,(col+1)));
end

%parse out the full names of every score, including T0
    isSearch = cellfun(@(x)isequal(x,'Score Name'),data);
    [row,col] = find(isSearch);
    for column = 1:numScores
        stageNames(column) = string(data(row,(col+1+column)));
    end
    
%parse out the abbreviations of every score
    isSearch = cellfun(@(x)isequal(x,'Abb.'),data);
    [row,col] = find(isSearch);
    for column = 1:numScores
        stageAbbs(column) = string(data(row,(col+1+column)));
    end

%%
%falahati vs our data
obsTemps = 1:length(Temps);
    meanAccByStage = zeros(sum(length(obsTemps)),length(obsScores));
    err = zeros(sum(length(obsTemps)),length(obsScores));
for temp = 1:length(obsTemps)
    %parse out score data for specific temp
        data = accumulated{obsTemps(temp)}(:,obsScores)/60;
    %determine number of measurable samples for each score
        N = sum(~isnan(data),1);
        meanAccByStage(temp,:) = mean(data,1,'omitnan');
        yStdErr = std(data,1,'omitnan');
        err(temp,:) = yStdErr;
        errTemp = 0.5;
        xneg(temp,:) = Temps(obsTemps(temp))-(Temps(obsTemps(temp))+errTemp);
        xpos = xneg;
end
yyaxis left
cr = scatter(Temps(1:12)-273.15,abs(meanAccByStage(1:12,1)),'LineWidth',1.5);
hold on 

crQ = fit(Temps(1:12)-273.15,abs(meanAccByStage(1:12,1)),'exp1');
%crVal = polyval(crQ,Temps-273.15);
%crFit = plot(Temps-273.15,crVal);
crFit = plot(crQ,Temps(1:12)-273.15,abs(meanAccByStage(1:12,1)))

fTimePixScale = [0,20,40,60,80,100;510,414,317,221,124,27];
fTimePix = [72,307,399,457,451];
fTiming = 100 - ((fTimePix-27)/(510-27))*100;
fTempPixScale = [5,10,15,20,25,30,35;76,129,182,234,287,340,393];
fTempPix = [93,137,180,256,393];
fTemps = 5 + ((fTempPix-76)/(393-76))*30;

fl = scatter(fTemps,fTiming,'mo','LineWidth',1.5);

fQ = fit(fTemps.',fTiming.','exp1');
%fVal = polyval(fQ,Temps-273.15);
%fFit = plot(Temps-273.15,fVal);
fFit = plot(fQ,Temps(1:12)-273.15,abs(meanAccByStage(1:12,1)));

title({"Falahati Timing vs Crapse Timing", "Syncitial Cleavage"})
ylabel("Timing (min)")
xlabel("Temperature (^oC)")
set(gca, 'LineWidth',1.5,'FontSize',16)

yyaxis right
ratio = crQ(Temps(1:12)-273.15)./fQ(Temps(1:12)-273.15);
rplot = plot(Temps(1:12)-273.15,ratio,'LineWidth',1.5)
ylim([0 5])
ylabel("Ratio Crapse/Falahati Timing")

legend([cr fl rplot],{"Crapse (NC14)","Falahati (NC11)","Ratio"})


%%
%Figure 1 C, non-log version
%as above, but using a personalized errorbar function to allow for
%transparency
%plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
%Raw Data Presentation

%generate as seperate figures

%array of temperatures to show (only being selective when there
%are temperatures very close together)
    obsTemps = 1:length(Temps);

%Figure 1, all data
    figureWidth = 1600;
    figureHeight = 675;
    figure('Position',[100 100 figureWidth figureHeight])
    aspectRatio = figureWidth/figureHeight;
%select what scores to look at, where 1 is T0
    obsScores = (1:numScores); 
    colors = distinguishable_colors(numScores);
    RGBa = horzcat(colors,zeros(numScores,1));
    RGBa = [RGBa(numScores-1,:);RGBa(1:numScores-2,:);RGBa(numScores,:)];
   
%initialize display data arrays for means and errors
    meanAccByStage = zeros(sum(length(obsTemps)),length(obsScores));
    err = zeros(sum(length(obsTemps)),length(obsScores));

%iterate over every observable temperature to condense multiple embryos in
%each score to a mean and calculate errors
    xneg = nan(length(obsTemps),1);
for temp = 1:length(obsTemps)
    %parse out score data for specific temp
        data = accumulated{obsTemps(temp)}(:,obsScores)/60;
    %determine number of measurable samples for each score
        N = sum(~isnan(data),1);
        meanAccByStage(temp,:) = mean(data,1,'omitnan');
        yStdErr = std(data,1,'omitnan');
        err(temp,:) = yStdErr;
        errTemp = 0.5;
        xneg(temp,:) = Temps(obsTemps(temp))-(Temps(obsTemps(temp))+errTemp);
        xpos = xneg;
end


%iterate over every score to to plot means of every temp
    ylimRange = [(min(Temps)-1.5) (max(Temps)+1.5)];
    xlimRange = [-250 6500];
    LineWidth = 2.5;
    MarkerSize = 50;
    CapSize = 0.003;
for score = 1:length(obsScores)
    hold on;
%plot basic for T0, specific with errors for T1 on
    if score == 2
        RGBa(obsScores(score),4) = 1;
        x = zeros(length(obsTemps),1)';
        y = Temps(obsTemps)';
        xError = zeros(1,length(x));
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x,y,RGBa(obsScores(score),:),xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
    elseif sum(score == [2 8 12]) ~= 1
    %check to see is any NaNs, index by not NaN. Make temporary variable to
    %index through
        idx = ~any(isnan(meanAccByStage(:,score)),2);
        devTime = meanAccByStage(:,score);
        errTime = err(:,score);
        t = Temps(obsTemps);  
        RGBa(obsScores(score),4) = 0.8;
        x = devTime';
        y = t';
        xError = errTime';
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,...
        %ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),...
            yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,...
            aspectRatio)
    else
    %check to see is any NaNs, index by not NaN. Make temporary variable to
    %index through
        idx = ~any(isnan(meanAccByStage(:,score)),2);
        devTime = meanAccByStage(:,score);
        errTime = err(:,score);
        t = Temps(obsTemps);  
        RGBa(obsScores(score),4) = 1;
        x = devTime';
        y = t';
        xError = errTime';
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,...
        %ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),...
            yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,...
            aspectRatio)
    end
end
%Aestetics
    %title('Fly stages developmental times over range of temperatures')
    ylim(ylimRange);
    ylabels = str2num(cell2mat(yticklabels));
    yticklabels(cellstr(num2str((ylabels+0.15)-273.15)));
    ylabel('Temperature (^oC)');
    xlim(xlimRange)
    xlabel('Time post 14^{th} Cleavage (min)');
    displayNames = stageAbbs(obsScores);
    %dummy figure to assign colors in legend
        for i = 1:numScores
            sp(i) = plot(1,1,'Color',RGBa(i,:),'LineWidth',LineWidth);
        end
        set(gca,'LineWidth', 3, 'FontSize', 35);
        legend(sp,cellstr(displayNames),'NumColumns',2, 'FontSize',35);
    
    
%%
%Appendix Figure 2, log version
%as above, but using a personalized errorbar function to allow for
%transparency
%plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
%Raw Data Presentation

%generate as seperate figures

%array of temperatures to show (only being selective when there
%are temperatures very close together)
    obsTemps = 1:length(Temps);

%Figure 1, all data
    figureWidth = 1600;
    figureHeight = 675;
    figure('Position',[100 100 figureWidth figureHeight])
    aspectRatio = figureWidth/figureHeight;
%select what scores to look at, where 1 is T0
    obsScores = (1:numScores); 
    colors = distinguishable_colors(numScores);
    RGBa = horzcat(colors,zeros(numScores,1));
    RGBa = [RGBa(numScores-1,:);RGBa(1:numScores-2,:);RGBa(numScores,:)];
   
%initialize display data arrays for means and errors
    meanAccByStage = zeros(sum(length(obsTemps)),length(obsScores));
    err = zeros(sum(length(obsTemps)),length(obsScores));

%iterate over every observable temperature to condense multiple embryos in
%each score to a mean and calculate errors
    xneg = nan(length(obsTemps),1);
for temp = 1:length(obsTemps)
    %parse out score data for specific temp
        data = accumulated{obsTemps(temp)}(:,obsScores)/60;
        data(:,3:numScores) = log(accumulated{obsTemps(temp)}(:,3:numScores)/60);
        data(:,1) = -log(abs(accumulated{obsTemps(temp)}(:,1)/60));
    %determine number of measurable samples for each score
        N = sum(~isnan(data),1);
    %take mean of each score (no NaNs) then calculate 95% CI and store
        meanAccByStage(temp,:) = mean(data,1,'omitnan');
        yStdErr = std(data,1,'omitnan');
        err(temp,:) = yStdErr;
        errTemp = 0.5;
        xneg(temp,:) = Temps(obsTemps(temp))-(Temps(obsTemps(temp))+errTemp);
        xpos = xneg;
end


%iterate over every score to to plot means of every temp
    ylimRange = [(min(Temps)-1.5) (max(Temps)+1.5)];
    xlimRange = [-log(120) log(6500)];
    LineWidth = 2.5;
    MarkerSize = 50;
    CapSize = 0.003;
for score = 1:length(obsScores)
    hold on;
%plot basic for T0, specific with errors for T1 on
    if score == 2
        RGBa(obsScores(score),4) = 1;
        x = zeros(length(obsTemps),1)';
        y = Temps(obsTemps)';
        xError = zeros(1,length(x));
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,...
        %ylimRange,xlimRange)
        plotPersonal(x,y,RGBa(obsScores(score),:),xError,yError,LineWidth,...
            MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
    elseif sum(score == [2 8 12]) ~= 1
    %check to see is any NaNs, index by not NaN. Make temporary variable to
    %index through
        idx = ~any(isnan(meanAccByStage(:,score)),2);
        devTime = meanAccByStage(:,score);
        errTime = err(:,score);
        t = Temps(obsTemps);  
        RGBa(obsScores(score),4) = 0.80;
        x = devTime';
        y = t';
        xError = errTime';
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,...
        %ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),...
            yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,...
            aspectRatio)
    else
    %check to see is any NaNs, index by not NaN. Make temporary variable to
    %index through
        idx = ~any(isnan(meanAccByStage(:,score)),2);
        devTime = meanAccByStage(:,score);
        errTime = err(:,score);
        t = Temps(obsTemps);  
        RGBa(obsScores(score),4) = 1;
        x = devTime';
        y = t';
        xError = errTime';
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,...
        %ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),...
            yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,...
            aspectRatio)
    end
end
%Aestetics
    %title('Fly stages developmental times over range of temperatures')
    ylim(ylimRange);
    ylabels = str2num(cell2mat(yticklabels));
    yticklabels(cellstr(num2str((ylabels+0.15)-273.15)));
    ylabel('Temperature (^oC)');
    xlabel('ln(time) post 14^{th} Cleavage [ln(min)]');
    displayNames = stageAbbs(obsScores);
    %dummy figure to assign colors in legend
        for i = 1:numScores
            sp(i) = plot(1,1,'Color',RGBa(i,:),'LineWidth',LineWidth);
        end
        set(gca,'LineWidth', 3, 'FontSize', 35);
        legend(sp,cellstr(displayNames),'NumColumns',2, 'FontSize',35);
     
%%
%Selected in depth CV Analyis
%Data Import (rerun actual data after this because this script overides the
%good data loaded in before)

figure('Position',[10 10 800 600]);

meanPerStage = cell2mat(cellfun(@(x)mean(x,1,'omitnan'),perStage,'UniformOutput',false));
stdPerStage = cell2mat(cellfun(@(x)std(x,1,'omitnan'),perStage,'UniformOutput',false));
cvMat = (stdPerStage./meanPerStage)*100;
cvMatNoT0 = cvMat(:,[1,3:12]);
plot(Temps-273.15,cvMatNoT0,'LineWidth',1.5)
title("Temporal Variance Over Temperatures in Flies")
xlabel("Temperature ^oC")
ylabel("Coefficent of Variance (CV%)")


labs = strcat(stageAbbs(1:11)," - ", stageAbbs(2:12));
legend(labs,'NumColumns',4,'Location','northwest')
set(gca,'FontSize', 25,'LineWidth',2);

%%
%Figure S3C
%Arrhenius Stage-wise
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    
    embryoMass = 8.99*10^-3;% mg

%Plot Stages
    figure('Position',[10 10 1350 850]);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    obsScores = (1:numScores);
    perStageInterval = cellfun(@(x) x(:,obsScores),perStage,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 sizeData-1 sizeData];
    viableTemps = [2:sizeData-1];
    
%Mass correction by t/m^(1/4)
    %perStageInterval = cellfun(@(x) x(:,obsScores)./embryoMass^0.25,perStage,'UniformOutput',false);    
    
    
%Initialize tables
    %Statistical Table
        stagRegM = zeros(length(obsScores),8);
        stagRegT = array2table(stagRegM);
        stagRegT.Properties.VariableNames = {'numPoints','SX','SY',...
            'SXX','SYY','SXY','AveX','AveY'};
        stagRegT.Properties.RowNames = stageNamesNoT0;

%convert to pseudo Ks and inverse temps
    pseudoKperStage = cell(length(perStageInterval),1);
    for temp = 1:length(perStageInterval)
        pseudoKperStage{temp} = 1./abs(perStageInterval{temp});
    end
    invertedTemps = 1./Temps;
    xmin = 1/(max(Temps)+1);
    xmax = 1/(min(Temps)-1);
    xrange = [xmin; xmax];
    
%convert to ln space
    lnPseudoKperStage = cell(length(pseudoKperStage),1);
    for temp = 1:length(pseudoKperStage)
        lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
    end
    
%iterate over stages to plot each stage
%initialize data arrays
    slope = zeros(1,length(lnPseudoKperStage{1}));
    intercept = zeros(1,length(lnPseudoKperStage{1}));
    yIntCI = zeros(2,length(lnPseudoKperStage{1}));
    stageActCi = zeros(2,length(lnPseudoKperStage{1}));

%initialize variable to hold error Statistic
    standErr = nan(length(invertedTemps),length(lnPseudoKperStage{1}));

%for all use -> length(lnPseudoKperStage{1})
    [numRow,numCol] = calcRectangle(length(lnPseudoKperStage{1}));
    pannelCount = 1;
    sampleSize = NaN(1,length(lnPseudoKperStage{1}));
    overAllSampleSize = NaN(1,length(lnPseudoKperStage{1}));
    numRuns = 5000;
    bsEas = NaN(numRuns,length(lnPseudoKperStage{1}));
    bsAs = NaN(numRuns,length(lnPseudoKperStage{1}));
for stage = [1,3:length(lnPseudoKperStage{1})]
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps
        fittabeEmNum = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0%1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                end
            end
            overAllEms = overAllEms + sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
        end
    
    %initialize and store subplot data
        sp(pannelCount) = subplot(numRow,numCol,pannelCount);   

    %initialize arrays to hold all fittable points    
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);  
        emNum = 1;

        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
    %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1); 
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
            standErr(temp,stage) = yStdErr;
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
            
        %plotting error bar and mean
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                viable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);
        size(allX)
        %record sample size
            sampleSize(stage) = length(allData);
            overAllSampleSize(stage) = length(overAllY);
    %store regression stats
        %individual points
        
            stagRegT.numPoints(stage) = length(allData);
            stagRegT.SX(stage) = sum(allData(:,1));
            stagRegT.SY(stage) = sum(allData(:,2));
            stagRegT.SXX(stage) = sum(allData(:,1).^2);
            stagRegT.SYY(stage) = sum(allData(:,2).^2);
            stagRegT.SXY(stage) = sum(allData(:,1).*allData(:,2));
            stagRegT.AveX(stage) = mean(allData(:,1));
            stagRegT.AveY(stage) = mean(allData(:,2));

        if stage == 1
            flyRegT = stagRegT(stage,:);
        end
        
    %bootstrapping
    %{
        for i = 1:numRuns
            rsamp = randsample(1:length(allData),length(allData),true);
            bsData = allData(rsamp,:);
            lineFit = fit(bsData(:,1), bsData(:,2), 'poly1');
            bsEas(i,stage) = (-lineFit.p1*R)*tokJConvert;
            bsAs(i,stage) = lineFit.p2;
        end
    %}
    %plot linear regression
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');
        
        xra = min(xrange):0.00001:max(xrange);

        curvOFit = polyfit(overAllX, overAllY,2);
        yOCurv = polyval(curvOFit,xra);
        ocm = plot(xra,yOCurv,'--r', 'LineWidth',2.5);
        
    %plot fit lines
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'k', 'LineWidth',2.5);
        text(1./(273.15+27.5),-5.25,strcat("E_a = ",num2str(round(...
            (-lineFit.p1*R)*tokJConvert,0))," kJ/mol"),'Color','k',...
            'FontSize',23); 

    %Plot parameters
        %int and slope
            slope(stage) = lineFit.p1;
            intercept(stage) = lineFit.p2;
        %get confidence intervals
            ci = confint(lineFit,0.68);
            stageActCi(:,stage) = ci(:,1);
            yIntCI(:,stage) = ci(:,2);
            
            
        CIs = confint(lineFit,0.68);    
        plusCI = lineFit.p1 - CIs(1,1);
        round((plusCI*R)*tokJConvert,0);
        
    %beautify
        %Set consistant limits to compare
            ylim([-13 -4.5]) 
            yticks([-12 -10 -8 -6 -4])
            xlim(xrange);

        %adjust xticks
            %xticks(1./[299.1500 285.1500]);
            %xticklabels({'1/26','1/12'});
            xticks(1./[303.1500 294.1500 285.1500]);
            xticklabels({'1/30','1/21','1/12'});

        %set labels on outer most subplots, remove from interior        
            if stage==1||stage==6||stage==10
                ylabel('ln(1/time [s^{-1}])');
            else
                set(gca,'yticklabel',[])
            end
            if stage>=9 && stage<=12
                xlabel('1/(T [^oC] + 273)');
            else
                set(gca,'xticklabel',[])
            end

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 23);

        %place titles as text inside plots
            if stage == 1
                text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(stage),...
                    " to ",(stageAbbsNoT0(stage+1))))},...
                    'FontSize', 24,'HorizontalAlignment','center');
            else
                text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(stage-1),...
                    " to ",(stageAbbsNoT0(stage))))}, 'FontSize', 24,...
                    'HorizontalAlignment','center');
            end
    pannelCount = pannelCount + 1;
end

%Adjust Subplot positions to tighten image
    %Shift left  
        for stage = 1:length(obsScores)-1
            if stage ~= 1 && stage ~= 5 && stage ~= 9
                prevPos = get(sp(stage-1),'Position');
                prevPos = [prevPos(1)+prevPos(3)+0.0075 prevPos(2) ...
                    prevPos(3) prevPos(4)];
                set(sp(stage),'Position', prevPos)
            end
        end

    %shift up
        for stage = 1:length(obsScores)-1
            
            if stage >= 5 && stage <=8
                prevPos = get(sp(1),'Position');
                curPos = get(sp(stage),'Position');
                curPos(2) = prevPos(2)-prevPos(4)-0.015;
                set(sp(stage),'Position', curPos)
            end
            if stage >= 9 && stage <=12
                prevPos = get(sp(5),'Position');
                curPos = get(sp(stage),'Position');
                curPos(2) = prevPos(2)-prevPos(4)-0.015;
                set(sp(stage),'Position', curPos)
            end
        end

%Arrhenius parameters calculation
    stageActEa = (-slope*R)*tokJConvert;
    stageActA = exp(intercept);
    stageActCi = -(stageActCi)*R*tokJConvert;
%%
%Reviewer Reponse
%"violin plots"
means = NaN(1,12);
medians = NaN(1,12);
%own version
figure('Position',[10 10 1500 900]);
for i = [1,3:12]
    edges = [40:0.25:110];
    if i == 1
        h = subplot(1,12,i);
    else
        h = subplot(1,12,i-1);
    end
    histogram(bsEas(:,i),edges,'EdgeAlpha',0);
       %set(h,'Position',[0.1308 0.1100 0.7890 0.1100])
    %calc and plot mean/median
        meanEas = mean(bsEas(:,i));
        means(i) = meanEas;
        medianEas = median(bsEas(:,i));
        standErr = std(bsEas(:,i));
    %plot median and mean
        black = xline(medianEas,'r--','linewidth',3);
            %correct these so the dont "spill over" the hist
        red = xline(medianEas-standErr,'k--','linewidth',3);
        red = xline(medianEas+standErr,'k--','linewidth',3);
    
    xlim([40 110])
    ylim([0 500])
    set(h,'view',[90 -90])
    axis off
    if i ==1
        legend([black red],["Median E_a","68% CI in E_a"],'Location','northeast',...
    'LineWidth',2.5,'FontSize',42)
    else
    end
end
i = i-1;
axes('position',[0.1308-0.01 0.1100-.00155 0.7890-0.1308+0.0571 0.8150],...
    'color','none') %[left bottom width height]
    %run once to determine proper axes size/positions, check h
    %subplot1pos = [0.1308 0.1100 0.0462 0.8150]
    %subplot2pos = [0.1959 0.1100 0.0471 0.8150]
    %subplot11pos = [0.7890 0.1100 0.0501 0.8150]
xticks([0.015:0.0922:1])
xlabs = strcat(stageAbbs(1:11)," - ", stageAbbs(2:12));
xticklabels(xlabs)
%xtickangle(45)
set(gca,'lineWidth',4.5,'FontSize',35)
xlabel("Fly intervals",'FontWeight','Bold','FontSize',38)
yticks([0:1/7:1])
yticklabels([40:10:110])
ylabel("Apparent activation energy (kJ/mol)",'FontWeight','Bold','FontSize',38)
%legend([black red],["Median E_a","68% CI in E_a"],'Location','northeast',...
 %   'LineWidth',2.5,'FontSize',30)

%title(["Distribution of Bootstrapped Fly Interval Activation Energies"],'FontSize',30)

((stageActEa-means)./stageActEa)*100
((stageActEa-medians)./stageActEa)*100
%%
%Figure 4A
%Arrhenius Stage-wise
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    
    embryoMass = 8.99*10^-3;% mg

%Plot Stages
    figure('Position',[10 10 1050 550]);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    obsScores = (1:numScores);
    perStageInterval = cellfun(@(x) x(:,obsScores),perStage,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 sizeData-1 sizeData];
    viableTemps = [2:sizeData-1];
    
%Mass correction by t/m^(1/4)
    %perStageInterval = cellfun(@(x) x(:,obsScores)./embryoMass^0.25,perStage,'UniformOutput',false);    
    
%Initialize tables
    %Statistical Table
        stagRegM = zeros(length(obsScores),8);
        stagRegT = array2table(stagRegM);
        stagRegT.Properties.VariableNames = {'numPoints','SX','SY',...
            'SXX','SYY','SXY','AveX','AveY'};
        stagRegT.Properties.RowNames = stageNamesNoT0;

%convert to pseudo Ks and inverse temps
    pseudoKperStage = cell(length(perStageInterval),1);
    for temp = 1:length(perStageInterval)
        pseudoKperStage{temp} = 1./abs(perStageInterval{temp});
    end
    invertedTemps = 1./Temps;
    xmin = 1/(max(Temps)+1);
    xmax = 1/(min(Temps)-1);
    xrange = [xmin; xmax];
    
%convert to ln space
    lnPseudoKperStage = cell(length(pseudoKperStage),1);
    for temp = 1:length(pseudoKperStage)
        lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
    end
    
%iterate over stages to plot each stage
%initialize data arrays
    slope = zeros(1,length(lnPseudoKperStage{1}));
    intercept = zeros(1,length(lnPseudoKperStage{1}));
    yIntCI = zeros(2,length(lnPseudoKperStage{1}));
    stageActCi = zeros(2,length(lnPseudoKperStage{1}));

%initialize variable to hold error Statistic
    standErr = nan(length(invertedTemps),length(lnPseudoKperStage{1}));

%for all use -> length(lnPseudoKperStage{1})
    [numRow,numCol] = calcRectangle(length(3:8));
    pannelCount = 1;
    sampleSize = NaN(1,length(lnPseudoKperStage{1}));
    overAllSampleSize = NaN(1,length(lnPseudoKperStage{1}));
    numRuns = 5000;
    bsEas = NaN(numRuns,length(lnPseudoKperStage{1}));
    bsAs = NaN(numRuns,length(lnPseudoKperStage{1}));
for stage = [3:8]
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps
        fittabeEmNum = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0%1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                end
            end
            overAllEms = overAllEms + sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
        end
    
    %initialize and store subplot data
        sp(pannelCount) = subplot(numRow,numCol,pannelCount);   

    %initialize arrays to hold all fittable points    
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);  
        emNum = 1;

        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
    %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1); 
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
            standErr(temp,stage) = yStdErr;
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
            
        %build all data for specific temp for linear regression, while
        %plotting error bar
            %plotting error bar and mean
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                viable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
            
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);
        %record sample size
            sampleSize(stage) = length(allData);
            overAllSampleSize(stage) = length(overAllY);
    %store regression stats
        %individual points
        
            stagRegT.numPoints(stage) = length(allData);
            stagRegT.SX(stage) = sum(allData(:,1));
            stagRegT.SY(stage) = sum(allData(:,2));
            stagRegT.SXX(stage) = sum(allData(:,1).^2);
            stagRegT.SYY(stage) = sum(allData(:,2).^2);
            stagRegT.SXY(stage) = sum(allData(:,1).*allData(:,2));
            stagRegT.AveX(stage) = mean(allData(:,1));
            stagRegT.AveY(stage) = mean(allData(:,2));

        if stage == 1
            flyRegT = stagRegT(stage,:);
        end
        
    %bootstrapping
    %{
        for i = 1:numRuns
            rsamp = randsample(1:length(allData),length(allData),true);
            bsData = allData(rsamp,:);
            lineFit = fit(bsData(:,1), bsData(:,2), 'poly1');
            bsEas(i,stage) = (-lineFit.p1*R)*tokJConvert;
            bsAs(i,stage) = lineFit.p2;
        end
    %}
    %plot linear regression
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');
        
        xra = min(xrange):0.00001:max(xrange);

        curvOFit = polyfit(overAllX, overAllY,2);
        yOCurv = polyval(curvOFit,xra);
        ocm = plot(xra,yOCurv,'--r', 'LineWidth',2.5);
        
    %plot fit lines
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'k', 'LineWidth',2.5);
        text(1./(273.15+26.5),-5.25,strcat("E_a = ",num2str(round(...
            (-lineFit.p1*R)*tokJConvert,0))," kJ/mol"),'Color','k',...
            'FontSize',23); 

    %Plot parameters
        %int and slope
            slope(stage) = lineFit.p1;
            intercept(stage) = lineFit.p2;
        %get confidence intervals
            ci = confint(lineFit,0.68);
            stageActCi(:,stage) = ci(:,1);
            yIntCI(:,stage) = ci(:,2);
            
            
        CIs = confint(lineFit,0.68);    
        plusCI = lineFit.p1 - CIs(1,1);
        round((plusCI*R)*tokJConvert,0);
        
    %beautify
        %Set consistant limits to compare
            ylim([-13 -4.5]) 
            yticks([-12 -10 -8 -6 -4])
            xlim(xrange);

        %adjust xticks
            %xticks(1./[299.1500 285.1500]);
            %xticklabels({'1/26','1/12'});
            xticks(1./[303.1500 294.1500 285.1500]);
            xticklabels({'1/30','1/21','1/12'});

        %set labels on outer most subplots, remove from interior        
            if stage==3||stage==6
                ylabel('ln(1/time [s^{-1}])');
            else
                set(gca,'yticklabel',[])
            end
            if stage>=6 && stage<=8
                xlabel('1/(T [^oC] + 273)');
            else
                set(gca,'xticklabel',[])
            end

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 23);

        %place titles as text inside plots
            if stage == 1
                text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(stage),...
                    " to ",(stageAbbsNoT0(stage+1))))},...
                    'FontSize', 24,'HorizontalAlignment','center');
            else
                text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(stage-1),...
                    " to ",(stageAbbsNoT0(stage))))}, 'FontSize', 24,...
                    'HorizontalAlignment','center');
            end
    pannelCount = pannelCount + 1;
end
%
%Adjust Subplot positions to tighten image
    %Shift left  
        for stage = 1:length(3:8)
            if stage ~= 1 && stage ~= 4
                prevPos = get(sp(stage-1),'Position');
                prevPos = [prevPos(1)+prevPos(3)+0.0075 prevPos(2) ...
                    prevPos(3) prevPos(4)];
                set(sp(stage),'Position', prevPos)
            end
        end

    %shift up
        for stage = 1:length(3:8)      
            if stage >= 4 && stage <=6
                prevPos = get(sp(1),'Position');
                curPos = get(sp(stage),'Position');
                curPos(2) = prevPos(2)-prevPos(4)-0.015;
                set(sp(stage),'Position', curPos)
            end
        end
%{
%Arrhenius parameters calculation
    stageActEa = (-slope*R)*tokJConvert;
    stageActA = exp(intercept);
    stageActCi = -(stageActCi)*R*tokJConvert;
%}
%%    
%Figure 2A
%Arrhenius Stage-wise, just 2
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;

%Plot Stages
    figure('Position',[10 10 1000 400]);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    obsScores = (1:numScores);
    perStageInterval = cellfun(@(x) x(:,obsScores),perStage,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 sizeData-1 sizeData];
    viableTemps = [2:sizeData-1];

%convert to pseudo Ks and inverse temps
    pseudoKperStage = cell(length(perStageInterval),1);
    for temp = 1:length(perStageInterval)
        pseudoKperStage{temp} = 1./abs(perStageInterval{temp});
    end
    invertedTemps = 1./Temps;
    xmin = 1/(max(Temps)+1);
    xmax = 1/(min(Temps)-1);
    xrange = [xmin; xmax];
    
%convert to ln space
    lnPseudoKperStage = cell(length(pseudoKperStage),1);
    for temp = 1:length(pseudoKperStage)
        lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
    end
    
%iterate over stages to plot each stage
%initialize data arrays
    slope = zeros(1,length(lnPseudoKperStage{1}));
%initialize variable to hold error Statistic
    standErr = nan(length(invertedTemps),length(lnPseudoKperStage{1}));

%for all use -> length(lnPseudoKperStage{1})
    pannelCount = 1;
for stage = [6,7]
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps
        fittabeEmNum = 0;
        fittableMean = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0%1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                    fittableMean = fittableMean + 1;                    
                end
            end
        end
    
    %initialize and store subplot data
        subplot(1,2,pannelCount);   

    %initialize arrays to hold all fittable points    
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);  
        emNum = 1;
        allMeans = nan(fittableMean,1);
        allXM = nan(fittableMean,1);
        mNum = 1;
        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
    %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1);
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
            standErr(temp,stage) = yStdErr;
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
            
        %build all data for specific temp for linear regression, while
        %plotting error bar
        %plotting error bar and mean
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                viable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);

    %plot linear regression
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');

    %plot fit lines
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'k', 'LineWidth',2.5);
        CIs = confint(lineFit,0.95);
        plusCI = lineFit.p1 - CIs(1,1);
        minusCI = lineFit.p1 - CIs(2,1);
        text(1./(273.15+20.5),-5.55,strcat("E_a = ",num2str(round(...
            (-lineFit.p1*R)*tokJConvert,2,'significant'))," +/- ", num2str(round(...
            (plusCI*R)*tokJConvert,2,'significant'))),'Color','k','FontSize',24)
        text(1./(273.15+15.5),-6.0,"kJ/mol",'Color','k','FontSize',24)

    %Plot parameters
        %int and slope
            slope(stage) = lineFit.p1;
            intercept(stage) = lineFit.p2;
        %get confidence intervals
            ci = confint(lineFit,0.68);
            stageActCi(:,stage) = ci(:,1);
            yIntCI(:,stage) = ci(:,2);
        
    %beautify
        %Set consistant limits to compare
            ylim([-10 -5]) 
                %for all plots, to compare
            yticks([-11 -9 -7 -5])
            xlim(xrange);

        %adjust xticks
            %xticks(1./[299.1500 285.1500]);
            %xticklabels({'1/26','1/12'});
            xticks(1./[303.1500 294.1500 285.1500]);
            xticklabels({'1/30','1/21','1/12'});

        %set labels on outer most subplots, remove from interior        
            ylabel('ln(1/time [s^{-1}])');
            xlabel('1/(T [^oC] + 273)');

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 28);

        %place titles as text inside plots
            title(strcat(stageAbbsNoT0(stage-1)," to ",...
                stageAbbsNoT0(stage)),'FontSize', 28);
         
    pannelCount = pannelCount + 1;
end

%%    
%Abstract
%Arrhenius Stage-wise, just 2
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;

%Plot Stages
    figure('Position',[10 10 1200 400]);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    obsScores = (1:numScores);
    perStageInterval = cellfun(@(x) x(:,obsScores),perStage,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 sizeData-1 sizeData];
    viableTemps = [2:sizeData-1];
    
 %Initialize fit comparison Table
	AICTable = table('Size',[1,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
	AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};

%convert to pseudo Ks and inverse temps
    pseudoKperStage = cell(length(perStageInterval),1);
    for temp = 1:length(perStageInterval)
        pseudoKperStage{temp} = 1./abs(perStageInterval{temp});
    end
    invertedTemps = 1./Temps;
    xmin = 1/(max(Temps)+1);
    xmax = 1/(min(Temps)-1);
    xrange = [xmin; xmax];
    
%convert to ln space
    lnPseudoKperStage = cell(length(pseudoKperStage),1);
    for temp = 1:length(pseudoKperStage)
        lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
    end
    
%iterate over stages to plot each stage
%initialize data arrays
    slope = zeros(1,length(lnPseudoKperStage{1}));
%initialize variable to hold error Statistic
    standErr = nan(length(invertedTemps),length(lnPseudoKperStage{1}));

%for all use -> length(lnPseudoKperStage{1})
    pannelCount = 1;
for stage = [6,7]
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps
        fittabeEmNum = 0;
        fittableMean = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0%1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                    fittableMean = fittableMean + 1;                    
                end
            end
        end
    
    %initialize and store subplot data
        subplot(1,2,pannelCount);   

    %initialize arrays to hold all fittable points    
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);  
        emNum = 1;
        allMeans = nan(fittableMean,1);
        allXM = nan(fittableMean,1);
        mNum = 1;
        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
    %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1);
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
            standErr(temp,stage) = yStdErr;
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
            
        %build all data for specific temp for linear regression, while
        %plotting error bar
        %plotting error bar and mean
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                viable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(1,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);

    %plot linear regression
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');

    %plot fit lines
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'Color',[34,139,34]/255, 'LineWidth',2.5);
        CIs = confint(lineFit,0.95);
        plusCI = lineFit.p1 - CIs(1,1);
        minusCI = lineFit.p1 - CIs(2,1);
        text(1./(273.15+20.5),-5.55,"E_a = ",'Color','k','FontSize',28)
        text(1./(273.15+16.525),-5.43,strcat(num2str(round(...
            (-lineFit.p1*R)*tokJConvert,2,'significant'))," +/- ", num2str(round(...
            (plusCI*R)*tokJConvert,2,'significant'))),'Color',[34,139,34]/255,'FontSize',28)
        text(1./(273.15+15.5),-6.0,"kJ/mol",'Color','k','FontSize',28)
        
        oLinFit = polyfit(overAllX, overAllY, 1);
        oLin = polyval(oLinFit,xrange);
        %lm = plot(xrange,oLin, 'r--', 'LineWidth',2.5);
        oCurvFit = polyfit(overAllX, overAllY, 2);
        resTemps = min(xrange):0.00001:max(xrange);
        oCurv = polyval(oCurvFit,resTemps);
        curv = plot(resTemps,oCurv, 'r', 'LineWidth',2.5);          
        
        AICTable.numPoints = length(overAllY);
        AICTable.yi = num2cell(overAllY,1);
        AICTable.xi = num2cell(overAllX,1);
        
        BICcomp = fitBICCompCalc(AICTable,oLinFit,oCurvFit,1)
    %Plot parameters
        %int and slope
            slope(stage) = lineFit.p1;
            intercept(stage) = lineFit.p2;
        %get confidence intervals
            ci = confint(lineFit,0.68);
            stageActCi(:,stage) = ci(:,1);
            yIntCI(:,stage) = ci(:,2);
        
    %beautify
        %Set consistant limits to compare
            ylim([-10 -5]) 
                %for all plots, to compare
            yticks([])
            xlim(xrange);

        %adjust xticks
            %xticks(1./[299.1500 285.1500]);
            %xticklabels({'1/26','1/12'});
            xticks([]);

        %set labels on outer most subplots, remove from interior        
            ylabel('ln(1/time)');
            xlabel('1/T (K^{-1})');

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 32);

        %place titles as text inside plots
            title({strcat(stageNames(stage-1)),strcat(" to ",...
                stageNames(stage))},'FontSize', 28);
        
        %report BIC
            %text(1./(273.15+33),-9.25,strcat(...
                %'$ln(\frac{Likelihood_{Quadratic}}{Likelihood_{Linear}}) = ',num2str(round(log(BICcomp),...
                %2,'significant')),'$'),'Color','k','FontSize',37,'Interpreter','latex')
            text(1./(273.15+34),-9.25,"ln(",'Color','k','FontSize',37,'Interpreter','latex')
            text(1./(273.1+15.5),-9.25,strcat(...
                "$)= ",num2str(round(log(BICcomp),2,'significant')),'$'),...
                'Color','k','FontSize',37,'Interpreter','latex')
            text(1./(273.15+30.85),-8.95,"Likelihood Quadratic",'Color','r','FontSize',25,'Interpreter','latex')
            text(1./(273.15+29.5),-9.55,"Likelihood Linear",'Color',[34,139,34]/255,'FontSize',25,'Interpreter','latex')
            plot([1./(273.15+31) 1./(273.1+15.20)],[-9.25 -9.25],'k','LineWidth',2.5)
         
    pannelCount = pannelCount + 1;
end


%%
%Figure 3A
%Arrhenius Accumulated
   figure('Position',[100 100 1250 450]);
    
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    
    embryoMass = 8.99*10^-3; %mg
    
%remove T0 time point from analysis
    obsScores = (1:numScores);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [91];%[1 sizeData-1 sizeData];%
        %[91] for all data 
    viableTemps = [2:sizeData-1];
    
%Mass correction by t/m^(1/4)
    %accumulatedInterval = cellfun(@(x) x(:,obsScores)./embryoMass^0.25,accumulated,'UniformOutput',false); 
    
%Initialize tables
    
%initialize tables
    %stats table
        accRegM = zeros(length(obsScores),8);
        accRegT = array2table(accRegM);
        accRegT.Properties.VariableNames = {'numPoints','SX','SY','SXX','SYY','SXY','AveX','AveY'};
        accRegT.Properties.RowNames = stageNamesNoT0;
        
%Initialize fit comparison Table
	AICTable = table('Size',[1,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
	AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};
        
%convert to pseudo Ks and inverse time
    pseudoKperStage = cell(length(accumulatedInterval),1);
    for temp = 1:length(accumulatedInterval)
        pseudoKperStage{temp} = 1./accumulatedInterval{temp};
    end
    invertedTemps = 1./Temps;
    xmin = 1/(max(Temps)+1);
    xmax = 1/(min(Temps)-1);
    xrange = [xmin; xmax];
    
%convert to ln space
    lnPseudoKperStage = cell(length(pseudoKperStage),1);
    for temp = 1:length(pseudoKperStage)
        lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
    end
    
%iterate over stages to plot each stage
%initialize data arrays
    slope = zeros(1,length(lnPseudoKperStage{1}));
    intercept = zeros(1,length(lnPseudoKperStage{1}));
    accActCi = zeros(2,length(lnPseudoKperStage{1}));
    
    justThree = [8,12];
    justThreeCounter = 1;

%for all use -> length(lnPseudoKperStage{1}) 
for stage = justThree
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps only
        fittabeEmNum = 0;
        fittableMean = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) >0
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                        fittableMean = fittableMean + 1;                    
                end
            end
        end
        
    %initialize arrays to hold all data X and Y, and for means
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);
        emNum = 1;
        
    %initialize subplots
        subplot(1,2,justThreeCounter)
        
   %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1); 
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
                %standard Error
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
            
        %build all data for specific temp for linear regression, while
        %plotting error bar
            meanLnPseudoK = mean(data);
            %plotting error bar and mean
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                viable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end

            %{
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %fit to individuals
                fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);            
            elseif length(data)==1
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(3,:),'LineWidth',2.5);                
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);
            end
            %}
    end
    
    %combine allX/Y
        allData = horzcat(allX,allY);
        
    %store regression stats
        
        accRegT.numPoints(stage) = length(allData);
        accRegT.SX(stage) = sum(allData(:,1));
        accRegT.SY(stage) = sum(allData(:,2));
        accRegT.SXX(stage) = sum(allData(:,1).^2);
        accRegT.SYY(stage) = sum(allData(:,2).^2);
        accRegT.SXY(stage) = sum(allData(:,1).*allData(:,2));
        accRegT.AveX(stage) = mean(allData(:,1));
        accRegT.AveY(stage) = mean(allData(:,2));             
        
        AICTable.numPoints = length(allY);
        AICTable.yi = num2cell(allY,1);
        AICTable.xi = num2cell(allX,1);
    
    %plot linear regression
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');
        mdl = polyfit(allData(:,1), allData(:,2), 1);
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'k', 'LineWidth',2.5);
        
        curvFit = polyfit(allData(:,1), allData(:,2),2);
        xra = min(xrange):0.00001:max(xrange);
        yCurv = polyval(curvFit,xra);
        cm = plot(xra,yCurv,'--r', 'LineWidth',2.5);
    
    %calculate slope
        slope(stage) = lineFit.p1;
        intercept(stage) = lineFit.p2;
        
    %get confidence intervals
        ci = confint(lineFit,0.68);
        accActCi(:,stage) = ci(:,1);
            
        hold on;
        
        BICcomp = fitBICCompCalc(AICTable,mdl,curvFit,1);
    %beautify
        %ylabels
            ylabel('ln(1/time [s^{-1}])');
            meanY = round(real(mean(allY,'omitnan')),0);
            theoryMax = polyval([lineFit.p1 lineFit.p2],xrange(1));
            theoryMin = polyval([lineFit.p1 lineFit.p2],xrange(2));
            ylim([theoryMin-0.5,theoryMax]);
            ylabs = yticklabels;
            ticks = {ylabs{cellfun(@(x) floor(str2num(x)) == str2num(x),...
                ylabs)}};
            yticks(cell2mat({cellfun(@str2num,ticks)}));
            yticklabels(ticks);
            text(1./(273.15+18.5),theoryMax-0.35,strcat(...
                '$ln(\frac{L_{Q}}{L_{L}}) = ',num2str(round(log(BICcomp),...
                2,'significant')),'$'),'Color','k','FontSize',30,'Interpreter','latex')

        %xlabels
            xlabel('1/(T [^oC] + 273)');
            xlim(xrange);
            %xticks(1./[303.15 293.15 283.15]);
            %xticklabels({'1/30','1/20','1/10'});
            xticks(1./[303.1500 294.1500 285.1500]);
            xticklabels({'1/30','1/21','1/12'});
            
        %gca, titles, legend
            set(gca,'linewidth', 3, 'FontSize', 28);
            title(strcat("A to ",stageAbbsNoT0(stage)),'FontSize',28);
            legend([lm cm], {"Linear Fit","Quadratic"},'Location',...
                'southwest','FontSize',25);
    
    justThreeCounter = justThreeCounter + 1;
end

%Ea calculation
    accActEa = -(slope)*R*tokJConvert;
    accActA = exp(intercept);
    accActCi = -(accActCi)*R*tokJConvert;   
    
%%
%Statistical Significance calculations for stage-wise
    %run the stages section first
figure('Position',[100 100 1100 900]);
%h = heatmap(tbl,xvar,yvar)
%create a table
    %each plot against eachother 1:10x1:10
    %fill table with pval = ancova(numPoints,SX,SY,SXX,SYY,SXY,AveX,AveY)
        %numPoints,SX,SY,SXX,SYY,SXY,AveX,AveY: comes from two regressions
        %to compare
        %feed regressionAnalysis table of numPoints,SX,SY,SXX,SYY,SXY,AveX,
        %AveY made of just concatination of 2 regressions in question

%adjjust stageNames to include only obsScores, from above code chunk
    stageNamesNoT0 = stageNames(obsScores);
    stageAbbsNoT0 = stageAbbs(obsScores);  

%Initialize presentation data matrix
    pValM = nan(length(obsScores)-1,length(obsScores)-1);
    pwrVals = nan(length(obsScores)-1,length(obsScores)-1);
    
for x = 2:length(obsScores)
    for y = x:length(obsScores)
        if x == y
            pValM(y-1,x-1) = 1;
                %stages compared against same stage should = 1
            pwrVals(y-1,x-1) = 0;
        elseif x == 2
            compStages = {char(stageNamesNoT0(x-1)),char(stageNamesNoT0(y))};
            compRegs = stagRegT(compStages,:);
        %run stage regression table through home-made ancova function
            pValM(y-1,x-1) = ancova(compRegs);
            pwrVals(y-1,x-1) = calcPower(compRegs);
        elseif y == 2
            compStages = {char(stageNamesNoT0(x)),char(stageNamesNoT0(y-1))};
            compRegs = stagRegT(compStages,:);
        %run stage regression table through home-made ancova function
            pValM(y-1,x-1) = ancova(compRegs);
            pwrVals(y-1,x-1) = calcPower(compRegs);
        else
            compStages = {char(stageNamesNoT0(x)),char(stageNamesNoT0(y))};
            compRegs = stagRegT(compStages,:);
        %run stage regression table through home-made ancova function
            pValM(y-1,x-1) = ancova(compRegs);
            pwrVals(y-1,x-1) = calcPower(compRegs);
        end
    end  
end

%plot 4 plots, one of p-values, then color coded for thhreshold sig at
%0.05,0.01,0.001
%pvals
    %subplot(2,2,1)
    h1 = heatmap(strcat(stageAbbsNoT0(1:11),'-',stageAbbsNoT0(2:12)),...
        strcat(stageAbbsNoT0(1:11),'-',stageAbbsNoT0(2:12)),...
        round(pValM,4),'CellLabelColor','none','CellLabelFormat','%.1e');

    %h1.Title = {'\fontsize{40}Fly ANCOVA p-values',...
        %'\fontsize{40}comparing individual intervals'};
    h1.Title = {'\fontsize{40}Fly ANCOVA p-values'};
    h1.FontSize = 15;
    h1.MissingDataColor = [1 1 1];
    h1.GridVisible = 'off';
    %colormap
        caxis([0 1])
        noSig = horzcat(repmat(0.6,950,1),repmat(0.6,950,1),ones(950,1));
        singSig = horzcat(repmat(0.5,40,1),repmat(0,40,1),repmat(0.5,40,1));
        doubSig = horzcat(ones(9,1),repmat(0.3,9,1),ones(9,1));
        tripSig = horzcat(ones(2,1),repmat(0.3,2,1),repmat(0.1,2,1));
        map = vertcat(tripSig,doubSig,singSig,noSig);
        colormap(h1,map)
        
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Interval 1}')
    ylabel({'\fontsize{38}\bf{Interval 2}'})
    
    
    %adjust color bar by covering up
        hmp = h1.Position; 
        cbax = axes('Position',[sum(hmp([1,3])), 0, 1-sum(hmp([1,3]))+0.2, 1],...
            'XTick',[], 'YTick', [], 'Color',[1 1 1]);
        cbax.XAxis.Visible = 'off';
        cbax.YAxis.Visible = 'off';
    % Set the new axis color map to match the 
    % heatmap's colormap, cutoff at 0.1 to create trunkated version
        noSig = horzcat(repmat(0.6,50,1),repmat(0.6,50,1),ones(50,1));
        newMap = vertcat(tripSig,doubSig,singSig,noSig);
        cbax.Colormap = newMap; 
    % Add a colorbar the same vertical position as the heatmap
        cbh = colorbar(cbax,'Position',[.90, hmp(2), .03, hmp(4)],'AxisLocation','in'); 
    % Set the limits to 0:1 and set the ticks 
        cbh.Limits = [0,1]; 
        nColors = size(h1.Colormap,1); 
        cbh.Ticks = [2/100,11/100,50.5/100];
    % Set the new labels for each tick
        cbh.TickLabels = {'1E-3','1E-2','5E-2'}; 
    % Set the colorbar fontsize to the same value as heatmap fontsize
        cbh.FontSize = h1.FontSize;
        cbh.LineWidth = 4;
        
        prevPos = get(h1,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(h1,'Position', prevPos)
        
        prevPos = get(cbh,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(cbh,'Position', prevPos)
        cbheight = prevPos(4);
        
        prevPos = get(cbax,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(cbax,'Position', prevPos)

    IPx = h1.InnerPosition(1);
    IPy = h1.InnerPosition(2);
    textAxis = axes('Position',[0 0 1 1],'visible','off');
    scaleX = h1.InnerPosition(3)/11;
    scaleY = h1.InnerPosition(4)/11;
    
%
for x = 1:length(obsScores)-1
    for y = x:length(obsScores)-1
        txt = sprintf('%0.0E',pValM(y,x));
        if str2num(txt) < 1 
            txt = replace(txt,'-0','-');
        elseif str2num(txt) > 1 || str2num(txt) == 1 
            txt = replace(txt,'+0','+');
        end
        xShift = IPx - 0.5*scaleX;
        yShift = IPy + h1.InnerPosition(4) + 0.5*scaleY;
        text(xShift + scaleX*x,yShift - scaleY*y,txt,'Color','k','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
end
hold on
plot([IPx h1.InnerPosition(3)+IPx],[IPy IPy],'k','LineWidth',4)
plot([IPx IPx],[IPy h1.InnerPosition(4)+IPy],'k','LineWidth',4)
plot([IPx-0.01 h1.InnerPosition(3)+IPx],[h1.InnerPosition(4)+IPy h1.InnerPosition(4)+IPy],'w','LineWidth',4)
ylim([0 1])
xlim([0 1])

%} 
%{
%pwrValues
figure('Position',[100 100 1100 900]);        

    h1 = heatmap(strcat(stageAbbsNoT0(1:11),'-',stageAbbsNoT0(2:12)),strcat(stageAbbsNoT0(1:11),...
        '-',stageAbbsNoT0(2:12)),round(pwrVals,4),'CellLabelColor','black');
    h1.Title = {'\fontsize{40}Fly ANCOVA power analysis', '\fontsize{40}comparing individual intervals'};
    h1.FontSize = 15;
    h1.MissingDataColor = [1 1 1];
    h1.GridVisible = 'off';
    caxis([0 1])
    noSig = horzcat(repmat(0.6,799,1),repmat(0.6,799,1),ones(799,1));
    singSig = horzcat(repmat(0.5,200,1),repmat(0,200,1),repmat(0.5,200,1));
    map = vertcat(noSig,singSig);
    colormap(h1,map)
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Interval 1}')
    ylabel({'\fontsize{38}\bf{Interval 2}'})
       
    %adjust color bar by covering up
        hmp = h1.Position; 
        cbax = axes('Position',[sum(hmp([1,3])), 0, 1-sum(hmp([1,3]))+0.2, 1],...
            'XTick',[], 'YTick', [], 'Color',[1 1 1]);
        cbax.XAxis.Visible = 'off';
        cbax.YAxis.Visible = 'off';
    % Set the new axis color map to match the 
    % heatmap's colormap, cutoff at 0.1 to create trunkated version
        noSig = horzcat(repmat(0.6,800,1),repmat(0.6,800,1),ones(800,1));
        newMap = vertcat(noSig,singSig);
        cbax.Colormap = newMap; 
    % Add a colorbar the same vertical position as the heatmap
        cbh = colorbar(cbax,'Position',[.90, hmp(2), .03, hmp(4)],'AxisLocation','in'); 
    % Set the limits to 0:1 and set the ticks 
        cbh.Limits = [0,1]; 
        nColors = size(h1.Colormap,1); 
        cbh.Ticks = [0/100,80/100,100/100]; 
   % Set the new labels for each tick
        cbh.TickLabels = [0,0.8,1]; 
   % Set the colorbar fontsize to the same value as heatmap fontsize
        cbh.FontSize = h1.FontSize;
        
        prevPos = get(h1,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(h1,'Position', prevPos)
        
        prevPos = get(cbh,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(cbh,'Position', prevPos)
        
        prevPos = get(cbax,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(cbax,'Position', prevPos)
%}
%%
%Reviewer Reponse
%Power vs. N
fh = figure('Position',[100 100 800 650]);

%adjust stageNames to include only obsScores, from above code chunk
    stageNamesNoT0 = stageNames(obsScores);
    stageAbbsNoT0 = stageAbbs(obsScores);
%Initialize presentation data matrix
    numRuns = 50;
    pwrVals = nan(1,numRuns);
    pValM = nan(1,numRuns);

%choose stages to compare
    %stages are labeled by end stage
    stage1 = 6;
    stage2 = 7;
    compStages = {char(stageNamesNoT0(stage1)),char(stageNamesNoT0(stage2))};
    compRegs = stagRegT(compStages,:);
%run stage regression table through home-made ancova function
    for i = 1:numRuns
        pValM(i) = ancovaSpec(compRegs,i);
        pwrVals(i) = calcSpecPower(compRegs,i); 
    end
    
    %calculate actual statistics
        pValMStat = ancova(compRegs);
        pwrValStat = calcPower(compRegs); 
    
    hold on
    yyaxis right
    plot(1:numRuns,pwrVals,'linewidth',2)
    yline(pwrValStat,"LineWidth",2)
    ylabel("Power")
    
    yyaxis left
    plot(1:numRuns,pValM,'linewidth',2)
    yline(pValMStat,"LineWidth",2)
    ylabel("p-Values")
    
    xlabel("Sample Size per Condition Compared")
    title({"Fly Power and p-Value vs. Sample Size",strcat(stageAbbsNoT0(...
        stage1-1)," - ",stageAbbsNoT0(stage1)," against ",stageAbbsNoT0(...
        stage2-1)," - ",stageAbbsNoT0(stage2))})
    set(gca,'linewidth', 3, 'FontSize', 28);
 
%%
%Figure 2C
%Comparing Eas with statistical significance between stage-wise
%Must run the above stages code chunk first
figure('Position',[100 100 900 700]);

%parse Ea data and CIs to present
    spos = stageActCi(1,:) - stageActEa;
    sneg = stageActEa - stageActCi(2,:);
    nonScaled = [1,3:length(obsScores)];
 
%plot errorbar of data
    stageEaPlot = errorbar(1:length(obsScores)-1, stageActEa(nonScaled),...
        sneg(nonScaled), spos(nonScaled),'o','LineWidth',2.5);

%Non Scaled Aestetics   
    %even ticks
        xlim([0 12])
        xticks([1:12])
    %tick aestethics
        stageEaPlot.Color = 'b';
        xticklabels(strcat(stageAbbsNoT0(1:length(obsScores)-1),"-",...
            stageAbbsNoT0(2:length(obsScores))))
        xtickangle(45)
        
    set(gca,'linewidth', 3, 'FontSize', 28)
    %title([{"Apparent activation energy", " Fly developmental intervals"}],...
        %"FontSize", 28)

    ylabel("Apparent activation energy (kJ/mol)", "FontSize", 28,'FontWeight','Bold')
    xlabel("Fly interval", "FontSize", 28,'FontWeight','Bold')

%hard code brackets/braces based on the stats above
    hold on;
%use Eas+errorbar to determine bracket start point
	compare = [stageActEa(3)-spos(3) stageActEa(5)-spos(5)];
	largerEa = max(compare);
	smallerEa = min(compare);

%Draw Braces between significantly different Eas, Use ancova calc above
%DRAWBRACE([X1,Y1], [X2,Y2], W)
%new data

%upside down
	%between 8 an 9
    	y1 = stageActEa(6)-spos(6)-1;
    	drawbrace([6,y1], [5,y1], 10,'LineWidth',2.5,'Color','k')
         	text(5.27,y1-4.25,'\fontsize{25}***');
            
    ylim([40 120])

%%
%Response to Reviewers
%Comparing Eas with statistical significance between stage-wise
%Must run the above stages code chunk first
figure('Position',[0 0 1700 500]);

%parse Ea data and CIs to present
    spos = stageActCi(1,:) - stageActEa;
    sneg = stageActEa - stageActCi(2,:);
    nonScaled = [1 3:length(obsScores)];
    %scale by 22.2 C data
    scaleTempIndx = 7;
    tempScale = mean(accumulated{scaleTempIndx},1,'omitnan');
    normalizedScale = tempScale/tempScale(length(tempScale));
 
%plot errorbar of data
    %scaled
    meanXs = nan(length(nonScaled),1);
    for i = 1:length(nonScaled)
        meanX = mean([normalizedScale(i), normalizedScale(i+1)]);
        meanXs(i) = meanX;
        errs = diff([meanX,normalizedScale(i)]);
        stageEaPlot = errorbar(meanX, stageActEa(nonScaled(i)),...
            sneg(nonScaled(i)), spos(nonScaled(i)),'bo','LineWidth',2.5);
        hold on
    end 
    for i = 1:length(normalizedScale)
        xline(normalizedScale(i),'k--','LineWidth',1.5)
    end
%Scaled Aestetics    
    %developmentally scaled ticks
        xlim([0 1])
        xticks(meanXs)
    %tick aestethics
        stageEaPlot.Color = 'b';
        labs = strcat(stageAbbs(1:12-1),"-",stageAbbs(2:12));
        innerLabs = labs([4 6]);
        labs([4 6]) = "";
        xticklabels(labs)
        xtickangle(45)
        text(normalizedScale(4), 38, innerLabs(1),'Rotation',45,'HorizontalAlignment','left','FontSize',20)
        text(normalizedScale(6)+800/tempScale(length(tempScale)), 38, innerLabs(2),'Rotation',45,'HorizontalAlignment','left','FontSize',20)

        %labels
        set(gca,'linewidth', 3, 'FontSize', 20)
        %title([{"Apparent activation energies: Scaled by proportion of fly development"}],...
            %"FontSize", 25)

        ylabel("Apparent activation energy (kJ/mol)", "FontSize", 23,'FontWeight','Bold')
        title(strcat("Fly interval activation energies (Scaled by developmental interval times at ",...
            num2str(Temps(scaleTempIndx)-273.15), " ^oC)"), "FontSize", 23,'FontWeight','Bold')   
        ylim([40 120]);

%hard code brackets/braces based on the stats above
    hold on;
%use Eas+errorbar to determine bracket start point
	compare = [stageActEa(3)-spos(3) stageActEa(5)-spos(5)];
	largerEa = max(compare);
	smallerEa = min(compare);

%Draw Braces between significantly different Eas, Use ancova calc above
%DRAWBRACE([X1,Y1], [X2,Y2], W)
%new data
%upside down
	%between 8 an 9
    	y1 = stageActEa(6)-spos(6)-1;
    	drawbrace([meanXs(6),y1], [meanXs(5),y1], 10,'LineWidth',2.5,'Color','k')
         	text(mean([meanXs(6) meanXs(5)])-1000/tempScale(length(tempScale)),y1-7.3,'\fontsize{25}***');

    ylim([35 110])

%%
%Figure 3B & S4C
%Copied from frog code
%BIC Calculations for every stage interval
%Heatmap for AIC calculations
    
    obsScores = (1:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);    
    
%Initialize fit comparison Table
	AICTable = table('Size',[1,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
	AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};

%initialize matrix to hold errors for monte carlo
    monteTimeErrorsLinear = nan(numScores);
    intervalAllXLinear = cell(numScores);
    monteTimeErrorsExtremes = nan(numScores);
    intervalAllXExtremes = cell(numScores);

%initialize fit comparison arrays
	AICValues = zeros(numScores,4);
        
%initialize Data for Arrhenius space x and y
    lnPseudoKperStage = cell(length(accumulatedInterval),1);
    invertedTemps = 1./Temps;

%initialize presentation data matrecies    
    BICMat = nan(numScores);
    EaMat = nan(numScores);
    sampleSizeMatLin = nan(numScores);
    sampleSizeMatExt = nan(numScores);
    isConcaveBioFly = nan(numScores);
    isConvexBioFly = nan(numScores);
    isLinearBioFly = nan(numScores);
    
%control what temperatures to exclude from the linear fit BIC comparison
%calculation
    counter = 1;

%calculate poly n 1 and 2 fits for every developmental interval then compare 
%them using BIC
for p = 1:2
    if p == 1
        %here excluded temps = all temps - viable temps, so only both most
        %extreme
        excludedTemps = [1 length(lnPseudoKperStage)];
    else
        excludedTemps = [76];
    end
for startStage = 1:length(obsScores) 
    for endStage = startStage:length(obsScores) 
    %calculate interval times by diff between acc(stages)
        dataRange = cellfun(@(x) x(:,[startStage endStage]),accumulatedInterval,'UniformOutput',false);
        diffDataRange = cellfun(@(x) diff(x,1,2),dataRange,'UniformOutput',false);
    %change any zeros to Nans
        convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'[^0123456789.]0','NaN')).',diffDataRange,'UniformOutput',false);
        convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'^0[^0123456789.]','NaN')).',convertZeros,'UniformOutput',false);        
    %convert times to Arrhenius ln(k)
        lnPseudoKperStage = cellfun(@(x) log(1./x),convertZeros,'UniformOutput',false);
    %determine number of fitting embryos or means
        fittabeEmNum = 0;
        for temp = 1:length(lnPseudoKperStage)
        %if not one of the excluded temperature and more then one embryo
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp})) > 1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}));
                end
            end
        end
        
    %initialize arrays to hold all fittable points X,Y, and means        
    	allX = nan(fittabeEmNum,1);
    	allY = nan(fittabeEmNum,1);
        emNum = 1;
        
    %initialize tempararyError array
        tempararyError = nan(1,length(invertedTemps));
        
	%For every temperature build up the fitting data
        for temp = 1:length(invertedTemps)
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp};
            data = data(~isnan(data));
        %copy over temp errors
            N = size(data,1); 
                %number non-nan experiments in data
            if N > 1
                yStdErr = std(data)/sqrt(N-1); 
                tempararyError(temp) = yStdErr;
            end            
        %build fitting data
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0
                %if not one of the excluded temperature and more then one embryo
            %individual data points
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;

            elseif length(data)==1              
            else
            end
        end    
        %copy over errors for monte carlo
        if p ==1
            %linear
            monteTimeErrorsLinear(endStage,startStage) = mean(tempararyError,'omitnan'); 
            intervalAllXLinear{endStage,startStage} = 1./allX;
            sampleSizeMatLin(endStage,startStage) = length(allY);
        else
            monteTimeErrorsExtremes(endStage,startStage) = mean(tempararyError); 
            intervalAllXExtremes{endStage,startStage} = 1./allX;
            sampleSizeMatExt(endStage,startStage) = length(allY);
        end
        %Build BIC comparision variables, comparing models, then run
        %through BIC program
        if startStage~=endStage
        %store data points for AIC
        %all points as replicates
        	
        	AICTable.numPoints = length(allY);
            AICTable.yi = num2cell(allY,1);
            AICTable.xi = num2cell(allX,1);
                     	
        %AIC Comparison calculations, n=1 is assumed best fit (row stages, 
    	%columns are poly compared to) 
          	stage = 1;
                
          	%points as replicates
          	
            	polyNMdl1 = polyfit(allX, allY, 1);
             	polyNMdl2 = polyfit(allX, allY, 2);
                
           	BICComps = fitBICCompCalc(AICTable,polyNMdl1,polyNMdl2,stage);
   
        %Store BIC comp values and Eas (just in case)
            preDec = split(num2str(abs(log(BICComps))),'.');
            BICMat(endStage,startStage) = log(BICComps);
        %{
        if length(preDec{1}) >= 2
            BICMat(endStage,startStage) = round(log(BICComps),0);
        else
            BICMat(endStage,startStage) = round(log(BICComps),1);
        end
            %}
            EaMat(endStage,startStage) = (polyNMdl1(1)*-8.3144598)/1000; 
                %adjust slope to Ea with convertion Ea = -slope*R(in kJ)
        if p == 2 %2 extreme, 1 linear
        %Check concavity        
            isConcaveBioFly(endStage,startStage) = polyNMdl2(1) > 0; 
            isConvexBioFly(endStage,startStage) = polyNMdl2(1) < 0; 
            isLinearBioFly(endStage,startStage) = polyNMdl2(1) == 0; 
        end
        else
        end
    end
end

%plot the relevent BIC data range and the associated stage abbreviations
    fh = figure('position',[100 100 1000 800]);
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    %exp space
        h = heatmap(xlabels, ylabels, round(BICMat(2:12,1:11),2,'significant'),'CellLabelColor','black');

%Aestetics
    caxis([-5 5])
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    %colormap
        topx = flip(0.0:0.001:1);
        topmap = horzcat(ones(length(topx),1),(1*(topx)).',(1*(topx)).');
        botx = flip(0.0:0.001:1);
        botmap = horzcat((1*flip(botx)).',(1*flip(botx)).',ones(length(botx),1));
        map = vertcat(botmap,topmap);
        colormap(h,map)
        
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Start score code}')
    ylabel({'\fontsize{38}\bf{End score code}'})
    
    %determine temperatures used and adjust title accordingly
        annotation('textbox',[0.950, 0.942, 0.00005, 0.00005], 'string', '+','FontSize',32,'LineWidth',0.01)
        annotation('textbox',[0.964, 0.1435, 0.00005, 0.00005], 'string', '-','FontSize',32,'LineWidth',0.01)
            %reduce textbox size to near zero to get rid of it
    if p == 1
        title('\fontsize{40}Fly ln(L_{Q}/L_{L}): viable temp range');
    else
        title('\fontsize{40}Fly ln(L_{Q}/L_{L})');
    end
    
    
    %adjust color bar by covering up
        hmp = h.Position; 
        cbax = axes('Position',[sum(hmp([1,3])), 0, 1-sum(hmp([1,3])), 1],...
            'XTick',[], 'YTick', [], 'Color',[1 1 1]);
        cbax.XAxis.Visible = 'off';
        cbax.YAxis.Visible = 'off';
    % Set the new axis color map to match the 
    % heatmap's colormap
        cbax.Colormap = h.Colormap; 
    % Add a colorbar the same vertical position as the heatmap
        cbh = colorbar(cbax,'Position',[.90, hmp(2), .03, hmp(4)],'AxisLocation','in'); 
    % Set the limits to 0:1 and set the ticks 
        cbh.Limits = [0,1]; 
        nColors = size(h.Colormap,1); 
        cbh.Ticks = (0 : 1/10: 1); 
    % Set the new labels for each tick
        cbh.TickLabels = [-5:5]; 
    % Set the colorbar fontsize to the same value as heatmap fontsize
        cbh.FontSize = h.FontSize;
end
 
%%
%Figure 4A
%combining In silico and biological
%run only after running stage-wise code-chunk
figure('Position',[100 100 675 500]);

%remove T0 time point from analysis
    obsScores = (3:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    
%select stage of interest
    stage = 6;

%calculate and plot silico prediction first

    %parse constants and data from stage figure
        rangeT = [Temps(2) Temps(length(Temps)-2)];
            %only pull the temperatures we fit above
        xmin = 1/(max(Temps)+1);
        xmax = 1/(min(Temps)-1);
        xrange = [xmin; xmax];
        Ea = stageActEa(3:stage+2).'; %+2 cause of -b and a
            %only pull the Eas of stages we will predict on
        lnA = log(stageActA(3:stage+2)).';
        
    %set up fuction to predict data
        resolutionT = 160;
        temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';
        fcn = writelnKFcnVectorized(Ea,lnA);
        
    %calculate predicted data and plot
        lnKs = double(fcn(temps));
        
    %fit silico data and plot
        mdl = polyfit(1./temps,lnKs,1);
        xrangeHighRes = (xrange(1):range(xrange)/(resolutionT-1):xrange(2)).';
        fitEval = polyval(mdl,xrangeHighRes);
        %actual plotting done below
    
    %parse biological experiment data
        pseudoKperStage = cell(length(accumulatedInterval),1);
        for temp = 1:length(accumulatedInterval)
            pseudoKperStage{temp} = 1./accumulatedInterval{temp};
        end
        invertedTemps = 1./Temps;
        xmin = 1/(max(Temps)+1);
        xmax = 1/(min(Temps)-1);
        xrange = [xmin; xmax];
        
    %convert to ln space
        lnPseudoKperStage = cell(length(pseudoKperStage),1);
        for temp = 1:length(pseudoKperStage)
            lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
        end
        
    excludedTemps = [1 length(lnPseudoKperStage)-1 length(lnPseudoKperStage)];
    viableTemps = [2:length(lnPseudoKperStage)-1];
    
    %____________________________
    fittabeEmNum = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                end
            end
            overAllEms = overAllEms + sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
        end
    
    %initialize arrays to hold all fittable points    
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);  
        emNum = 1;

        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
    %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1); 
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
            standErr(temp,stage) = yStdErr;
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
        
        %plot the errorbars and mean
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                viable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end
        %{
        %build all data for specific temp for linear regression, while
        %plotting error bar
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&...
                    sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %fit to individuals
                eb = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'bo','LineWidth',2.5);
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,xneg,...
                    xpos,'r*','LineWidth',2.5);
            end 
            %}
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
                
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);
    %____________________________
    %{
    %calc total num of embryos
        fittableMean = 0;
        for temp = 1:length(lnPseudoKperStage)
            if (temp ~= 1)&&(temp ~= length(lnPseudoKperStage)-1)&&...
                    (temp ~= length(lnPseudoKperStage))
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 1
                    fittableMean = fittableMean + 1;
                end
            end
        end
        
    %initialize empty vectors to hold all data for fitting
    	allMeans = nan(fittableMean,1);
    	allXM = nan(fittableMean,1);
     	mNum = 1;
        
        sampSize = nan(length(invertedTemps),1);
    %replot the experimental results
        for temp = 1:length(invertedTemps)
            hold on;
            data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            sampSize(temp) = length(data);
      %set/calculate errors
            N = size(data,1); 
                %number experiments in data
            yStdErr = std(data)/sqrt(N); 
                %standard Error
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5;
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
        %plot error bar
            meanLnPseudoK = mean(data);
       	%build all data for specific temp for linear regression
            if (temp ~= 1)&&(temp ~= length(lnPseudoKperStage)-1)&&...
                    (temp ~= length(lnPseudoKperStage))&&length(data)>1
                eb = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'bo','LineWidth',2.5);
                allMeans(mNum) = meanLnPseudoK;
                allXM(mNum) = invertedTemps(temp);
                mNum = mNum + 1;
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,xneg,...
                    xpos,'r*','LineWidth',2.5);
            end 
         end
      %}  
    %calculate and plot biological linear regression
     	%mdl2 = polyfit(allXM, allMeans,1);
        mdl2 = polyfit(allX, allY,1);
      	yLin = polyval(mdl2,xrange);
      	lm = plot(xrange,yLin, 'k', 'LineWidth',2.5);
        lF = plot(xrangeHighRes,fitEval,'--c','LineWidth',2.5); 
        
    %Eas
        text(1./(273.15+16.5),-9.95,strcat("E_a = ",num2str(round(...
            (-mdl2(1)*R)*tokJConvert,0))," kJ/mol"),'Color','k','FontSize',24)
        text(1./(273.15+16.5),-10.25,strcat("E_a = ",num2str(round(...
            (-mdl(1)*R)*tokJConvert,0))," kJ/mol"),'Color','c','FontSize',24)
        
    %beautify
        ylim([mean(allY)-2 mean(allY)+1.25])
        ytLabels = yticks;
        yticks(ytLabels(floor(ytLabels)==ytLabels));
            %present only the interger yticks
        ylabel('ln(1/time [s^{-1}])');
        xlabel('1/(T [^oC] + 273)');
        %xticks(1./[303.1500 293.1500 283.1500]);
        xticks(1./[303.1500 294.1500 285.1500]);
        xticklabels({'1/30','1/21','1/12'});
        xlim(xrange);
        %xticklabels(strcat('1/',strsplit(num2str(round((1./xticks-273.15),1)))));
        set(gca,'linewidth', 3, 'FontSize', 30);   
        
    %seperate to avoid being affected by gca
        %title({'Fly Biological vs Prediction';strcat("A to ",stageAbbs(stage+2))}, 'FontSize', 25) %+2 for same reason
        legend([lm lF],{'Biological Fit 14.3-27^oC','Predictive Fit'},'FontSize',24)

%%
%Figure 4C
%Comparing worst case to biological data to a linear fit
%run only after running stage-wise code-chunk
figure('Position',[100 100 800 600]);

%display biological data
%remove T0 time point from analysis
    obsScores = (3:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    
%select stage of interest
    stage = 6;

%calculate and plot silico prediction first
    yyaxis left 
    
    %parse constants and data from stage figure
        rangeT = [Temps(2) Temps(length(Temps)-2)];
            %only pull the temperatures we fit above
        xrange = 1./[307.15; 280.15];
        Ea = stageActEa(3:stage+2).'; %+2 cause of -b and a
            %only pull the Eas of stages we will predict on
        lnA = log(stageActA(3:stage+2)).';
        
    %set up fuction to predict data
        resolutionT = 160;
        temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';
        fcn = writelnKFcnVectorized(Ea,lnA);
        
    %calculate predicted data and plot
        lnKs = double(fcn(temps));
        
    %fit silico data and plot
        mdl = polyfit(1./temps,lnKs,1);
        fitEval = polyval(mdl,xrange);
        %actual plotting done below
    
    %parse biological experiment data
        pseudoKperStage = cell(length(accumulatedInterval),1);
        for temp = 1:length(accumulatedInterval)
            pseudoKperStage{temp} = 1./accumulatedInterval{temp};
        end
        invertedTemps = 1./Temps;
        xmin = 1/(max(Temps)+1);
        xmax = 1/(min(Temps)-1);
        xrange = [xmin; xmax];
        
    %convert to ln space
        lnPseudoKperStage = cell(length(pseudoKperStage),1);
        for temp = 1:length(pseudoKperStage)
            lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
        end
    
    excludedTemps = [1 length(lnPseudoKperStage)-1 length(lnPseudoKperStage)];
    viableTemps = [2:length(lnPseudoKperStage)-1];
    
    %_________________
    fittabeEmNum = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                end
            end
            overAllEms = overAllEms + sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
        end
    
    %initialize arrays to hold all fittable points    
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);  
        emNum = 1;

        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
    %iterate over temperaters to average embryos for each temp then plot
    for temp = 1:length(invertedTemps)
        hold on;
        
        %parse temperary 'data' of nans
        	data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
            
        %set/calculate errors
            N = size(data,1); 
                %number non-nan experiments in data
            yStdErr = std(data)/sqrt(N-1); 
            standErr(temp,stage) = yStdErr;
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5; 
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
            
       %plot the errorbars and mean   
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp == viableTemps) == 1&&length(data)>0
            %fit to individuals
                blue = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            else
                nonViable = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            end
        %build all data for specific temp for linear regression,
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == ...
                    1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            end
        %{
        %build all data for specific temp for linear regression, while
        %plotting error bar
            meanLnPseudoK = mean(data);
            if isnan(meanLnPseudoK)==1
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&...
                    sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %fit to individuals
                blue = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'bo','LineWidth',2.5);
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,xneg,...
                    xpos,'r*','LineWidth',2.5);
            end 
            %}
            %copy over overAllData
                nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);
    %_________________
    %{
    %calc total num of embryos
        fittableMean = 0;
        for temp = 1:length(lnPseudoKperStage)
            if (temp ~= 1)&&(temp ~= length(lnPseudoKperStage)-1)&&...
                    (temp ~= length(lnPseudoKperStage))
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 1
                    fittableMean = fittableMean + 1;
                end
            end
        end
        
    %initialize empty vectors to hold all data for fitting
    	allMeans = nan(fittableMean,1);
    	allXM = nan(fittableMean,1);
     	mNum = 1;
        
    %replot the experimental results
        for temp = 1:length(invertedTemps)
            hold on;
            data = lnPseudoKperStage{temp}(:,stage);
            data = data(~isnan(data));
      %set/calculate errors
            N = size(data,1); 
                %number experiments in data
            yStdErr = std(data)/sqrt(N); 
                %standard Error
            yci95 = repmat(yStdErr,2,1);
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5;
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
        %plot error bar
            meanLnPseudoK = mean(data);
       	%build all data for specific temp for linear regression
            if (temp ~= 1)&&(temp ~= length(lnPseudoKperStage)-1)&&...
                    (temp ~= length(lnPseudoKperStage))&&length(data)>1
                blue = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'bo','LineWidth',2.5);
                allMeans(mNum) = meanLnPseudoK;
                allXM(mNum) = invertedTemps(temp);
                mNum = mNum + 1;
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,xneg,...
                    xpos,'r*','LineWidth',2.5);
            end 
         end
    %}    
    %calculate and plot biological linear regression
     	mdl2 = polyfit(allX, allY,1);
      	yLin = polyval(mdl2,xrange);
      	lm = plot(xrange,yLin, 'k', 'LineWidth',2.5);
        
        ylim([mean(allY)-2 max(yLin)])
        %ylim([mean(allMeans)-2 mean(allMeans)+1.25])
        
        lowerYCorrection = (min(yLin) - (mean(allY)-2))/(max(yLin)-min(yLin));
        
        ylabs = cell2mat(cellfun(@str2num,yticklabels,'un',0));
                yticks(ylabs(ylabs == floor(ylabs)));
                yticklabels(ylabs(ylabs == floor(ylabs)));
                ylabel('ln(1/time [s^{-1}]) Fly data');
        
    %beautify

    %{
        ytLabels = yticks;
        yticks(ytLabels(floor(ytLabels)==ytLabels));
            %present only the interger yticks
        xlim(xrange);
        xlabel('1/(T [^oC] + 273.15)');
        xticks(1./[303.1500 293.1500 283.1500]);
        xlim(xrange);
        xticklabels(strcat('1/',strsplit(num2str(round((1./xticks-273.15),1)))));
        set(gca,'linewidth', 3, 'FontSize', 30);   
    %}
        

%worst case
    R = 8.3144598*(10^(-3));
    netSize = 1000;
    
    %Set temperature range and initialize temps
    	rangeT = [280.15 307.65];
    	resolutionT = 30;
    	temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';
        above =  temps<(Temps(length(Temps)-2));
        below =  temps<(Temps(2));
        fitInd = find(above~=below);
        fitTs = temps(fitInd);
        
    %worst case
        worstEa = 20 - log(netSize-1)*(R*295.15);
        Ea = [100; repmat(worstEa,netSize-1,1)];
        lnA = [37.2926197009076; repmat(4.79899926729709,netSize-1,1)];

        %plot "mean"
            %write function with NO noise
                fcnWorst = writelnKFcnVectorized(Ea,lnA);
            %calculate lnKs using temp range, resolution, and function(Ea,lnA)
                worstlnKs = double(fcnWorst(temps));

                yyaxis right 
                magenta = plot(1./temps,worstlnKs,'--','LineWidth',2.5);
                    hold on

                linFitWorst = polyfit(1./fitTs',double(fcnWorst(fitTs))',1);
                fitEval = polyval(linFitWorst,xrange);
                black = plot(xrange,fitEval,'k','LineWidth',2.5);
                    hold on

                ylabel('ln(k) worst case','FontSize',28)
                ylim([min(fitEval)-lowerYCorrection*(max(fitEval)-min(fitEval)) max(fitEval)]); 
        
        %aesthetics
            %control position so can focus on curvature
                xticks(1./([30 20 10]+273.15))

            %Conditional titles and legends
                %title({'Fly A - G vs 1000x Worst Case'},'FontSize',25)

            %legend
                ylabs = cell2mat(cellfun(@str2num,yticklabels,'un',0));
                yticks(ylabs(ylabs == floor(ylabs)));
                yticklabels(ylabs(ylabs == floor(ylabs)));
                xlabel('1/(T [^oC] + 273)','FontSize',28)
                xticks(1./[303.1500 294.1500 285.1500]);
                xticklabels({'1/30','1/21','1/12'});
                xlim(xrange)
                %xticklabels(strcat('1/',strsplit(num2str(round((1./xticks-273.15),1)))));
            %convert decimal xticks to 1/T in K
                set(gca,'linewidth', 3, 'FontSize', 28);
                legend([blue nonViable magenta black], {'','Biological Data: A - G','1000x Worst Case Model',strcat("Linear Fit ",num2str(Temps(2)-273.15)," - ", num2str(Temps(length(Temps)-2)-273.15)," ^oC")},'FontSize',24);    
%%
%Comparing old N and new N
%load in new N data first
%Data Import

%New
%Destination
    %filename = 'FlyScoresFinalCombined.xlsx';
    filename = 'FlyScoresFinal.xlsx';
    folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Fly Scores/';

%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    %[~,sheet_name]=xlsfinfo(strcat(filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{1});
    %data = readcell(strcat(filename),'Sheet', sheet_name{1});

%find number of events scored ("scores") in file
    isSearch = cellfun(@(x)isequal(x,'Score Number'),data);
    [row,col] = find(isSearch);
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+2;
        %+1 to include T0 as a score
    
%Initialize arrays to store different data types        
    accumulatedNew = cell(sizeData,1);
    perStageNew = cell(sizeData,1);

%iterate over each sheet to parse data types into their respective arrays        
for i = 1:sizeData
%read in specific data of sheet i
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{i});
    %data = readcell(strcat(filename),'Sheet', sheet_name{i});
    
%parse out "accumulated" data rows for each embryo of sheet i and store in 
%cell array
    isSearch = cellfun(@(x)isequal(x,'Accumulated'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize cell i to contain each accum data row from a single sheet
        accumulatedNew{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
        for embryo = 1:length(row)
        %iterate over all of the columns (or stages of interest) from sheet
            for column = 1:numScores
             %check if data is not missing, if so, copy over into matrix
                if ismissing(data{row(embryo),(col+1+column)})==0
                    accumulatedNew{i}(embryo,column) = cell2mat(data(row...
                        (embryo),(col+1+column)))*60; %convert to seconds
                else
                end
            end
        end
        
%parse out per stage data for each sheet in same manner
    isSearch = cellfun(@(x)isequal(x,'Per stage'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize a cell to contain each perStage data row from a single sheet
        perStageNew{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
    for embryo = 1:length(row)
    %iterate over all of the columns (or stages of interest) from sheet
        for column = 1:numScores
        %check is data is not missing, if so, copy over into matrix
            if ismissing(data{row(embryo),(col+1+column)})==0
                perStageNew{i}(embryo,column) = cell2mat(data(row...
                    (embryo),(col+1+column)))*60; %convert to seconds
            else
            end
        end
    end
end    
    

%Old    
%Destination
    filename = 'FlyScoresFinalNewScoring.xlsx';

   % filename = 'FlyScoresFinal.xlsx';
    folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Fly Scores/';

%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    %[~,sheet_name]=xlsfinfo(strcat(filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{1});
    %data = readcell(strcat(filename),'Sheet', sheet_name{1});

%find number of events scored ("scores") in file
    isSearch = cellfun(@(x)isequal(x,'Score Number'),data);
    [row,col] = find(isSearch);
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+2;
        %+1 to include T0 as a score
    
%Initialize arrays to store different data types        
    accumulatedOld = cell(sizeData,1);
    perStageOld = cell(sizeData,1);
    Temps = zeros(sizeData,1);
    stageNames = strings(1,numScores);
    stageAbbs = strings(1,numScores);

%iterate over each sheet to parse data types into their respective arrays        
for i = 1:sizeData
%read in specific data of sheet i
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{i});
    %data = readcell(strcat(filename),'Sheet', sheet_name{i});
    
%parse out "accumulated" data rows for each embryo of sheet i and store in 
%cell array
    isSearch = cellfun(@(x)isequal(x,'Accumulated'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize cell i to contain each accum data row from a single sheet
        accumulatedOld{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
        for embryo = 1:length(row)
        %iterate over all of the columns (or stages of interest) from sheet
            for column = 1:numScores
             %check if data is not missing, if so, copy over into matrix
                if ismissing(data{row(embryo),(col+1+column)})==0
                    accumulatedOld{i}(embryo,column) = cell2mat(data(row...
                        (embryo),(col+1+column)))*60; %convert to seconds
                else
                end
            end
        end
        
%parse out per stage data for each sheet in same manner
    isSearch = cellfun(@(x)isequal(x,'Per stage'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize a cell to contain each perStage data row from a single sheet
        perStageOld{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
    for embryo = 1:length(row)
    %iterate over all of the columns (or stages of interest) from sheet
        for column = 1:numScores
        %check is data is not missing, if so, copy over into matrix
            if ismissing(data{row(embryo),(col+1+column)})==0
                perStageOld{i}(embryo,column) = cell2mat(data(row...
                    (embryo),(col+1+column)))*60; %convert to seconds
            else
            end
        end
    end

%parse out temperature data for each sheet and store in array
    isSearch = cellfun(@(x)isequal(x,'Temp'),data);
    [row,col] = find(isSearch);
    Temps(i) = cell2mat(data(row,(col+1)));
end

%parse out the full names of every score, including T0
    isSearch = cellfun(@(x)isequal(x,'Score Name'),data);
    [row,col] = find(isSearch);
    for column = 1:numScores
        stageNames(column) = string(data(row,(col+1+column)));
    end
    
%parse out the abbreviations of every score
    isSearch = cellfun(@(x)isequal(x,'Abb.'),data);
    [row,col] = find(isSearch);
    for column = 1:numScores
        stageAbbs(column) = string(data(row,(col+1+column)));
    end    
%%    
%actual comparison of above data import
%Arrhenius Stage-wise
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;

%Plot Stages
    figure('Position',[10 10 1350 850]);
    %figure('Position',[10 10 1350 750]);
    %colors = distinguishable_colors(3);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    %perStageIntervalNew = cellfun(@(x) x(:,obsScores),perStageNew,'UniformOutput',false);
    %perStageIntervalOld = cellfun(@(x) x(:,obsScores),perStageOld,'UniformOutput',false);
    
    obsScores = (1:numScores);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 sizeData-1 sizeData];
    
%Initialize tables
    %Statistical Table
        compRegs = zeros(2,8); %2 cause comparing 2 data sets for each stage
        compRegs = array2table(compRegs);
        compRegs.Properties.VariableNames = {'numPoints','SX','SY',...
            'SXX','SYY','SXY','AveX','AveY'};
        compRegs.Properties.RowNames = ["New" "Old"];

%invert temps and temp range        
    invertedTemps = 1./Temps;
    xmin = 1/307.15;
    xmax = 1/280.15;
    xrange = [xmin; xmax];

%convert raw data to log(inverse(t))
    lnPseudoKperStageNew = cellfun(@(x) real(log(1./x)),perStageNew,'UniformOutput',false);
    lnPseudoKperStageOld = cellfun(@(x) real(log(1./x)),perStageOld,'UniformOutput',false);
%{    
%convert to pseudo Ks and inverse temps
    pseudoKperStage = cell(length(perStageInterval),1);
    for temp = 1:length(perStageInterval)
        pseudoKperStage{temp} = 1./abs(perStageInterval{temp});
    end
    
%convert to ln space
    lnPseudoKperStage = cell(length(pseudoKperStage),1);
    for temp = 1:length(pseudoKperStage)
        lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
    end
%}    
%iterate over stages to plot each stage
%initialize variable to hold error Statistic
    standErr = nan(length(invertedTemps),length(lnPseudoKperStageNew{1}));

%for all use -> length(lnPseudoKperStage{1})
    [numRow,numCol] = calcRectangle(length(lnPseudoKperStageNew{1}));
    %[numRow,numCol] = calcRectangle(6);
    pannelCount = 1;
for stage = [1,3:length(lnPseudoKperStageNew{1})]
    %initialize and store subplot data
        sp(pannelCount) = subplot(numRow,numCol,pannelCount);  
    
    %Plot new N data (want new ploted first so old shows up overlapping)
        %calc total num of fittable embryos per stage for all temps, 2 or more
        %reps
            fittabeEmNum = 0;
            fittableMean = 0;
            overAllEms = 0;
            for temp = 1:length(lnPseudoKperStageNew)
                if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                    if sum(~isnan(lnPseudoKperStageNew{temp}(:,stage))) > 0% 1
                        fittabeEmNum = fittabeEmNum + ...
                            sum(~isnan(lnPseudoKperStageNew{temp}(:,stage)));
                        fittableMean = fittableMean + 1;                    
                    end
                end
                overAllEms = overAllEms + sum(~isnan(lnPseudoKperStageNew{temp}(:,stage)));
            end 

        %initialize arrays to hold all fittable points    
            allX = nan(fittabeEmNum,1);
            allY = nan(fittabeEmNum,1);  
            emNum = 1;

        %iterate over temperaters to average embryos for each temp then plot
        for temp = 1:length(invertedTemps)
            hold on;

            %parse temperary 'data' of nans
                data = lnPseudoKperStageNew{temp}(:,stage);
                data = data(~isnan(data));

            %set/calculate errors
                N = size(data,1); 
                    %number non-nan experiments in data
                yStdErr = std(data)/sqrt(N-1); 
                standErr(temp,stage) = yStdErr;
                    %standard Error
                ci95 = tinv([0.025 0.975], N-1); 
                    %95 confidence interval
                yci95 = yStdErr.*ci95;
                yneg = yci95(1);
                ypos = yci95(2);
                errTemp = 1.96*0.5; 
                    %1.96 is multiplier associated with a two sided confidence interval
                xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
                xpos = xneg;

            %build all data for specific temp for linear regression, while
            %plotting error bar
                meanLnPseudoK = mean(data);
                if isnan(meanLnPseudoK)==1
                elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0%1
                    nextEmNum = emNum+size(data,1);
                    allY(emNum:nextEmNum-1,1) = data;
                    allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                        size(data,1),1);
                    emNum = nextEmNum;
                %fit to individuals
                    fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                        ypos,xneg,xpos,'o','Color',[0 1 1],'LineWidth',2.5);
                elseif length(data)==1
                    errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                        xneg,xpos,'*','Color',[1 0.5 0],'LineWidth',2.5);                
                else
                    errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                        xneg,xpos,'*','Color',[1 0.5 0],'LineWidth',2.5);
                end

        end
        %combine allX/Y
            allData = horzcat(allX,allY);

        %store regression stats
            %individual points

                compRegs.numPoints(1) = length(allData);
                compRegs.SX(1) = sum(allData(:,1));
                compRegs.SY(1) = sum(allData(:,2));
                compRegs.SXX(1) = sum(allData(:,1).^2);
                compRegs.SYY(1) = sum(allData(:,2).^2);
                compRegs.SXY(1) = sum(allData(:,1).*allData(:,2));
                compRegs.AveX(1) = mean(allData(:,1));
                compRegs.AveY(1) = mean(allData(:,2));

        %plot linear regression
            lineFit = fit(allData(:,1), allData(:,2), 'poly1');
            Ea1 = round((-lineFit.p1*R)*tokJConvert,0);

        %plot fit lines
            yLin = polyval([lineFit.p1 lineFit.p2],xrange);
            lm = plot(xrange,yLin, 'c', 'LineWidth',2.5);


    %Plot old data
    %Plot new N data (want new ploted first so old shows up overlapping)
        %calc total num of fittable embryos per stage for all temps, 2 or more
        %reps
            fittabeEmNum = 0;
            fittableMean = 0;
            overAllEms = 0;
            for temp = 1:length(lnPseudoKperStageOld)
                if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                    if sum(~isnan(lnPseudoKperStageOld{temp}(:,stage))) > 0%1
                        fittabeEmNum = fittabeEmNum + ...
                            sum(~isnan(lnPseudoKperStageOld{temp}(:,stage)));
                        fittableMean = fittableMean + 1;                    
                    end
                end
                overAllEms = overAllEms + sum(~isnan(lnPseudoKperStageOld{temp}(:,stage)));
            end 

        %initialize arrays to hold all fittable points    
            allX = nan(fittabeEmNum,1);
            allY = nan(fittabeEmNum,1);  
            emNum = 1;

        %iterate over temperaters to average embryos for each temp then plot
        for temp = 1:length(invertedTemps)
            hold on;

            %parse temperary 'data' of nans
                data = lnPseudoKperStageOld{temp}(:,stage);
                data = data(~isnan(data));

            %set/calculate errors
                N = size(data,1); 
                    %number non-nan experiments in data
                yStdErr = std(data)/sqrt(N-1); 
                standErr(temp,stage) = yStdErr;
                    %standard Error
                ci95 = tinv([0.025 0.975], N-1); 
                    %95 confidence interval
                yci95 = yStdErr.*ci95;
                yneg = yci95(1);
                ypos = yci95(2);
                errTemp = 1.96*0.5; 
                    %1.96 is multiplier associated with a two sided confidence interval
                xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
                xpos = xneg;

            %build all data for specific temp for linear regression, while
            %plotting error bar
                meanLnPseudoK = mean(data);
                if isnan(meanLnPseudoK)==1
                elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0%1
                    nextEmNum = emNum+size(data,1);
                    allY(emNum:nextEmNum-1,1) = data;
                    allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                        size(data,1),1);
                    emNum = nextEmNum;
                %fit to individuals
                    fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                        ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
                elseif length(data)==1
                    errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                        xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
                else
                    errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                        xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);
                end

        end
        %combine allX/Y
            allData = horzcat(allX,allY);

        %store regression stats
            %individual points

                compRegs.numPoints(2) = length(allData);
                compRegs.SX(2) = sum(allData(:,1));
                compRegs.SY(2) = sum(allData(:,2));
                compRegs.SXX(2) = sum(allData(:,1).^2);
                compRegs.SYY(2) = sum(allData(:,2).^2);
                compRegs.SXY(2) = sum(allData(:,1).*allData(:,2));
                compRegs.AveX(2) = mean(allData(:,1));
                compRegs.AveY(2) = mean(allData(:,2));

        %plot linear regression
            lineFit = fit(allData(:,1), allData(:,2), 'poly1');
            Ea2 = round((-lineFit.p1*R)*tokJConvert,0);

        %plot fit lines
            yLin = polyval([lineFit.p1 lineFit.p2],xrange);
            lm = plot(xrange,yLin, 'b', 'LineWidth',2.5);
               
 
        
    %perform statistical analysis            
        pVal = ancova(compRegs);
        text(1./(273.15+30),-5.25,strcat("E_a = "),'Color','k','FontSize',19)
        text(1./(273.15+24),-5.1,strcat(num2str(Ea1)),'Color','c','FontSize',19)
        text(1./(273.15+21),-5.1,strcat(",",num2str(Ea2)),'Color','b','FontSize',19)
        text(1./(273.15+18),-5.1,strcat(" kJ/mol"),'Color','k','FontSize',19)
        text(1./(273.15+19.5),-6,strcat("pVal = ",num2str(round(pVal,3))),'Color','k','FontSize',19)
        
        
        
    %beautify
        %Set consistant limits to compare
            ylim([-13 -4.5]) 
                %for all plots, to compare
                %ylim([-13.5 -3])
            %yticks([-12 -10 -8 -6])
                yticks([-12 -10 -8 -6 -4])
            xlim(xrange);

        %adjust xticks
            xticks(1./[299.1500 285.1500]);
            xticklabels({'1/26','1/12'});

        %set labels on outer most subplots, remove from interior        
            if stage==1||stage==6||stage==10%stage==3||stage==6%
                ylabel('ln(1/time [s^{-1})]');
            else
                set(gca,'yticklabel',[])
            end
            if stage>=9 && stage<=12%stage>=6 && stage<=8%
                xlabel('1/(T [^oC] + 273.15)');
            else
                set(gca,'xticklabel',[])
            end

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 23);
                %28

        %place titles as text inside plots
            if stage == 1
                text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(stage),...
                    " to ",(stageAbbsNoT0(stage+1))))},...
                    'FontSize', 24,'HorizontalAlignment','center');
                    %'FontSize', 28,'HorizontalAlignment','center');
            else
                text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(stage-1),...
                    " to ",(stageAbbsNoT0(stage))))}, 'FontSize', 24,...
                    'HorizontalAlignment','center');
                    %" to ",(stageAbbsNoT0(stage))))}, 'FontSize', 28,...
                    %'HorizontalAlignment','center');
            end
            
    pannelCount = pannelCount + 1;
end



%Adjust Subplot positions to tighten image
    %Shift left  
        for stage = 1:length(obsScores)-1
            if stage ~= 1 && stage ~= 5 && stage ~= 9
                prevPos = get(sp(stage-1),'Position');
                prevPos = [prevPos(1)+prevPos(3)+0.0075 prevPos(2) ...
                    prevPos(3) prevPos(4)];
                set(sp(stage),'Position', prevPos)
            end
        end

    %shift up
        for stage = 1:length(obsScores)-1
            
            if stage >= 5 && stage <=8
                prevPos = get(sp(1),'Position');
                curPos = get(sp(stage),'Position');
                curPos(2) = prevPos(2)-prevPos(4)-0.015;
                set(sp(stage),'Position', curPos)
            end
            if stage >= 9 && stage <=12
                prevPos = get(sp(5),'Position');
                curPos = get(sp(stage),'Position');
                curPos(2) = prevPos(2)-prevPos(4)-0.015;
                set(sp(stage),'Position', curPos)
            end
        end
        
%%
%All CV Analyis
%Data Import (rerun actual data after this because this script overides the
%good data loaded in before)

%Destination
    filename = 'UnfliteredFlyScoresFinal.xlsx';
    folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Fly Scores/';

%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{1});

%find number of events scored ("scores") in file
    isSearch = cellfun(@(x)isequal(x,'Score Number'),data);
    [row,col] = find(isSearch);
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+2;
        %+1 to include T0 as a score
    
%Initialize arrays to store different data types        
    accumulated = cell(sizeData,1);
    perStage = cell(sizeData,1);
    Temps = zeros(sizeData,1);
    stageNames = strings(1,numScores);
    stageAbbs = strings(1,numScores);

%iterate over each sheet to parse data types into their respective arrays        
for i = 1:sizeData
%read in specific data of sheet i
    data = readcell(strcat(folder,filename),'Sheet', sheet_name{i});
    
%parse out "accumulated" data rows for each embryo of sheet i and store in 
%cell array
    isSearch = cellfun(@(x)isequal(x,'Accumulated'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize cell i to contain each accum data row from a single sheet
        accumulated{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
        for embryo = 1:length(row)
        %iterate over all of the columns (or stages of interest) from sheet
            for column = 1:numScores
             %check if data is not missing, if so, copy over into matrix
                if ismissing(data{row(embryo),(col+1+column)})==0
                    accumulated{i}(embryo,column) = cell2mat(data(row...
                        (embryo),(col+1+column)))*60; %convert to seconds
                else
                end
            end
        end
        
%parse out per stage data for each sheet in same manner
    isSearch = cellfun(@(x)isequal(x,'Per stage'),data);
    [row,col] = find(isSearch);
    col = unique(col);
    
    %initialize a cell to contain each perStage data row from a single sheet
        perStage{i} = nan(length(row),numScores);
        
    %iterate over number of embryos (rows) in specific sheet
    for embryo = 1:length(row)
    %iterate over all of the columns (or stages of interest) from sheet
        for column = 1:numScores
        %check is data is not missing, if so, copy over into matrix
            if ismissing(data{row(embryo),(col+1+column)})==0
                perStage{i}(embryo,column) = cell2mat(data(row...
                    (embryo),(col+1+column)))*60; %convert to seconds
            else
            end
        end
    end

%parse out temperature data for each sheet and store in array
    isSearch = cellfun(@(x)isequal(x,'Temp'),data);
    [row,col] = find(isSearch);
    Temps(i) = cell2mat(data(row,(col+1)));
end

%parse out the full names of every score, including T0
    isSearch = cellfun(@(x)isequal(x,'Score Name'),data);
    [row,col] = find(isSearch);
    for column = 1:numScores
        stageNames(column) = string(data(row,(col+1+column)));
    end
    
%parse out the abbreviations of every score
    isSearch = cellfun(@(x)isequal(x,'Abb.'),data);
    [row,col] = find(isSearch);
    for column = 1:numScores
        stageAbbs(column) = string(data(row,(col+1+column)));
    end


%Build CV heatmap
    %Grid of stages by stages
    %each box represents the start (xaxis) and end point (yaxis)
    %value shows the CV (%) for that duration averaged over all temperatures
fh = figure('position',[100 100 1000 800])

%initialize matrix to hold presentation data
    CVMat = nan(numScores);
    tempSampleSizeMat = nan(numScores);
    embryoSampleSizeMat = nan(floor((numScores^2)/2),2);
    
%iterate over ever (forward) stage interval to calculate co-variance
counter = 1;
for startStage = 1:numScores
%only necessary to look at forward intervals, so end stage begins at
%current start stage
    for endStage = startStage:numScores
    %calculate interval times by diff between acc(stages) times
        dataRange = cellfun(@(x) x(:,[startStage endStage]),accumulated,'UniformOutput',false);
        diffDataRange = cellfun(@(x) diff(x,1,2),dataRange,'UniformOutput',false);
    %convert any zeros to NaNs so as to avoid messing up mean calculation 
    %(should be no zeros)    
        convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'[^0123456789.]0','NaN')).',diffDataRange,'UniformOutput',false);
        convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'^0[^0123456789.]','NaN')).',convertZeros,'UniformOutput',false);        
    %calculate mean over each embryo(row) at each stage(loop) for all temperatures(cells)   
        meanDataPerTemp = cell2mat(cellfun(@(x) mean(x,1,'omitnan'),convertZeros,'UniformOutput',false));
    %same for std
        stdDataPerTemp = cell2mat(cellfun(@(x) std(x,1,'omitnan'),convertZeros,'UniformOutput',false));
    %calculate all CVs for this interval over all temps
        CVPerTemp = stdDataPerTemp./meanDataPerTemp;
    %copy into matrix
        CVMat(endStage,startStage) = mean(CVPerTemp,'omitnan')*100;
        tempSampleSizeMat(endStage,startStage) = length(convertZeros);
        embryoSampleSizeMat(counter,:) = [min(cell2mat(...
            cellfun(@(x) length(x),convertZeros,'UniformOutput',false)))...
            max(cell2mat(cellfun(@(x) length(x),convertZeros,...
            'UniformOutput',false)))];
        
        counter = counter + 1;
    end
end

%determine heatmap presentation area, slightly limited here
ylabels = stageAbbs(2:numScores);
xlabels = stageAbbs(1:numScores-1);
h = heatmap(xlabels, ylabels, round(CVMat(2:numScores,1:numScores-1),2,'significant'),'CellLabelColor','none')

%Aestetics
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    caxis([0 31])

    x = [0:0.001:1];
    map = horzcat((1*flip(x.^3)).',1-flip(x.^3).', zeros(length(x),1));
        %non-linear gradient
    topx = flip(0.0:0.001:1);
    topmap = horzcat(ones(length(topx),1),zeros(length(topx),1),zeros(length(topx),1));
    botx = flip(0.0:0.001:1);
    botmap = horzcat((1*flip(botx.^2)).',ones(length(topx),1),(1*flip(botx)).');
    map = vertcat(botmap,topmap);
        
    colormap(h,map)
    set(gca,'FontSize', 32);
    title('\fontsize{38}Mean fly interval CV');
    xlabel('\fontsize{38}\bf{Start score code}')
    ylabel({'\fontsize{38}\bf{End score code}'})
    
    %adjust color bar by covering up
        hmp = h.Position; 
        cbax = axes('Position',[sum(hmp([1,3])), 0, 1-sum(hmp([1,3])), 1],...
            'XTick',[], 'YTick', [], 'Color',[1 1 1]);
        cbax.XAxis.Visible = 'off';
        cbax.YAxis.Visible = 'off';
    % Set the new axis color map to match the 
    % heatmap's colormap
        cbax.Colormap = h.Colormap; 
    % Add a colorbar the same vertical position as the heatmap
        cbh = colorbar(cbax,'Position',[.90, hmp(2), .03, hmp(4)],'AxisLocation','in'); 
    % Set the limits to 0:1 and set the ticks 
        cbh.Limits = [0,1]; 
        nColors = size(h.Colormap,1); 
        cbh.Ticks = (0:1/2:1); 
    % Set the new labels for each tick
        cbh.TickLabels = [0:31/2:31]; 
    % Set the colorbar fontsize to the same value as heatmap fontsize
        cbh.FontSize = h.FontSize;
        
        
    IPx = h.InnerPosition(1);
    IPy = h.InnerPosition(2);
    textAxis = axes('Position',[0 0 1 1],'visible','off');
    scaleX = h.InnerPosition(3)/(numScores-1);
    scaleY = h.InnerPosition(4)/(numScores-1);
    
%
for x = 1:numScores-1
    for y = x:numScores-1
        txt = num2str(round(CVMat(y+1,x),2,'significant'));
        if str2num(txt) < 10 && contains(txt,'.') == 0
            txt = strcat(txt,'.0');
        else
        end
        xShift = IPx - 0.5*scaleX;
        yShift = IPy + h.InnerPosition(4) + 0.5*scaleY;
        text(xShift + scaleX*x,yShift - scaleY*y,txt,'Color','k','FontSize',21,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
end
%} 