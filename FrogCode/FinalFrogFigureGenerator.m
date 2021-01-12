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
    %filename = 'FrogScoresFinal_4.xlsx';
    filename = 'FrogScoresFinalCombinedData.xlsx';
    %filename = 'FrogScoresFinalNewData.xlsx';
    folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Frog Scores/';

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
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+1;
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
    %data = readcell(strcat(filename),'Sheet', sheet_name{i});
    
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
%Raw Data Presentation

%generate as seperate figures

%array of temperatures to show (only being selective when there
%are temperatures very close together)
    obsTemps = [1,2,3,5,8,9,10,12,13,16];

%Figure 1, all data
    figure('Position',[100 100 1400 650])

%select what scores to look at, where 1 is T0
    obsScores = (1:12);
        %drop muscle twitching,13
    colors = distinguishable_colors(numScores);

%initialize display data arrays for means and errors
    meanAccByStage = zeros(sum(length(obsTemps)),length(obsScores));
    err = zeros(sum(length(obsTemps)),length(obsScores));

%iterate over every observable temperature to condense multiple embryos in
%each score to a mean and calculate errors
for temp = 1:length(obsTemps)
%parse out score data for specific temp
    data = accumulated{obsTemps(temp)}(:,obsScores);
%determine number of measurable samples for each score
    N = sum(~isnan(data),1);
%take mean of each score (no NaNs) then calculate 95% CI and store
    meanAccByStage(temp,:) = mean(data,1,'omitnan');
    yStdErr = std(data,1,'omitnan')./sqrt(N); 
        %standard Error
    ci95 = tinv(0.025, N-1); 
        %95 confidence interval
    yci95 = yStdErr.*ci95;
    err(temp,:) = yci95;
%calculate 95% CI for temeprature
	errTemp = 1.96*0.5;
    xneg(temp,:) = Temps(obsTemps(temp))-(Temps(obsTemps(temp))+errTemp);
    xpos = xneg;
end


%iterate over every score to to plot means of every temp
for score = 1:length(obsScores)
    hold on;
%plot basic for T0, specific with errors for T1 on
    if score == 1
        errorbar(zeros(length(obsTemps),1),Temps(obsTemps),xneg,...
            'vertical','-o','Color',colors(obsScores(score),:),'LineWidth',2.5);
    else
    %check to see is any NaNs, index by not NaN. Make temporary variable to
    %index through
        idx = ~any(isnan(meanAccByStage(:,score)),2);
        devTime = meanAccByStage(:,score);
        errTime = err(:,score);
        t = Temps(obsTemps);       
        errorbar(devTime(idx),t(idx),xneg(idx),xpos(idx),errTime(idx),...
            errTime(idx),'-o','Color',colors(obsScores(score),:),'LineWidth',2.5);
    end
end

%generate box around zoomed area
    plot([-40,-40],[(max(Temps)+1.25),(min(Temps)-1.25)],'LineWidth',2,'Color','k');
    plot([-40,181.3],[(max(Temps)+1.25),(max(Temps)+1.25)],'LineWidth',2,'Color','k');
    plot([800,181.3],[(min(Temps)-1.25),(max(Temps)+1.25)],'LineWidth',2,'Color','k');
    plot([800,-40],[(min(Temps)-1.25),(min(Temps)-1.25)],'LineWidth',2,'Color','k');

%Aestetics
    title('Frog stages developmental times over range of temperatures')
    ylim([(min(Temps)-1.5) (max(Temps)+1.5)]);
    ylabels = str2num(cell2mat(yticklabels));
    yticks(ylabels([1 3 5 7 9 11])+0.15);
    yticklabels(cellstr(num2str((ylabels([1 3 5 7 9 11])+0.15)-273.15)));
    ylabel('Temperature (C)');
    xlim([-120 3750])
    xlabel('Time post 3rd cleav (mins)');
    displayNames = stageAbbs(obsScores);
    legend(cellstr(displayNames),'NumColumns',2);
    set(gca,'LineWidth', 2.5, 'FontSize', 33);
    
%Figure 2, replot data from figure 1, but zoomed in
    figure('Position',[100 100 1400 650])

%select what scores to look at, where 1 is T0
    obsScores = (1:8);
    colors = distinguishable_colors(numScores);

%iterate over every score to to plot means of every temp
for score = 1:length(obsScores)
    hold on;
%plot basic for T0, specific with errors for T1 on
    if score == 1
        errorbar(zeros(length(obsTemps),1),Temps(obsTemps),xneg,...
            'vertical','-o','Color',colors(obsScores(score),:),'LineWidth',2.5);
    else
    %check to see is any NaNs, index by not NaN. Make temporary variable to
    %index through
        idx = ~any(isnan(meanAccByStage(:,score)),2);
        devTime = meanAccByStage(:,score);
        errTime = err(:,score);
        t = Temps(obsTemps);       
        errorbar(devTime(idx),t(idx),xneg(idx),xpos(idx),errTime(idx),...
            errTime(idx),'-o','Color',colors(obsScores(score),:),'LineWidth',2.5);
    end
end

%generate box around zoomed area
    plot([-10,-10],[(max(Temps)+1.25),(min(Temps)-1.25)],'LineWidth',2,'Color','k');
    plot([-10,151.3],[(max(Temps)+1.25),(max(Temps)+1.25)],'LineWidth',2,'Color','k');
    plot([780,151.3],[(min(Temps)-1.25),(max(Temps)+1.25)],'LineWidth',2,'Color','k');
    plot([780,-10],[(min(Temps)-1.25),(min(Temps)-1.25)],'LineWidth',2,'Color','k');
    
%Aestetics
    ylim([(min(Temps)-1.5) (max(Temps)+1.5)]);
    ylim([(min(Temps)-1.5) (max(Temps)+1.5)]);
    ylabels = str2num(cell2mat(yticklabels));
    yticks(ylabels([1 3 5 7 9 11])+0.15);
    yticklabels(cellstr(num2str((ylabels([1 3 5 7 9 11])+0.15)-273.15)));
    ylabel('Temperature (C)');
    xlim([-20 750])
    xlabel('Time post 3rd cleav (mins)');
    displayNames = stageAbbs(obsScores);
    legend(cellstr(displayNames),'NumColumns',2);
    set(gca,'linewidth', 2.5, 'FontSize', 33);
%%
%Figure 1 D, non-log version
%as above, but using a personalized errorbar function to allow for
%transparency
%plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
%Raw Data Presentation

%generate as seperate figures

%array of temperatures to show (only being selective when there
%are temperatures very close together)
    obsTemps = [1:10 12:13 15:length(Temps)];%[1:length(accumulated)];%
        %for combined data [1:10 12:13 17:length(Temps)]; 11 is new data,
        %very bad
%Figure 1, all data
    figureWidth = 1600;
    figureHeight = 675;
    figure('Position',[100 100 figureWidth figureHeight])
    aspectRatio = figureWidth/figureHeight;
%select what scores to look at, where 1 is T0
    obsScores = (1:12);
        %drop muscle twitching,13
    colors = distinguishable_colors(numScores);
    RGBa = horzcat(colors,zeros(numScores,1));
   
%initialize display data arrays for means and errors
    meanAccByStage = zeros(sum(length(obsTemps)),length(obsScores));
    err = zeros(sum(length(obsTemps)),length(obsScores));

%iterate over every observable temperature to condense multiple embryos in
%each score to a mean and calculate errors
    xneg = nan(length(obsTemps),1);
    %plottedTemps = %
for temp = 1:length(obsTemps)%
%parse out score data for specific temp
    data = accumulated{obsTemps(temp)}(:,obsScores)/60;
%determine number of measurable samples for each score
    N = sum(~isnan(data),1);
%take mean of each score (no NaNs) then calculate 95% CI and store
    meanAccByStage(temp,:) = mean(data,1,'omitnan');
    yStdErr = std(data,1,'omitnan');%./sqrt(N); 
        %standard Error
    ci95 = tinv(0.025, N-1); 
        %95 confidence interval
    yci95 = yStdErr.*ci95;
    err(temp,:) = yStdErr;%*yci95;
%calculate 95% CI for temeprature
	errTemp = 0.5;
    xneg(temp,:) = Temps(obsTemps(temp))-(Temps(obsTemps(temp))+errTemp);
    xpos = xneg;
end

%iterate over every score to to plot means of every temp
    ylimRange = [(min(Temps)-1.5) (max(Temps)+1.5)];
    xlimRange = [-120 3750];
    %xlimRange = [-24 750];
    LineWidth = 2.5;
    MarkerSize = 50;
    CapSize = 0.003;
    
for score = 1:length(obsScores)
    hold on;
%plot basic for T0, specific with errors for T1 on
    if score == 1
        RGBa(obsScores(score),4) = 1;
        x = zeros(length(obsTemps),1)';
        y = Temps(obsTemps)';
        xError = zeros(1,length(x));
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x,y,RGBa(obsScores(score),:),xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
        %er = errorbar(zeros(length(obsTemps),1),Temps(obsTemps),xneg,...
         %   'vertical','o','Color',oColors(obsScores(score),:),'LineWidth',2.5);
        %plot(zeros(length(obsTemps),1),Temps(obsTemps),'-','Color',oColors(obsScores(score),:),'LineWidth',2.5);
    elseif sum(score == [1 8 12]) ~= 1
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
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
        %FILTER OUT SINGLE POINTS WITH NO ERROR CAUSE ONLY ONE EMBRYO?
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
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
        %idx = ~any(isnan(meanAccByStage(:,score)),2);
        %devTime = meanAccByStage(:,score);
        %errTime = err(:,score);
        %t = Temps(obsTemps);    
        %RGBa(obsScores(score),4) = 1;
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        %errorbar(devTime(idx),t(idx),xneg(idx),xpos(idx),errTime(idx),...
         %   errTime(idx),'-o','Color',oColors(obsScores(score),:),'LineWidth',2.5);
        %FILTER OUT SINGLE POINTS WITH NO ERROR CAUSE ONLY ONE EMBRYO?
    end
end
%Aestetics
    title('New + Old: Frog stages developmental times over range of temperatures')
    %title('New + Old Frog stages developmental times over range of temperatures')
    %title('New Frog stages developmental times over range of temperatures')

    ylim(ylimRange);
    ylabels = str2num(cell2mat(yticklabels));
    %yticks(ylabels([1 3 5 7 9 11])+0.15);
    yticklabels(cellstr(num2str((ylabels+0.15)-273.15)));
    ylabel('Temperature (^oC)');
    xlim(xlimRange)
    xlabel('Time post 3^{rd} Cleavage (min)');
    %displayNames = stageAbbs(obsScores);
    displayNames = stageAbbs(obsScores);
    %dummy figure to assign colors in legend
    clear sp
    for i = 1:length(obsScores)
        sp(i) = plot(1,1,'Color',RGBa(i,:),'LineWidth',LineWidth);
    end
    legend(sp,cellstr(displayNames),'NumColumns',2);
    set(gca,'LineWidth', 3, 'FontSize', 35);
%%
%Figure 1 D Log version
%as above, but using a personalized errorbar function to allow for
%transparency
%plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
%Raw Data Presentation

%generate as seperate figures

%array of temperatures to show (only being selective when there
%are temperatures very close together)
    obsTemps = [1:length(accumulated)];

%Figure 1, all data
    figureWidth = 1600;
    figureHeight = 675;
    figure('Position',[100 100 figureWidth figureHeight])
    aspectRatio = figureWidth/figureHeight;
%select what scores to look at, where 1 is T0
    obsScores = (1:12);
        %drop muscle twitching,13
    colors = distinguishable_colors(numScores);
    RGBa = horzcat(colors,zeros(numScores,1));
   
%initialize display data arrays for means and errors
    meanAccByStage = zeros(sum(length(obsTemps)),length(obsScores));
    err = zeros(sum(length(obsTemps)),length(obsScores));

%iterate over every observable temperature to condense multiple embryos in
%each score to a mean and calculate errors
    xneg = nan(length(obsTemps),1);
for temp = 1:length(obsTemps)
%parse out score data for specific temp
    data = log(accumulated{obsTemps(temp)}(:,obsScores)/60);
    data(:,1)=0;
%determine number of measurable samples for each score
    N = sum(~isnan(data),1);
%take mean of each score (no NaNs) then calculate 95% CI and store
    meanAccByStage(temp,:) = mean(data,1,'omitnan');
    yStdErr = std(data,1,'omitnan')./sqrt(N); 
        %standard Error
    ci95 = tinv(0.025, N-1); 
        %95 confidence interval
    yci95 = yStdErr.*ci95;
    err(temp,:) = yStdErr;%yci95;
%calculate 95% CI for temeprature
	errTemp = 0.5;%1.96*0.5;
    xneg(temp,:) = Temps(obsTemps(temp))-(Temps(obsTemps(temp))+errTemp);
    xpos = xneg;
end

%iterate over every score to to plot means of every temp
    ylimRange = [(min(Temps)-1.5) (max(Temps)+1.5)];
    xlimRange = [-log(120) log(3750)];%[-120 3750];
    %xlimRange = [-log(24) log(750)];%[-24 750];
    LineWidth = 2.5;
    MarkerSize = 50;
    CapSize = 0.003;
    
for score = 1:length(obsScores)
    hold on;
%plot basic for T0, specific with errors for T1 on
    if score == 1
        RGBa(obsScores(score),4) = 1;
        x = zeros(length(obsTemps),1)';
        y = Temps(obsTemps)';
        xError = zeros(1,length(x));
        yError = xneg';
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x,y,RGBa(obsScores(score),:),xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
        %er = errorbar(zeros(length(obsTemps),1),Temps(obsTemps),xneg,...
         %   'vertical','o','Color',oColors(obsScores(score),:),'LineWidth',2.5);
        %plot(zeros(length(obsTemps),1),Temps(obsTemps),'-','Color',oColors(obsScores(score),:),'LineWidth',2.5);
    elseif sum(score == [1 8 12]) ~= 1
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
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
        %FILTER OUT SINGLE POINTS WITH NO ERROR CAUSE ONLY ONE EMBRYO?
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
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        plotPersonal(x(idx),y(idx),RGBa(obsScores(score),:),xError(idx),yError(idx),LineWidth,MarkerSize,CapSize,ylimRange,xlimRange,aspectRatio)
        %idx = ~any(isnan(meanAccByStage(:,score)),2);
        %devTime = meanAccByStage(:,score);
        %errTime = err(:,score);
        %t = Temps(obsTemps);    
        %RGBa(obsScores(score),4) = 1;
        %plotPersonal(x,y,RGBa,xError,yError,LineWidth,MarkerSize,CapSize,ylimRange,xlimRange)
        %errorbar(devTime(idx),t(idx),xneg(idx),xpos(idx),errTime(idx),...
         %   errTime(idx),'-o','Color',oColors(obsScores(score),:),'LineWidth',2.5);
        %FILTER OUT SINGLE POINTS WITH NO ERROR CAUSE ONLY ONE EMBRYO?
    end
end
%Aestetics
    title('Frog stages developmental times over range of temperatures')
    ylim(ylimRange);
    ylabels = str2num(cell2mat(yticklabels));
    %yticks(ylabels([1 3 5 7 9 11])+0.15);
    yticklabels(cellstr(num2str((ylabels+0.15)-273.15)));
    ylabel('Temperature (^oC)');
    %xlim(xlimRange)
    xlabel('log_e(time) post 3^{rd} Cleavage [log_e[(min)]');
    %displayNames = stageAbbs(obsScores);
    displayNames = stageAbbs(obsScores);
    %dummy figure to assign colors in legend
    clear sp
    for i = 1:length(obsScores)
        sp(i) = plot(1,1,'Color',RGBa(i,:),'LineWidth',LineWidth);
    end
    legend(sp,cellstr(displayNames),'NumColumns',2);
    set(gca,'LineWidth', 3, 'FontSize', 35);
%%
%Raw (1/T) Exponential BIC
figure('Position',[100 100 650 1400])

%array of temperatures to show (only being selective when there
%are temperatures very close together)
    obsTemps = [1,2,3,5,8,9,10,12,13,16];
    colors = distinguishable_colors(numScores+1);

%Initialize Fit comparison Table
    AICTable = table('Size',[numScores,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
    AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};
    AICTable.Properties.RowNames = stageNames;
    
%initialize fit comparison arrays
    BICComps = nan(numScores,1);

%iterate over stages to calculate means for every temperature to fit
for stage = 1:(length(accumulated{1})-1)
    %drop muscle twitching 

    hold on;

%calc total num of fittable embryos and means per stage for all temps, 
%2 or more reps
	fittabeEmNum = 0;
	fittableMean = 0;
	for temp = obsTemps
        if (temp ~= 1)&&(temp ~= length(accumulated)-1)&&...
            (temp ~= length(accumulated))
            if sum(~isnan(accumulated{temp}(:,stage))) > 1
                fittabeEmNum = fittabeEmNum + ...
                sum(~isnan(accumulated{temp}(:,stage)));
            	fittableMean = fittableMean + 1;
            end
        end
	end

        
%initialize arrays to hold all fittable points and means X,Y        
	allX = nan(fittabeEmNum,1);
	allY = nan(fittabeEmNum,1);
	allMeans = nan(fittableMean,1);
	allXM = nan(fittableMean,1);
	emNum = 1;
	mNum = 1;  
        
%iterate over temperaters to average embryos for each temp then plot
    for temp = obsTemps
    %parse temperary 'data' of nans
    	data = accumulated{temp}(:,stage);
    	data = data(~isnan(data));
            
    %set/calculate errors
    	N = size(data,1); 
            %number non-nan experiments in data
    	yStdErr = std(data)/sqrt(N-1); 
            %standard Error
    	ci95 = tinv([0.025 0.975], N-1); 
            %95 confidence interval
    	yci95 = yStdErr.*ci95;
    	yneg = yci95(1);
    	ypos = yci95(2);
     	errTemp = 1.96*0.5; 
            %1.96 is multiplier associated with a two sided confidence interval
     	xneg = 1./Temps(temp)-1./(Temps(temp)+errTemp);
     	xpos = xneg;
            
    %build all data for specific temp for linear regression, while
    %plotting error bar
    	meanLnPseudoK = mean(data);
        if isnan(meanLnPseudoK)==1
    	elseif (temp ~= 1)&&(temp ~= length(accumulated)-1)&&...
        	(temp ~= length(accumulated))&&length(data)>1
        
        	nextEmNum = emNum+size(data,1);
        	allY(emNum:nextEmNum-1,1) = data;
        	allX(emNum:nextEmNum-1,1) = 1./repmat(Temps(temp),...
            	size(data,1),1);
            emNum = nextEmNum;
        %plot
            fitted = errorbar(1./Temps(temp),meanLnPseudoK,yneg,ypos,xneg,...
            	xpos,'o','Color',colors(stage,:),'LineWidth',2);
        %sort means and associate temp for fitting
            allMeans(mNum) = mean(data);
            allXM(mNum) = 1./Temps(temp);
            mNum = mNum + 1;
      	elseif length(data)==1
        %plot non-fitted data with different color, too few reps
        	errorbar(1./Temps(temp),meanLnPseudoK,yneg,ypos,xneg,xpos,...
                '*','Color',colors(numScores+1,:),'MarkerEdgeColor',...
                colors(stage,:),'MarkerFaceColor',colors(stage,:),'LineWidth',2);                
        else
        %plot non-fitted data with different color, non-Arrhenius range
        	errorbar(1./Temps(temp),meanLnPseudoK,yneg,ypos,xneg,xpos,...
            	'*','Color',colors(numScores+1,:),'MarkerEdgeColor',...
                colors(stage,:),'MarkerFaceColor',colors(stage,:),'LineWidth',2);
        end    
    end
    
%combine allX/Y
	allData = horzcat(allX,allY);
        
%store data points for AIC, each embryo as individual
    %{
	AICTable.numPoints(stage) = length(allData(:,2));
	AICTable.yi(stage) = num2cell(allData(:,2),1);
	AICTable.xi(stage) = num2cell(allData(:,1),1);

%plot linear regression
	f = fit(allData(:,1),allData(:,2), 'exp1');
	mdl = fitlm(allData(:,1), allData(:,2));

%AIC Comparison calculations, exp is assumed best fit (row stages, 
%columns are poly compared to) 
    if stage == 1
    else
    	polyNMdl2 = polyfit(allData(:,1), allData(:,2),1);
    	polyNMdl1 = fit(allData(:,1), allData(:,2),'exp1');
    	BICComps(stage,1) = BICExpCompCalc(AICTable,polyNMdl1,polyNMdl2,stage);
    end
    %}

%each temp as a replicate
	AICTable.numPoints(stage) = length(allMeans);
    AICTable.yi(stage) = num2cell(allMeans,1);
    AICTable.xi(stage) = num2cell(allXM,1);
        
%plot linear regression
	f = fit(allXM,allMeans, 'exp1');
	mdl = fitlm(allXM,allMeans);

%AIC Comparison calculations, exp is assumed best fit (row stages, 
%columns are poly compared to) 
    if stage == 1
    else
    	polyNMdl2 = polyfit(allXM,allMeans,1);
    	polyNMdl1 = fit(allXM,allMeans,'exp1');
    	BICComps(stage,1) = BICExpCompCalc(AICTable,polyNMdl1,polyNMdl2,stage);
    end
    
%plot fit lines
    if stage == 1
    	lm = plot([1./303.15 1./283.15],[0 0]);
    	lm.LineWidth = 2;
    	lm.Color = colors(stage,:);
    else
    	lm = plot(f);
    	lm.LineWidth = 2;
    	lm.Color = colors(stage,:);
    end
	ex(stage) = lm;
	legend off
end  

hold on;
    %beautify  
        ylabel('time (min)');
        xlabel('1/Temperature (C^{-1})');
            
        %Set consistant limits to compare
            xlim([1./303.15 1./283.15])
            ylim([-120 4000])

        %adjust xticks
            xticks(1./[303.15 293.15 283.15]);
            xticklabels({'1/30','1/20','1/10'});

        %set gca to be legible    
            set(gca,'linewidth', 2.5, 'FontSize', 28);
            title({'Frog stage developmental times';'post 3rd cleavage'})
            legend(ex,stageAbbs(1:12),'NumColumns',2,'Location','northwest')

    %AIC graph
        figure('Position',[100 100 850 800])

        bar(1:12,BICComps(1:12))
        xticklabels(cellstr(stageAbbs(1:12)))
        xtickangle(45)
        title({'BIC comparison between prefered'; 'exponential and linear fits'})
        ylabel({'BIC comparison value','= exp((BICexp???BIClinear)/2)'})
        set(gca,'linewidth', 2.5, 'FontSize', 28);

%%
%Figure 2E
%CV Analyis
%Data Import

%Destination
    filename = 'UnfilteredFrogScoresFinal.xlsx';
    folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Frog Scores/';

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
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+1;
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
    %data = readcell(strcat(filename),'Sheet', sheet_name{i});
    
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
figure('position',[100 100 1000 800])

%initialize matrix to hold presentation data
    CVMat = nan(numScores);
    
%iterate over ever (forward) stage interval to calculate co-variance
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
    end
end

%determine heatmap presentation area, slightly limited here
ylabels = stageAbbs(2:numScores);
xlabels = stageAbbs(1:numScores-1);
h = heatmap(xlabels, ylabels, round(CVMat(2:numScores,1:numScores-1),1),'CellLabelColor','black');

%Aestetics
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    %h.ColorbarVisible = 'off';

    caxis([0 20])

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
    title('\fontsize{38}Mean frog interval CV');
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
        cbh.TickLabels = [0:20/2:20]; 
        % Set the colorbar fontsize to the same value as heatmap fontsize
        cbh.FontSize = h.FontSize;
%%
%Selected in depth CV Analyis
%Data Import (rerun actual data after this because this script overides the
%good data loaded in before)

%
figure('Position',[10 10 800 600]);

meanPerStage = cell2mat(cellfun(@(x)mean(x,1,'omitnan'),perStage,'UniformOutput',false));
stdPerStage = cell2mat(cellfun(@(x)std(x,1,'omitnan'),perStage,'UniformOutput',false));
cvMat = (stdPerStage./meanPerStage)*100;
notEmpty = sum(isnan(cvMat),2)~=13;
cvMatNotEmpty = cvMat(notEmpty,:)
cvMatNoT0 = cvMatNotEmpty(:,[2:12]);
plot(Temps(notEmpty)-273.15,cvMatNoT0,'LineWidth',1.5)
title("Temporal Variance Over Temperatures in Frogs")
xlabel("Temperature ^oC")
ylabel("Coefficent of Variance (CV%)")


labs = strcat(stageAbbs(1:11)," - ", stageAbbs(2:12));
legend(labs,'NumColumns',4,'Location','northwest')
set(gca,'FontSize', 25,'LineWidth',2);
%%
%Figure S3D
%Arrhenius Stage-wise
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    
    embryoMass = 1.834; %mg

%Plot Stages
    figure('Position',[10 10 1350 850]);
    %colors = distinguishable_colors(3);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    obsScores = (2:12);
    perStageInterval = cellfun(@(x) x(:,obsScores),perStage,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 16 17 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];%[12 13]
        %for original Final data [1 16 17 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];
        %for combined data [1 11 14:16 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];
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
        pseudoKperStage{temp} = 1./perStageInterval{temp};
    end
    invertedTemps = 1./Temps;
    xmin = 1/305.15;
    xmax = 1/281.15;
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
    stageRou = zeros(length(lnPseudoKperStage{1}),1);
    stageQuadEa = zeros(length(lnPseudoKperStage{1}),1);
    stageQuadA = zeros(length(lnPseudoKperStage{1}),1);
    
%just 3, use -> justThree = (7,10,12);
    justThree = [4,11];
    justThreeCounter = 1;

%for all use -> length(lnPseudoKperStage{1})
    [numRow,numCol] = calcRectangle(length(lnPseudoKperStage{1}));
    sampleSize = NaN(1,length(lnPseudoKperStage{1}));
    overAllSampleSize = NaN(1,length(lnPseudoKperStage{1}));
    numRuns = 5000;
    bsEas = NaN(numRuns,length(lnPseudoKperStage{1}));
    bsAs = NaN(numRuns,length(lnPseudoKperStage{1}));
for stage = 1:length(lnPseudoKperStage{1})
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
            elseif temp ~= 12 && temp ~= 13
                overAllEms = overAllEms + sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
            end
        end
    
    %initialize and store subplot data
        sp(stage) = subplot(numRow,numCol,stage);   

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
                %standard Error
            ci95 = tinv([0.025 0.975], N-1); 
                %95 confidence interval
            yci95 = repmat(yStdErr,1,2);%yStdErr.*ci95;
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 0.5;%1.96*0.5; 
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
                allMeans(mNum) = mean(data);
                allXM(mNum) = invertedTemps(temp);
            %fit to individuals
                %fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
            %fit to means
                fitted = errorbar(allXM(mNum),allMeans(mNum),yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
                mNum = mNum + 1;
            elseif length(data)==1
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);                
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);
            end
            %copy over overAllData
            if temp ~= 12 && temp ~= 13
            nextOEmNum = oEmNum+size(data,1);
                overAllY(oEmNum:nextOEmNum-1,1) = data;
                overAllX(oEmNum:nextOEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                oEmNum = nextOEmNum;
            end
     
    end
    %combine allX/Y
        allData = horzcat(allX,allY);
        %record sample size
            sampleSize(stage) = length(allData);
            overAllSampleSize(stage) = length(overAllY);
    %store regression stats
        %individual points
        %
            stagRegT.numPoints(stage) = length(allData);
            stagRegT.SX(stage) = sum(allData(:,1));
            stagRegT.SY(stage) = sum(allData(:,2));
            stagRegT.SXX(stage) = sum(allData(:,1).^2);
            stagRegT.SYY(stage) = sum(allData(:,2).^2);
            stagRegT.SXY(stage) = sum(allData(:,1).*allData(:,2));
            stagRegT.AveX(stage) = mean(allData(:,1));
            stagRegT.AveY(stage) = mean(allData(:,2));
        %{
        %means
            stagRegT.numPoints(stage) = length(allMeans);
            stagRegT.SX(stage) = sum(allXM);
            stagRegT.SY(stage) = sum(allMeans);
            stagRegT.SXX(stage) = sum(allXM.^2);
            stagRegT.SYY(stage) = sum(allMeans.^2);
            stagRegT.SXY(stage) = sum(allXM.*allMeans);
            stagRegT.AveX(stage) = mean(allXM);
            stagRegT.AveY(stage) = mean(allMeans);
            %}
        if stage == 1
            frogRegT = stagRegT(stage,:);
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
        %lineFit = fit(allXM, allMeans, 'poly1');
            %lineFit2 = polyfit(allData(:,1), allData(:,2),2);
            %lineFit3 = polyfit(allData(:,1), allData(:,2),3);
            %[lineFit4,S4,mu4] = polyfit(allData(:,1), allData(:,2),8);
    
    %plot fit lines
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'b', 'LineWidth',2.5);
        text(1./(273.15+22.5),-6.25,strcat("E_a = ",num2str(round((-lineFit.p1*R)*tokJConvert,0))," kJ/mol"),'Color','b','FontSize',19)
        %yLin2 = polyval(lineFit2,invertedTemps);
        %lm2 = plot(invertedTemps,yLin2, 'm', 'LineWidth',2);
        %yLin3 = polyval(lineFit3,invertedTemps);
        %lm3 = plot(invertedTemps,yLin3, 'g', 'LineWidth',2);
        %yLin4 = polyval(lineFit4,invertedTemps,S4,mu4);
        %lm4 = plot(invertedTemps,yLin4, 'r', 'LineWidth',2);
        %calculate slope
        curvFit = polyfit(allData(:,1), allData(:,2),2);
        stageRou(stage) = -curvFit(1);
            %negative cause mdl fits to +nx2+mx+b
            %but our equation form is lnA - Ea/RT - p/T2
            %therefore b = lnA, m = -Ea, n = -p
            %n = -p means curveFit(1) is negative of p
        stageQuadEa(stage) = -(curvFit(2)*R)*tokJConvert;
        stageQuadA(stage) = curvFit(3);  
        curvOFit = polyfit(overAllX, overAllY,2);
        xra = min(xrange):0.00001:max(xrange);
        yOCurv = polyval(curvOFit,xra);
        ocm = plot(xra,yOCurv,'--r', 'LineWidth',2.5);
        xra = min(xrange):0.00001:max(xrange);
        yCurv = polyval(curvFit,xra);
        cm = plot(xra,yCurv,'--b', 'LineWidth',2.5);
        
        
    %Plot parameters
        %int and slope
            slope(stage) = lineFit.p1;
            intercept(stage) = lineFit.p2;
        %get confidence intervals
            ci = confint(lineFit,0.68);%confint(lineFit,0.95);
            stageActCi(:,stage) = ci(:,1);
            yIntCI(:,stage) = ci(:,2);
            
        CIs = confint(lineFit,0.95);    
        plusCI = lineFit.p1 - CIs(1,1);
        round((plusCI*R)*tokJConvert,0);
        
    %beautify
        %Set consistant limits to compare
            ylim([-13.4 -5.6]) 
                %for all plots, to compare
            yticks([-12 -10 -8 -6])
            xlim(xrange);

        %adjust xticks
            xticks(1./[299.1500 285.1500]);
            xticklabels({'1/26','1/12'});

        %set labels on outer most subplots, remove from interior        
            if stage==1||stage==5||stage==9
                ylabel('ln(1/time [s^{-1}])');
            else
                set(gca,'yticklabel',[])
            end
            if stage>=8 && stage<=12
                xlabel('1/(T [^oC] + 273.15)');
            else
                set(gca,'xticklabel',[])
            end

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 23);

        %place titles as text inside plots
            if stage == 1
                text(mean(xrange),-12.5,{strcat("A to ", stageAbbsNoT0(stage))},...
                    'FontSize', 24,'HorizontalAlignment','center');
                %{
            elseif stage>=2 && stage<=7
                lStage = split(stageAbbsNoT0(stage-1),' ');
                text(mean(xrange),-12.5,{strcat(lStage(1)," to ",...
                    (stageAbbsNoT0(stage)))}, 'FontSize', 24,...
                    'HorizontalAlignment','center');
                %}
            else
                text(mean(xrange),-12.5,{(strcat(stageAbbsNoT0(stage-1),...
                    " to ",(stageAbbsNoT0(stage))))}, 'FontSize', 24,...
                    'HorizontalAlignment','center');
            end

    justThreeCounter = justThreeCounter + 1;
end

%Adjust Subplot positions to tighten image
    %Shift left
        for stage = 1:length(obsScores)
            if stage ~= 1 && stage ~= 5 && stage ~= 9
                prevPos = get(sp(stage-1),'Position');
                prevPos = [prevPos(1)+prevPos(3)+0.0075 prevPos(2) ...
                    prevPos(3) prevPos(4)];
                set(sp(stage),'Position', prevPos)
            end
        end
    %shift up
        for stage = 1:length(obsScores)
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
%violin plots

%violin(bsEas([1,3:12]),[])
means = NaN(1,11);
%own version
figure('Position',[10 10 1500 800]);
for i = [1:11]
    if i == 1
       h = subplot(1,12,i);
       histogram(bsEas(:,i),50,'EdgeAlpha',0);
       %set(h,'Position',[0.1308 0.1100 0.7890 0.1100])
    else
       h = subplot(1,12,i-1);
       histogram(bsEas(:,i),50,'EdgeAlpha',0);
    end
    %calc and plot mean/median
        meanEas = mean(bsEas(:,i));
        means(i) = meanEas;
        medianEas = median(bsEas(:,i));
    %plot median and mean
        black = xline(meanEas,'k--','linewidth',2);
            %correct these so the dont "spill over" the hist
        red = xline(medianEas,'r--','linewidth',2);
    
    xlim([40 110])
    ylim([0 400])
    set(h,'view',[90 -90])
    axis off
end
axes('position',[0.1308-0.01 0.1100-.00155 0.7890-0.1308+0.0571 0.8150],'color','none') %[left bottom width height]
    %run once to determine proper axes size/positions, check h
    %subplot1pos = [0.1308 0.1100 0.0462 0.8150]
    %subplot2pos = [0.1959 0.1100 0.0471 0.8150]
    %subplot11pos = [0.7890 0.1100 0.0501 0.8150]
xticks([0.015:0.0922:1])
xlabs = strcat(stageAbbs(1:11)," - ", stageAbbs(2:12));
xticklabels(xlabs)
xlabel("Stage intervals")
yticks([0:1/7:1])
yticklabels([40:10:110])
ylabel("Apparent activation energy (kJ)")
legend([black red],["mean E_a","median E_a"],'Location','northwest','LineWidth',2.5,'FontSize',20)
set(gca,'lineWidth',2.5,'FontSize',20)
title(["Distribution of Bootstrapped Activation Energies";"Frog Developmental Stages"],'FontSize',25)

%%    
%Figure 2B
%Arrhenius Stage-wise, just 2
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;

%Plot Stages
    figure('Position',[10 10 1000 400]);
    %colors = distinguishable_colors(3);
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    obsScores = (1:numScores);
    perStageInterval = cellfun(@(x) x(:,obsScores),perStage,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [1 16 17 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];

%convert to pseudo Ks and inverse temps
    pseudoKperStage = cell(length(perStageInterval),1);
    for temp = 1:length(perStageInterval)
        pseudoKperStage{temp} = 1./abs(perStageInterval{temp});
    end
    invertedTemps = 1./Temps;
    xmin = 1/307.15;
    xmax = 1/280.15;
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
for stage = [8,12]
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps
        fittabeEmNum = 0;
        fittableMean = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0
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
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %fit to individuals
                fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
            %fit to means
                allMeans(mNum) = mean(data);
                allXM(mNum) = invertedTemps(temp);
                %fitted = errorbar(allXM(mNum),allMeans(mNum),yneg,...
                    %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
                mNum = mNum + 1;
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

    %plot linear regression
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');
        %lineFit = fit(allXM, allMeans, 'poly1');

    %plot fit lines
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'b', 'LineWidth',2.5);
        CIs = confint(lineFit,0.95);
        plusCI = lineFit.p1 - CIs(1,1);
        minusCI = lineFit.p1 - CIs(2,1);
        text(1./(273.15+20.5),-6.75,strcat("E_a = ",num2str(round((-lineFit.p1*R)*tokJConvert,0))," +/- ", num2str(round((plusCI*R)*tokJConvert,0))),'Color','b','FontSize',24)
        text(1./(273.15+15.5),-7.25,"kJ/mol",'Color','b','FontSize',24)

        
    %Plot parameters
        %int and slope
            slope(stage) = lineFit.p1;
            intercept(stage) = lineFit.p2;
        %get confidence intervals
            ci = confint(lineFit,0.95);
            stageActCi(:,stage) = ci(:,1);
            yIntCI(:,stage) = ci(:,2);
        
    %beautify
        %Set consistant limits to compare
            ylim([-12 -6]) 
                %for all plots, to compare
            yticks([-11 -9 -7])
            xlim(xrange);

        %adjust xticks
            xticks(1./[299.1500 285.1500]);
            xticklabels({'1/26','1/12'});

        %set labels on outer most subplots, remove from interior        
                ylabel('ln(1/time [s^{-1}])');
                xlabel('1/(T [^oC] + 273.15)');

        %set gca to be legible    
            set(gca,'linewidth', 3, 'FontSize', 28);

        %place titles as text inside plots
                %text(mean(xrange),-12.5,{strcat(stageAbbsNoT0(stage),...
                 %   " to ",stageAbbsNoT0(stage+1))},...
                  %  'FontSize', 22,'HorizontalAlignment','center');
                title(strcat(stageAbbsNoT0(stage-1),...
                    " to ",stageAbbsNoT0(stage)),...
                    'FontSize', 28);
         
    pannelCount = pannelCount + 1;
end
%%
%Figure 3C
%Arrhenius Accumulated
   figure('Position',[100 100 1250 450]);
    
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    
    embryoMass = 1.834; %mg
    
%remove T0 time point from analysis
    obsScores = (2:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedTemps = [16 17];%[1 16 17 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];
    
%Mass correction by t/m^(1/4)
    %accumulatedInterval = cellfun(@(x) x(:,obsScores)./embryoMass^0.25,accumulated,'UniformOutput',false);
    
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
        %AICTable.Properties.RowNames = stageNames;         
        
%convert to pseudo Ks and inverse time
    pseudoKperStage = cell(length(accumulatedInterval),1);
    for temp = 1:length(accumulatedInterval)
        pseudoKperStage{temp} = 1./accumulatedInterval{temp};
    end
    invertedTemps = 1./Temps;
    xmin = 1/303.15;
    xmax = 1/283.15;
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
    
    %just 3, use -> justThree = (7,10,12);
    justThree = [7,11];
    justThreeCounter = 1;

%for all use -> length(lnPseudoKperStage{1}) 
    [numRow,numCol] = calcRectangle(length(lnPseudoKperStage{1}));
for stage = justThree
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps only
        fittabeEmNum = 0;
        fittableMean = 0;
        for temp = 1:length(invertedTemps)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) >0%>= 2
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
        allMeans = nan(fittableMean,1);
        allXM = nan(fittableMean,1);
        mNum = 1;
        
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
            elseif sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0%>=2
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %fit to individuals
                %fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
            %fit to means
                allMeans(mNum) = mean(data);
                allXM(mNum) = invertedTemps(temp);
            %fit to means
                fitted = errorbar(allXM(mNum),allMeans(mNum),yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
                mNum = mNum + 1;            
            elseif length(data)==1
                %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                 %   xneg,xpos,'*','Color',colors(3,:),'LineWidth',2);                
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);
            end
    end
    
    %combine allX/Y
        allData = horzcat(allX,allY);
        
    %store regression stats
        %
        accRegT.numPoints(stage) = length(allData);
        accRegT.SX(stage) = sum(allData(:,1));
        accRegT.SY(stage) = sum(allData(:,2));
        accRegT.SXX(stage) = sum(allData(:,1).^2);
        accRegT.SYY(stage) = sum(allData(:,2).^2);
        accRegT.SXY(stage) = sum(allData(:,1).*allData(:,2));
        accRegT.AveX(stage) = mean(allData(:,1));
        accRegT.AveY(stage) = mean(allData(:,2));
        %
        AICTable.numPoints = length(allY);
        AICTable.yi = num2cell(allY,1);
        AICTable.xi = num2cell(allX,1);
        %}
        %{
	%means
        stagRegT.numPoints(stage) = length(allMeans);
        stagRegT.SX(stage) = sum(allXM);
        stagRegT.SY(stage) = sum(allMeans);
        stagRegT.SXX(stage) = sum(allXM.^2);
        stagRegT.SYY(stage) = sum(allMeans.^2);
        stagRegT.SXY(stage) = sum(allXM.*allMeans);
        stagRegT.AveX(stage) = mean(allXM);
        stagRegT.AveY(stage) = mean(allMeans); 
        
        
        AICTable.numPoints = length(allMeans);
        AICTable.yi = num2cell(allMeans,1);
        AICTable.xi = num2cell(allXM,1);
        %}
    %plot linear regression
        %lineFit = fit(allXM, allMeans, 'poly1');
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        mdl = polyfit(allData(:,1), allData(:,2), 1);
        %mdl = polyfit(allXM, allMeans, 1);
        lm = plot(xrange,yLin, 'b', 'LineWidth',2.5);
        %curvFit = polyfit(allData(:,1), allData(:,2),2);
        curvFit = polyfit(allXM, allMeans,2);
        xra = min(xrange):0.00001:max(xrange);
        yCurv = polyval(curvFit,xra);
        cm = plot(xra,yCurv,'--r', 'LineWidth',2.5);
    %calculate slope
        slope(stage) = lineFit.p1;
        intercept(stage) = lineFit.p2;
        
    %get confidence intervals
        ci = confint(lineFit,0.95);
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
            ylabs = cell2mat(cellfun(@str2num,yticklabels,'un',0));
            yticks(ylabs(ylabs == floor(ylabs)));
            yticklabels(ylabs(ylabs == floor(ylabs)));
            %text(1./(273.15+16.5),theoryMax-0.25,strcat("E_a = ",num2str(round((-slope(stage)*R)*tokJConvert,0))),'Color','b','FontSize',18)
            %text(1./(273.15+18.75),theoryMax-0.25,strcat('$ln(\frac{L_{Q}}{L_{L}}) = ',num2str(round(log(BICcomp),1)),'$'),'Color','m','FontSize',30,'Interpreter','latex')
            text(1./(273.15+18.75),theoryMax-0.25,strcat('$ln(\frac{L_{Q}}{L_{L}}) = ',num2str(round(log(BICcomp),1)),'$'),'Color','k','FontSize',30,'Interpreter','latex')

        %xlabels
            xlabel('1/(T [^oC] + 273.15)');
            xlim(xrange);
            xticks(1./[303.15 293.15 283.15]);
            xticklabels({'1/30','1/20','1/10'});
        %gca, titles, legend
            set(gca,'linewidth', 3, 'FontSize', 28);
            title(strcat("A to ",stageAbbsNoT0(stage)),'FontSize',28);
            legend([lm cm], {"Linear Fit","Quadratic"},'Location','southwest','FontSize',25);
    
    justThreeCounter = justThreeCounter + 1;
end

%Ea calculation
    accActEa = -(slope)*R*tokJConvert;
    accActA = exp(intercept);
    accActCi = -(accActCi)*R*tokJConvert;
%%
%Accumulated with cross validation
   %figure('Position',[100 100 1200 1100]);
   figure('Position',[100 100 1300 450]);
    
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    
%remove T0 time point from analysis
    obsScores = (2:numScores);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    
%initialize tables
    %stats table
        accRegM = zeros(length(obsScores),8);
        accRegT = array2table(accRegM);
        accRegT.Properties.VariableNames = {'numPoints','SX','SY','SXX','SYY','SXY','AveX','AveY'};
        accRegT.Properties.RowNames = stageNamesNoT0;
        
%convert to pseudo Ks and inverse time
    pseudoKperStage = cell(length(accumulatedInterval),1);
    for temp = 1:length(accumulatedInterval)
        pseudoKperStage{temp} = 1./accumulatedInterval{temp};
    end
    invertedTemps = 1./Temps;
    xmin = 1/307.15;
    xmax = 1/280.15;
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
    
    %just 3, use -> justThree = (7,10,12);
    justThree = [7,11];
    justThreeCounter = 1;

%for all use -> length(lnPseudoKperStage{1}) 
    [numRow,numCol] = calcRectangle(length(lnPseudoKperStage{1}));
for stage = justThree
    %calc total num of fittable embryos per stage for all temps, 2 or more
    %reps only
        fittabeEmNum = 0;
        fittableMean = 0;
        overAllEms = 0;
        for temp = 1:length(lnPseudoKperStage)
            if (temp ~= 1)&&(temp ~= length(lnPseudoKperStage)-1)&&...
                    (temp ~= length(lnPseudoKperStage))
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
                        fittableMean = fittableMean + 1;                    
                end
            end
            overAllEms = overAllEms + sum(~isnan(lnPseudoKperStage{temp}(:,stage)));
        end
        
    %initialize arrays to hold all data X and Y, and for means
        allX = nan(fittabeEmNum,1);
        allY = nan(fittabeEmNum,1);
        emNum = 1;
        allMeans = nan(fittableMean,1);
        allXM = nan(fittableMean,1);
        mNum = 1;
        overAllX = nan(overAllEms,1);
        overAllY = nan(overAllEms,1);
        oEmNum = 1;
        
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
            elseif (temp ~= 1)&&(temp ~= length(lnPseudoKperStage)-1)&&...
                    (temp ~= length(lnPseudoKperStage))&&length(data)>1
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %fit to individuals
                fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
            %fit to means
                allMeans(mNum) = mean(data);
                allXM(mNum) = invertedTemps(temp);
            %fit to means
                %fitted = errorbar(allXM(mNum),allMeans(mNum),yneg,...
                    %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
                mNum = mNum + 1;            
            elseif length(data)==1
                %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    %xneg,xpos,'*','Color',colors(3,:),'LineWidth',2);                
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'*','Color',colors(2,:),'LineWidth',2);
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
        %lineFit = fit(allXM, allMeans, 'poly1');
        lineFit = fit(allData(:,1), allData(:,2), 'poly1');
        yLin = polyval([lineFit.p1 lineFit.p2],xrange);
        lm = plot(xrange,yLin, 'b', 'LineWidth',2);
        if stage == 11
            curvFit = polyfit(allData(:,1), allData(:,2),2);
            xra = min(xrange):0.00001:max(xrange);
            yCurv = polyval(curvFit,xra);
            cm = plot(xra,yCurv,'--b', 'LineWidth',2);
        elseif stage == 7
            curvOFit = polyfit(overAllX, overAllY,2);
            xra = min(xrange):0.00001:max(xrange);
            yOCurv = polyval(curvOFit,xra);
            ocm = plot(xra,yOCurv,'--m', 'LineWidth',2);
        end
    %calculate slope
        slope(stage) = lineFit.p1;
        intercept(stage) = lineFit.p2;
        
    %get confidence intervals
        ci = confint(lineFit,0.95);
        accActCi(:,stage) = ci(:,1);
            
        hold on;
    %beautify
        %ylabels
            ylabel('ln(1/time)');
            meanY = round(real(mean(allY,'omitnan')),0);
            theoryMax = polyval([lineFit.p1 lineFit.p2],xrange(1));
            theoryMin = polyval([lineFit.p1 lineFit.p2],xrange(2));
            ylim([theoryMin-0.5,theoryMax]);
            ylabs = yticklabels;
            ylabs = cellfun(@str2num,ylabs);
            ylabs = ylabs(ylabs == floor(ylabs));
            yticks(ylabs)
            yticklabels(ylabs)
            %text(1./(273.15+16.5),theoryMax-0.5,strcat("E_a = ",num2str(round((-slope(stage)*R)*tokJConvert,0))),'Color','b','FontSize',18)

        %xlabels
            xlabel('1/T (C^{-1})');
            xlim(xrange);
            xticks(1./[303.15 293.15 283.15]);
            xticklabels({'1/30','1/20','1/10'});
        %gca, titles, legend
            set(gca,'linewidth', 2.5, 'FontSize', 28);
            title(lower(strcat("3rd Cleave. to ",stageAbbsNoT0(stage))),'FontSize',28);
            if stage == 7
                legend([lm ocm], {'Linear Fit','Quadratic: All'},'Location','southwest');
            elseif stage == 11
                legend([lm cm], {'Linear Fit','Quadratic: Lin'},'Location','southwest');
            end
%{        
    %Cross Val
    %initialize data vectors
    %linRange
    %xData = allX;
    %yData = allY;
    %all Range
    xData = overAllX;
    yData = overAllY;
    numRuns = length(5:length(xData)-4);
    meanTestMSEPerRun = zeros(4,numRuns);
    meanTrainMSEPerRun = zeros(4,numRuns);
    %for 2 to dataSize/2
    for testSize = 5:length(xData)-4
        %initialize parfor vectors
        parForTestMSE = nan(4,1000);
        parForTrainMSE = nan(4,1000);
        degrees = [1 2 3 4];
        %1000 parfor runs
        parfor i = 1:1000
            %Randomly split Data into two even parts (2 up to half thhe data)
            [tempXTrainData,idx] = datasample(xData,testSize,'Replace',false);
            tempYTrainData = yData(idx');
            %train data
            %tempXTrainData = tempXData(1:testSize);
            %tempYTrainData = tempYData(1:testSize);
            %test data
            tempXTestData = xData(setdiff(1:length(xData),idx));%tempXData(testSize+1:length(tempXData));
            tempYTestData = yData(setdiff(1:length(xData),idx));%tempYData(testSize+1:length(tempYData));
            mdlsTrainMSE = zeros(4,1);
            mdlsTestMSE = zeros(4,1); 
            if length(unique(tempXTrainData)) >= 5
            for degree = degrees
                %Calculate the various fits on training see (1,2,3,4,8) and RSS
                [mdl,S,mu] = polyfit(tempXTrainData,tempYTrainData,degree);
                esti = polyval(mdl,tempXTrainData,S,mu);
                RSS = sum((tempYTrainData - esti).^2);
                mdlsTrainMSE(degree) = RSS/length(tempXTrainData);
                %Use other data to test the validity of the fits
                esti = polyval(mdl,tempXTestData,S,mu);
                RSS = sum((tempYTestData - esti).^2);
                mdlsTestMSE(degree) = RSS/length(tempXTestData);
                    %Store the RSS for each data fit
            end
            parForTrainMSE(:, i) = mdlsTrainMSE;
            parForTestMSE(:, i) = mdlsTestMSE;
            end
        end
        %Store the mean for total 1000 run
       meanTestMSEPerRun(:,testSize-4) = median(parForTestMSE,2,'omitnan');
       meanTrainMSEPerRun(:,testSize-4) = median(parForTrainMSE,2,'omitnan');
    end
    %Report the mean RSS for each sample size
    colors = distinguishable_colors(4);
    plotHandles = nan(1,4);
    subplot(2,2,justThreeCounter+2)
    for i = [1 2 3 4]
        hold on
        scatter(5:(length(xData)-4),meanTestMSEPerRun(i,:),[],colors(i,:),'LineWidth',2)
        plotHandles(i) = plot(5:(length(xData)-4),meanTestMSEPerRun(i,:),'Color',colors(i,:),'LineWidth',2);
        %if (i/2) ~= floor(i/2)
        scatter(5:(length(xData)-4),meanTrainMSEPerRun(i,:),[],colors(i,:),'*','LineWidth',2)
        %end
        plot(5:(length(xData)-4),meanTrainMSEPerRun(i,:),'--','Color',colors(i,:),'LineWidth',2)
    end
    ylim([-0.00125 0.02])
    xlim([5 length(xData)-4])
    ylabel('MSE (RSS/n)')
    xlabel('Training Size')
    set(gca,'linewidth', 2.5, 'FontSize', 28);
    title({strcat("Cross Validation: ", lower(strcat("3rd Cleav to ",stageAbbsNoT0(stage))));'Linear Range Data'},'FontSize',18)
    title({strcat("Cross Validation: ", lower(strcat("3rd Cleav to ",stageAbbsNoT0(stage))));'All Data'},'FontSize',18)
    %title(num2str(i))
    %end
    legend(plotHandles([1 2 3 4]),{'1','2','3','4'},'FontSize',18);
%}
    justThreeCounter = justThreeCounter + 1;
end

%%
%Statistical Significance calculations for stage-wise
    %run the stages section first
%figure('Position',[100 100 1550 900]);
fh = figure('Position',[100 100 1100 900]);
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
    pValM = nan(length(obsScores),length(obsScores));
    pwrVals = nan(length(obsScores),length(obsScores));
    
    
for x = 1:length(obsScores)
    for y = x:length(obsScores)
        if x == y
            pValM(y,x) = 1;
            pwrVals(y,x) = 0;
                %stages compared against same stage should = 1
        else
            compStages = {char(stageNamesNoT0(x)),char(stageNamesNoT0(y))};
            compRegs = stagRegT(compStages,:);
        %run stage regression table through home-made ancova function
            pValM(y,x) = ancova(compRegs);
            pwrVals(y,x) = calcPower(compRegs);
        end
    end   
end

%plot 4 plots, one of p-values, then color coded for thhreshold sig at
%0.05,0.01,0.001
%pvals
    %subplot(2,2,1)
    h1 = heatmap(strcat(stageAbbs(1:11),'-',stageAbbsNoT0),strcat(stageAbbs(1:11),...
        '-',stageAbbsNoT0),round(pValM,4),'CellLabelColor','black','CellLabelFormat','%.1e');
    h1.Title = {'\fontsize{38}Frog ANCOVA p-values', '\fontsize{38}comparing individual stages'};
    h1.FontSize = 15;
    h1.MissingDataColor = [1 1 1];
    h1.GridVisible = 'off';
    caxis([0 1])
    noSig = horzcat(repmat(0.6,950,1),repmat(0.6,950,1),ones(950,1));
    singSig = horzcat(repmat(0.5,40,1),repmat(0,40,1),repmat(0.5,40,1));
    doubSig = horzcat(ones(9,1),repmat(0.3,9,1),ones(9,1));
    tripSig = horzcat(ones(2,1),repmat(0.3,2,1),repmat(0.1,2,1));
    map = vertcat(tripSig,doubSig,singSig,noSig);
    colormap(h1,map)
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Stage 1}')
    ylabel({'\fontsize{38}\bf{Stage 2}'})
    
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
        cbh.TickLabels = [0.001,0.01,0.05]; 
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
        
        
        
        
        
figure('Position',[100 100 1100 900]);        
%pwrValues
    %subplot(2,2,1)
    h1 = heatmap(strcat(stageAbbs(1:11),'-',stageAbbsNoT0),strcat(stageAbbs(1:11),...
        '-',stageAbbsNoT0),round(pwrVals,4),'CellLabelColor','black');
    h1.Title = {'\fontsize{38}Frog ANCOVA power analysis', '\fontsize{38}comparing individual stages'};
    h1.FontSize = 15;
    h1.MissingDataColor = [1 1 1];
    h1.GridVisible = 'off';
    caxis([0 1])
    noSig = horzcat(repmat(0.6,799,1),repmat(0.6,799,1),ones(799,1));
    singSig = horzcat(repmat(0.5,200,1),repmat(0,200,1),repmat(0.5,200,1));
    map = vertcat(noSig,singSig);
    colormap(h1,map)
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Stage 1}')
    ylabel({'\fontsize{38}\bf{Stage 2}'})
    
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

%{
%0.05
    subplot(2,2,2)
    signifDif = double(pValM<0.05);
    map = [0.25 .75 1
        1 .25 .5];
    h2 = heatmap(stageAbbsNoT0,stageAbbsNoT0,signifDif,'Colormap',map);
    h2.Title = 'Statistically Different Slopes, p>0.05';
    h2.FontSize = 15;
    
%0.01
    subplot(2,2,3)
    signifDif = double(pValM<0.01);
    map = [0.25 .75 1
        1 .25 .5];
    h2 = heatmap(stageAbbsNoT0,stageAbbsNoT0,signifDif,'Colormap',map);
    h2.Title = 'Statistically Different Slopes, p>0.01';
    h2.FontSize = 15;

%0.001
    subplot(2,2,4)
    signifDif = double(pValM<0.001);
    map = [0.25 .75 1
        1 .25 .5];
    h2 = heatmap(stageAbbsNoT0,stageAbbsNoT0,signifDif,'Colormap',map);
    h2.Title = 'Statistically Different Slopes, p>0.001';
    h2.FontSize = 15;
    set(gca,'FontSize', 18);
    %}
%%
%Figure 2D
%Comparing Eas with statistical significance between stage-wise
%Must run the above stages code chunk first
figure('Position',[100 100 900 700]);

%parse Ea data and CIs to present
    spos = stageActCi(1,:) - stageActEa;
    sneg = stageActEa - stageActCi(2,:);
    nonScaled = 1:length(obsScores);
    %scale by 21.1 C data
    tempScale = mean(accumulated{18},1,'omitnan');
 
%plot errorbar of data
    %{
    %scaled
        stageEaPlot = errorbar(tempScale(nonScaled), stageActEa(nonScaled), sneg(nonScaled), spos(nonScaled),'o','LineWidth',2.5);
    %}
        %nonsclaed
    %
    stageEaPlot = errorbar(nonScaled, stageActEa, sneg, spos,'o','LineWidth',2.5);
    %}
%{
%Scaled Aestetics    
    %developmentally scaled ticks
    xlim([min(tempScale)-2*abs(min(tempScale(2))) max(tempScale)+2*abs(min(tempScale(2)))])
    xticks([tempScale(nonScaled)])
    %tick aestethics
    stageEaPlot.Color = 'b';
    xticklabels(strcat(stageAbbs(1:12-1),"-",stageAbbs(2:12)))
    xtickangle(45)
    %xtickangle(45)
    set(gca,'linewidth', 3, 'FontSize', 28)
    title([{"Apparent activation energy", " Frog developmental stages"}], "FontSize", 28)
    %title(["Apparent E_a Fly: New Data"], "FontSize", 28)
    %title(["Apparent E_a Fly: Old Data"], "FontSize", 28)
    %title(["Apparent E_a Fly: Combined Data"], "FontSize", 28)

    ylabel("Activation energy (kJ/mol)", "FontSize", 28)
    xlabel({"Stage", "(Scaled by average developmental times between stages: 18.7 ^oC)"}, "FontSize", 28)   
%} 
%Non-scaled Aestetics    
    stageEaPlot.Color = 'b';
    xlim([0 12])
    xticks((1:11))
    xticklabels(stageAbbsNoT0)
    xticklabels(strcat(stageAbbs(1:12-1),"-",stageAbbs(2:12)))
    xtickangle(45)
    %xtickangle(45)
    set(gca,'linewidth', 3, 'FontSize', 28)
    title([{"Apparent activation energy", " Frog developmental stages"}], "FontSize", 28)
    %title(["Frog apparent activation: New Data"], "FontSize", 28)
    %title(["Frog apparent activation: Combined Data"], "FontSize", 28)
    ylabel("Activation energy (kJ/mol)", "FontSize", 28)
    xlabel("Stage", "FontSize", 28)
    ylim([30 130]);

%hard code brackets/braces based on the stats above
row = [1];
col = [7];
sigDifs = length([row,col]);
for i = 1:1
    hold on;
    rowCol = [row(i),col(i)];
%use Eas+errorbar to determine bracket start point
    compare = [stageActEa(row(i))-spos(row(i)) stageActEa(col(i))-spos(col(i))];
    largerEa = min(compare);
    smallerEa = max(compare);
%calculate row/col number associated with larger/smaller Ea
    largerTime = rowCol(compare==largerEa);
    smallerTime = rowCol(compare==smallerEa);
%line from larger Ea up to horizontal, allow 1 unit space between error
%bar and line
    line(repmat(nonScaled(largerTime),2,1), [largerEa-1 largerEa-2],'LineWidth',2.5,'Color','m');
%line from smaller Ea up to horizontal,allow 1 unit space between error
%bar and line
    line(repmat(nonScaled(smallerTime),2,1), [smallerEa-1 largerEa-2],'LineWidth',2.5,'Color','m');
%Horizontal line connecting above two lines
    line([nonScaled(largerTime) nonScaled(smallerTime)], repmat(largerEa-2,2,1),'LineWidth',2.5,'Color','m');
    text(mean([nonScaled(largerTime) nonScaled(smallerTime)]),smallerEa-6,'\fontsize{20}#','Color','m');

end

%Draw Braces between significantly different Eas, Use ancova calc above
%DRAWBRACE([X1,Y1], [X2,Y2], W)
%
%rightside up
    y1 = stageActEa(10)+spos(10)+1;
    drawbrace([4,y1], [10,y1], 10,'LineWidth',2.5,'Color','k');
        text(6.77,y1+4,'\fontsize{25}***');
    %y1 = stageActEa(9)+spos(9)+1;
    %testB = drawbrace([4,y1], [9,y1], 10,'LineWidth',2,'Color','k');
        %text(6.345,y1+4,'\fontsize{25}**');
%}
%upside down
    %y1 = stageActEa(8)-spos(8)-1;
    %drawbrace([9,y1], [8,y1], 10,'LineWidth',2,'Color','k');
        %text(8.27,y1-5.5,'\fontsize{25}**');
    %y1 = stageActEa(10)-spos(10)-1;
    %drawbrace([11,y1], [10,y1], 10,'LineWidth',2,'Color','m');
        %text(10.39,y1-5,'\fontsize{20}#','Color','m');

%%
%Figure 3D & S4F
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
        %AICTable.Properties.RowNames = stageNames;

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
    isConcaveBioFrog = nan(numScores);
    isConvexBioFrog = nan(numScores);
    isLinearBioFrog = nan(numScores);
    
%control what temperatures to exclude from the linear fit BIC comparison
%calculation
    counter = 1;

%calculate poly n 1 and 2 fits for every developmental interval then compare 
%them using BIC
for p = 1:2
    if p == 1
        excludedTemps = [1 11 14:16 length(invertedTemps)-4 length(invertedTemps)-3 length(invertedTemps)-2 length(invertedTemps)-1 length(invertedTemps)];%[12 13]
        %for original Final data [1 16 17 length(invertedTemps)-4 length(invertedTemps)-3 length(invertedTemps)-2 length(invertedTemps)-1 length(invertedTemps)];
        %for combined data [1 11 14:16 length(invertedTemps)-4 length(invertedTemps)-3 length(invertedTemps)-2 length(invertedTemps)-1 length(invertedTemps)];
    else
        excludedTemps = [61];
        %for original Final data [16 17];
        %for combined data [11 14:16];
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
        fittableMean = 0;
        for temp = 1:length(invertedTemps)
        %if not one of the excluded temperature and more then one embryo
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp})) > 0
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}));
                    fittableMean = fittableMean + 1;
                end
            end
        end
        
    %initialize arrays to hold all fittable points X,Y, and means        
    	allX = nan(fittabeEmNum,1);
    	allY = nan(fittabeEmNum,1);
    	allMeans = nan(fittableMean,1);
    	allXM = nan(fittableMean,1);
        emNum = 1;
     	mNum = 1;
        
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
            %means as data points
                allMeans(mNum) = mean(data);
                allXM(mNum) = invertedTemps(temp);
                mNum = mNum + 1;
                %fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
            elseif length(data)==1
                %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    %xneg,xpos,'*','Color',colors(3,:),'LineWidth',2);                
            else
                %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    %xneg,xpos,'*','Color',colors(2,:),'LineWidth',2);
            end
        end    
        %copy over errors for monte carlo
        if p ==1
            %linear
            monteTimeErrorsLinear(endStage,startStage) = mean(tempararyError,'omitnan'); 
            intervalAllXLinear{endStage,startStage} = 1./allX;
        else
            monteTimeErrorsExtremes(endStage,startStage) = mean(tempararyError); 
            intervalAllXExtremes{endStage,startStage} = 1./allX;
        end
        %Build BIC comparision variables, comparing models, then run
        %through BIC program
        if startStage~=endStage
        %store data points for AIC
        %all points as replicates
        	%
        	AICTable.numPoints = length(allY);
            AICTable.yi = num2cell(allY,1);
            AICTable.xi = num2cell(allX,1);
         	%{
        %each temp as a replicate
         	AICTable.numPoints = length(allMeans);
         	AICTable.yi = num2cell(allMeans,1);
          	AICTable.xi = num2cell(allXM,1);
            %}
        %AIC Comparison calculations, n=1 is assumed best fit (row stages, 
    	%columns are poly compared to) 
          	stage = 1;
                
          	%points as replicates
          	%
            	polyNMdl1 = polyfit(allX, allY, 1);%start = 12, end = 14, temp = 2
                	%problem is in scores that arent actually considered
             	polyNMdl2 = polyfit(allX, allY, 2);%there is some complex warning here...
           	%{
          	%means as replicates
            	polyNMdl1 = polyfit(allXM, allMeans, 1);%start = 12, end = 14, temp = 2
                	%problem is in scores that arent actually considered
               	polyNMdl2 = polyfit(allXM, allMeans, 2);
            %}
           	BICComps = fitBICCompCalc(AICTable,polyNMdl1,polyNMdl2,stage);
           	%AICComps = AICCompCalc(AICTable,polyNMdl1,S1,mu1,polyNMdl2,S2,mu2,stage);

            %Print out a specific interval to check fit comparison
            %{
                if startStage == 9 && endStage == 11
                %if (startStage == 11 && endStage == 11) || (startStage == 4

                    hold on
                    %scatter(allX, allY,'LineWidth',2);
                    scatter(allXM, allMeans,'LineWidth',2);
                    linVal = polyval([polyNMdl1.p1 polyNMdl1.p2],1./Temps);
                    %linSSR = sum((allY-(polyNMdl1.p1.*allX+polyNMdl1.p2)).^2);
                    linSSR = sum((allMeans-(polyNMdl1.p1.*allXM+polyNMdl1.p2)).^2);
                    lPlot = plot(1./Temps,linVal,'LineWidth',2);
                    curVal = polyval([polyNMdl2.p1 polyNMdl2.p2 polyNMdl2.p3],1./Temps);
                    %curSSR = sum((allY-(polyNMdl2.p1.*allX.^2+polyNMdl2.p2.*allX+polyNMdl2.p3)).^2);
                    curSSR = sum((allMeans-(polyNMdl2.p1.*allXM.^2+polyNMdl2.p2.*allXM+polyNMdl2.p3)).^2);
                    cPlot = plot(1./Temps,curVal,'LineWidth',2);
                    title({strcat(stageAbbsNoT0(startStage)," to ",stageAbbsNoT0(endStage))});
                    legend([lPlot cPlot],{strcat("linSSR = ",num2str(linSSR)),strcat("curSSR = ",num2str(curSSR))})
                    set(gca,'FontSize', 32,'LineWidth',2);
                    hold off
                end
            %}    
        %Store BIC comp values and Eas (just in case)
         preDec = split(num2str(abs(log(BICComps))),'.');
        if length(preDec{1}) >= 2
            BICMat(endStage,startStage) = round(log(BICComps),0);
        else
            BICMat(endStage,startStage) = round(log(BICComps),1);
        end
            %BICMat(endStage,startStage) = BICComps;
            EaMat(endStage,startStage) = (polyNMdl1(1)*-8.3144598)/1000; 
                %adjust slope to Ea with convertion Ea = -slope*R(in kJ)
        if p == 2 %2 extreme, 1 linear
        %Check concavity        
            isConcaveBioFrog(endStage,startStage) = polyNMdl2(1) > 0; 
            isConvexBioFrog(endStage,startStage) = polyNMdl2(1) < 0; 
            isLinearBioFrog(endStage,startStage) = polyNMdl2(1) == 0; 
        end
        else
        end
    end
end

%plot the relevent BIC data range and the associated stage abbreviations
    fh = figure('position',[100 100 1000 800])
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    %exp space
    %h = heatmap(xlabels, ylabels, round(BICMat(2:12,1:11),1),'CellLabelColor','none');
    h = heatmap(xlabels, ylabels, real(BICMat(2:12,1:11)),'CellLabelColor','black','CellLabelFormat','%.1f');
    %log space
    %h = heatmap(xlabels, ylabels, round(log(BICMat(2:12,1:11)),1));

%Aestetics
    %caxis([0 4])
    caxis([-5 5])
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    topx = flip(0.0:0.001:1);
    %map = horzcat((1*flip(sqrt(x))).', zeros(length(x),1),ones(length(x),1));
    %map = horzcat((1*flip(x)).', zeros(length(x),1),ones(length(x),1));
    topmap = horzcat(ones(length(topx),1),(1*(topx)).',(1*(topx)).');
    botx = flip(0.0:0.001:1);
    botmap = horzcat((1*flip(botx)).',(1*flip(botx)).',ones(length(botx),1));
    map = vertcat(botmap,topmap);
    colormap(h,map)
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Start score code}')
    ylabel({'\fontsize{38}\bf{End score code}'})
    
    %determine temperatures used and adjust title accordingly
        excul = setxor(1:16, excludedTemps);
        annotation('textbox',[0.950, 0.943, 0.00005, 0.00005], 'string', '+','FontSize',32,'LineWidth',0.01)
        annotation('textbox',[0.964, 0.1435, 0.00005, 0.00005], 'string', '-','FontSize',32,'LineWidth',0.01)
            %reduce textbox size to near zero to get rid of it
    if p == 1
    	%title({'\fontsize{38}BIC comparision polynomial 1 vs 2';'\fontsize{38}Temps = Linear Range'});
        %title({'\fontsize{38}BIC natural log likelihood ratio, frog';'\fontsize{38} quadratic vs liniar, linear range'});
        title('\fontsize{38}Frog ln(L_{Q}/L_{L}): Linear Range');
    else
    	%title({'\fontsize{38}BIC comparision polynomial 1 vs 2';'\fontsize{38}Temps = Extreme Range'});
        %title({'\fontsize{38}BIC natural log likelihood ratio, frog';'\fontsize{38} quadratic vs liniar, entire range'});
        title('\fontsize{38}Old: Frog ln(L_{Q}/L_{L})');
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
        %{
%As above, but for Eas
    figure('position',[100 100 1000 800])
    ylabels = stageAbbs(2:12);
    xlabels = stageAbbs(1:11);
    h = heatmap(xlabels, ylabels, round(EaMat(2:12,1:11),1));
   
%Aestetics
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    x = (0:0.001:1);
    map = horzcat((1*flip(x.^2)).',(1*x.^2).', zeros(length(x),1));
    colormap(h,map)
    set(gca,'FontSize', 32);
    title('\fontsize{38}Apparent Eas for same Temps');
    xlabel("Start Event")
    ylabel("End Event")
    
    figure('position',[100 100 1000 800])
    XY = round(BICMat(2:12,1:11));
    %clims = [0 2];
    im = imagesc(XY)%,clims)
    im.CData
    colorbar
    %colorbar('Ticks',[0,1,2],...
     %   'TickLabels',{'0','1','2+'})
    %}
%%
%Rank (1-3) ordered BIC for 1-8 degree polyies
%Colors
    %9, 8 for color coding, 1 for text
    colors = [0 0 1;...
        1 0 0;...
        0 1 0;...
        0 0 0.1724;...
        1 0.1034 0.7241;...
        1 0.8276 0;...
        0 0.3448 0;...
        0.5172 0.5172 1;...
        0.6207 0.3103 0.2759];

%Heatmap for AIC calculations
    
    obsScores = (1:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);    
    
%Initialize fit comparison Table
	AICTable = table('Size',[1,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
	AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};
        %AICTable.Properties.RowNames = stageNames;

%initialize fit comparison arrays
	AICValues = zeros(numScores,4);
        
%initialize Data for Arrhenius space x and y
    lnPseudoKperStage = cell(length(accumulatedInterval),1);
    invertedTemps = 1./Temps;

%initialize presentation data matrecies    
    rank1 = nan(numScores);
    rank2 = nan(numScores);
    rank3 = nan(numScores);
    
%control what temperatures to exclude from the linear fit BIC comparison
%calculation
    excludedTemps = [1 length(lnPseudoKperStage)-1 length(lnPseudoKperStage)];
    counter = 1;

%calculate poly n 1-8 fits for every developmental interval then rank 
%them using BIC values0
for startStage = 1:length(obsScores) 
    for endStage = startStage:length(obsScores) 
        BICs = nan(1,8);
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
        fittableMean = 0;
        for temp = 1:length(lnPseudoKperStage)
        %if not one of the excluded temperature and more then one embryo
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp})) > 1
                    fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}));
                    fittableMean = fittableMean + 1;
                end
            end
        end
        
    %initialize arrays to hold all fittable points X,Y, and means        
    	allX = nan(fittabeEmNum,1);
    	allY = nan(fittabeEmNum,1);
    	allMeans = nan(fittableMean,1);
    	allXM = nan(fittableMean,1);
        emNum = 1;
     	mNum = 1;
        
    %initialize tempararyError array
        tempararyError = zeros(1,length(invertedTemps));
        
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
                tempararyError = yStdErr;
            end            
        %build fitting data
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>1
                %if not one of the excluded temperature and more then one embryo
            %individual data points
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = repmat(invertedTemps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
            %means as data points
                allMeans(mNum) = mean(data);
                allXM(mNum) = invertedTemps(temp);
                mNum = mNum + 1;
                %fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                    %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2);
            elseif length(data)==1
                %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    %xneg,xpos,'*','Color',colors(3,:),'LineWidth',2);                
            else
                %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    %xneg,xpos,'*','Color',colors(2,:),'LineWidth',2);
            end
        end    
        %copy over errors for monte carlo
            monteTimeErrors(endStage,startStage) = mean(tempararyError); 
            intervalAllX{endStage,startStage} = 1./allXM;
        %Build BIC comparision variables, comparing models, then run
        %through BIC program
        if startStage~=endStage
        %store data points for AIC
        %all points as replicates
        	%{
        	AICTable.numPoints = length(allY);
            AICTable.yi = num2cell(allY,1);
            AICTable.xi = num2cell(allX,1);
         	%}
        %each temp as a replicate
         	AICTable.numPoints = length(allMeans);
         	AICTable.yi = num2cell(allMeans,1);
          	AICTable.xi = num2cell(allXM,1);

        %AIC Comparison calculations, n=1 is assumed best fit (row stages, 
    	%columns are poly compared to) 
          	stage = 1;
                
          	%points as replicates
          	%{
            	polyNMdl1 = fit(allX, allY, 'poly1');%start = 12, end = 14, temp = 2
                	%problem is in scores that arent actually considered
             	polyNMdl2 = fit(allX, allY, 'poly2');%there is some complex warning here...
           	%}
          	%means as replicates
            for fitDegree = 1:8
                    %initialize x (coefficients) as empty 1x9
                    lb = repmat(-Inf,fitDegree,1);
                    ub = [0;repmat(Inf,fitDegree-1,1)];
                        %limit to concave or convex, up or down so it can
                        %converge
                switch fitDegree
                    case 1 
                        fun = @(x,allXM)x(1)^1*allXM + x(2);
                    case 2 
                        fun = @(x,allXM)x(1)^2*allXM + x(2)^1*allXM + x(3);
                    case 3 
                        fun = @(x,allXM)x(1)^3*allXM + x(2)^2*allXM + x(3)^1*allXM + x(4);
                    case 4 
                        fun = @(x,allXM)x(1)^4*allXM + x(2)^3*allXM + x(3)^2*allXM + x(4)^1*allXM + x(5);
                    case 5 
                        fun = @(x,allXM)x(1)^5*allXM + x(2)^4*allXM + x(3)^3*allXM + x(4)^2*allXM + x(5)^1*allXM + x(6);     
                    case 6 
                        fun = @(x,allXM)x(1)^6*allXM + x(2)^5*allXM + x(3)^4*allXM + x(4)^3*allXM + x(5)^2*allXM + x(6)^1*allXM + x(7);                         
                    case 7 
                        fun = @(x,allXM)x(1)^7*allXM + x(2)^6*allXM + x(3)^5*allXM + x(4)^4*allXM + x(5)^3*allXM + x(6)^2*allXM + x(7)^1*allXM + x(8);    
                    case 8 
                        fun = @(x,allXM)x(1)^8*allXM + x(2)^7*allXM + x(3)^6*allXM + x(4)^5*allXM + x(5)^4*allXM + x(6)^3*allXM + x(7)^2*allXM + x(8)^1*allXM + x(9);                                    
                end
                    %will need to flip fun to work properly with calcBIC
                x0 = ones(fitDegree+1,1);
            	mdl = lsqcurvefit(fun,x0,allXM,allMeans);%start = 12, end = 14, temp = 2
                	%problem is in scores that arent actually considered           
                BICval = calcBIC(AICTable,mdl');
                BICs(fitDegree) = BICval;
            end 

        %Store BIC comp values and Eas (just in case)
            sortedBICs = sort(BICs);
                %sort then find 1, 2, 3 to rank order
            rank1(endStage,startStage) = find(BICs == sortedBICs(1));
            rank2(endStage,startStage) = find(BICs == sortedBICs(2));
            rank3(endStage,startStage) = find(BICs == sortedBICs(3));
        else
        end
    end
end

%plot the relevent BIC data range and the associated stage abbreviations
    figure('position',[100 100 1000 800])
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    h = heatmap(xlabels, ylabels, round(rank1(2:12,1:11),1));

%Aestetics
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    x = (0:0.001:1);
    map = colors;
    colormap(h,map)
    caxis([1 8])
    set(gca,'FontSize', 24);
    xlabel("Start Event")
    ylabel("End Event")
    
    %determine temperatures used and adjust title accordingly
        excul = setxor(1:16, excludedTemps);
        title({'\fontsize{24}BIC comparision polynomial Rank 1',strcat("\fontsize{38}Temps ",num2str(min(excul)),":",num2str(max(excul)))});

%plot the relevent BIC data range and the associated stage abbreviations
    figure('position',[100 100 1000 800])
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    h = heatmap(xlabels, ylabels, round(rank2(2:12,1:11),1));

%Aestetics
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    x = (0:0.001:1);
    map = colors;
    colormap(h,map)
    caxis([1 8])
    set(gca,'FontSize', 24);
    xlabel("Start Event")
    ylabel("End Event")
    
    %determine temperatures used and adjust title accordingly
        excul = setxor(1:16, excludedTemps);
        title({'\fontsize{24}BIC comparision polynomial Rank 2',strcat("\fontsize{38}Temps ",num2str(min(excul)),":",num2str(max(excul)))});

%plot the relevent BIC data range and the associated stage abbreviations
    figure('position',[100 100 1000 800])
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    h = heatmap(xlabels, ylabels, round(rank3(2:12,1:11),1));

%Aestetics
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    x = (0:0.001:1);
    map = colors;
    colormap(h,map)
    caxis([1 8])
    set(gca,'FontSize', 24);
    xlabel("Start Event")
    ylabel("End Event")
    
    %determine temperatures used and adjust title accordingly
        excul = setxor(1:16, excludedTemps);
        title({'\fontsize{24}BIC comparision polynomial Rank 3',strcat("\fontsize{38}Temps ",num2str(min(excul)),":",num2str(max(excul)))});

%%
%combining In silico and biological
%run only after running stage-wise code-chunk
figure('Position',[100 100 600 450]);

%remove T0 time point from analysis
    obsScores = (2:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    excludedTemps = [1 13 sizeData-3 sizeData-2 sizeData-1 sizeData];
%select stage of interest
    stage = 11;

%calculate and plot silico prediction first

    %parse constants and data from stage figure
        rangeT = [Temps(2) Temps(length(Temps)-2)];
            %only pull the temperatures we fit above
        xrange = 1./[303.15; 283.15];
        Ea = stageActEa(1:stage).';
            %only pull the Eas of stages we will predict on
        lnA = log(stageActA(1:stage)).';
        
    %set up fuction to predict data
        resolutionT = 160;
        temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';
        fcn = writelnKFcnVectorized(Ea,lnA);
        
    %calculate predicted data and plot
        lnKs = double(fcn(temps));
        
    %fit silico data and plot
        mdl = polyfit(1./temps,lnKs,1);
        fitEval = polyval(mdl,xrange);
        lF = plot(xrange,fitEval,'b','LineWidth',2.5); 
    
    %parse biological experiment data
        pseudoKperStage = cell(length(accumulatedInterval),1);
        for temp = 1:length(accumulatedInterval)
            pseudoKperStage{temp} = 1./accumulatedInterval{temp};
        end
        invertedTemps = 1./Temps;
        xmin = 1/303.15;
        xmax = 1/283.15;
        xrange = [xmin; xmax];
        
    %convert to ln space
        lnPseudoKperStage = cell(length(pseudoKperStage),1);
        for temp = 1:length(pseudoKperStage)
            lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
        end
        
    %calc total num of embryos
        fittableMean = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 0
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
            ci95 = tinv([0.025 0.975], N-1); 
                %95 confidence interval
            yci95 = yStdErr.*ci95;
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 1.96*0.5;
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
        %plot error bar
            meanLnPseudoK = mean(data);
       	%build all data for specific temp for linear regression
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>0
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
        
    %calculate and plot biological linear regression
     	mdl2 = polyfit(allXM, allMeans,1);
      	yLin = polyval(mdl2,xrange);
      	lm = plot(xrange,yLin, 'b', 'LineWidth',2.5);
        
        lF = plot(xrange,fitEval,'--c','LineWidth',2.5); 
        
    %Eas
        text(1./(273.15+17.35),-11.05,strcat("E_a = ",num2str(round((-mdl2(1)*R)*tokJConvert,0))," kJ/mol"),'Color','b','FontSize',24)
        text(1./(273.15+17.35),-11.35,strcat("E_a = ",num2str(round((-mdl(1)*R)*tokJConvert,0))," kJ/mol"),'Color','c','FontSize',24)
        
    %beautify
        %ylim([mean(allMeans)-1.25 mean(allMeans)+1.25])
        ytLabels = yticks;
        yticks(ytLabels(floor(ytLabels)==ytLabels));
            %present only the interger yticks
        ylabel('ln(1/time [s^{-1})]');
        xlim(xrange);
        xlabel('1/(T [^oC] + 273.15)');
        xticks(1./[303.1500 293.1500 283.1500]);
        xlim(xrange);
        xticklabels(strcat('1/',strsplit(num2str(round((1./xticks-273.15),1)))));
        set(gca,'linewidth', 3, 'FontSize', 30);   
        
    %seperate to avoid being affected by gca
        title({'Frog Biological vs Prediction';strcat("A to ",stageAbbs(stage))}, 'FontSize', 25)
        legend([eb lm lF],{'Data','Biological Fit','Predictive Fit'},'FontSize',24)
%%
%Quadratic
%combining In silico and biological
%run only after running stage-wise code-chunk
figure('Position',[100 100 600 450]);

%remove T0 time point from analysis
    obsScores = (2:12);
    accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
    excludedTemps = [12 13];%[1 length(lnPseudoKperStage)-1 length(lnPseudoKperStage)];% 0;
%select stage of interest
    stage = 8;

%calculate and plot silico prediction first

    %parse constants and data from stage figure
        rangeT = [Temps(2) Temps(length(Temps)-2)];
            %only pull the temperatures we fit above
        xrange = 1./[307.15; 280.15];
        xtemps = min(xrange):0.00001:max(xrange);
        Ea = stageQuadEa(1:stage);
            %only pull the Eas of stages we will predict on
        lnA = stageQuadA(1:stage);
        rou = stageRou(1:stage);
        
    %set up fuction to predict data
        resolutionT = 160;
        temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';
        fcn = writeQuadlnKFcnVectorized(Ea,lnA,rou);
        
    %calculate predicted data and plot
        lnKs = double(fcn(temps));
        
    %fit silico data and plot, quadratic
        mdl = polyfit(1./temps,lnKs,2);
        fitEval = polyval(mdl,xtemps);
        %extreme range
        lF = plot(xtemps,fitEval,'m','LineWidth',2);
        %linear range
    
    %parse biological experiment data
        pseudoKperStage = cell(length(accumulatedInterval),1);
        for temp = 1:length(accumulatedInterval)
            pseudoKperStage{temp} = 1./accumulatedInterval{temp};
        end
        invertedTemps = 1./Temps;
        xmin = 1/307.15;
        xmax = 1/280.15;
        xrange = [xmin; xmax];
        
    %convert to ln space
        lnPseudoKperStage = cell(length(pseudoKperStage),1);
        for temp = 1:length(pseudoKperStage)
            lnPseudoKperStage{temp} = log(pseudoKperStage{temp});
        end
        
    %calc total num of embryos
        fittabeEmNum = 0;
        fittableMean = 0;
        for temp = 1:length(lnPseudoKperStage)
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                if sum(~isnan(lnPseudoKperStage{temp}(:,stage))) > 1
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
            ci95 = tinv([0.025 0.975], N-1); 
                %95 confidence interval
            yci95 = yStdErr.*ci95;
            yneg = yci95(1);
            ypos = yci95(2);
            errTemp = 1.96*0.5;
            xneg = invertedTemps(temp)-1./(Temps(temp)+errTemp);
            xpos = xneg;
        %plot error bar
            meanLnPseudoK = mean(data);
       	%build all data for specific temp for linear regression
            if sum(temp ~= excludedTemps)/length(excludedTemps) == 1&&length(data)>1
                nextEmNum = emNum+size(data,1);
                allY(emNum:nextEmNum-1,1) = data;
                allX(emNum:nextEmNum-1,1) = 1./repmat(Temps(temp),...
                    size(data,1),1);
                emNum = nextEmNum;
                eb = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                    xneg,xpos,'bo','LineWidth',2);
                allMeans(mNum) = meanLnPseudoK;
                allXM(mNum) = invertedTemps(temp);
                mNum = mNum + 1;
            else
                errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,xneg,...
                    xpos,'r*','LineWidth',2);
            end 
         end
        
    %calculate and plot biological linear regression
     	mdl2 = polyfit(allX, allY,2);
      	yLin = polyval(mdl2,xtemps);
      	lm = plot(xtemps,yLin, '--b', 'LineWidth',2);
        
    %beautify
        %ylim([mean(allMeans)-2 mean(allMeans)+1.25])
        ytLabels = yticks;
        yticks(ytLabels(floor(ytLabels)==ytLabels));
            %present only the interger yticks
        ylabel('ln(k)');
        xlim(xrange);
        xlabel('1/T (C^{-1})');
        xticks(1./[303.1500 293.1500 283.1500]);
        xlim(xrange);
        xticklabels(strcat('1/',strsplit(num2str(round((1./xticks-273.15),1)))));
        set(gca,'linewidth', 2, 'FontSize', 30);   
        
    %seperate to avoid being affected by gca
        title({'Frog Biological vs Prediction';strcat("a to ",stageAbbs(stage+1))}, 'FontSize', 25)
        legend([lF lm eb],{'Predictive Fit','Biological Fit','Biological Data'},'FontSize',20)

%%
%Comparing old N and new N
%load in new N data first
%Data Import

%DataSet1 (New)
%Destination
    %filename1 = 'FrogScoresFinalNewData.xlsx';
    filename1 = 'FrogScoresFinal_4.xlsx';
    folder1 = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Frog Scores/';

%DataSet2 (Old)
%Destination      
    
    filename2 = 'FrogScoresFinalNewData.xlsx';
    folder2 = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Frog Scores/';
    
%DS1, new
%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder1,filename1));
    %[~,sheet_name]=xlsfinfo(strcat(filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);
    data = readcell(strcat(folder1,filename1),'Sheet', sheet_name{1});
    %data = readcell(strcat(filename),'Sheet', sheet_name{1});

%find number of events scored ("scores") in file
    isSearch = cellfun(@(x)isequal(x,'Score Number'),data);
    [row,col] = find(isSearch);
    newTemps = zeros(sizeData,1);
    numScores = max(cell2mat(data(row,(col+2:size(data,2)-3))))+2;
        %+1 to include T0 as a score
    
%Initialize arrays to store different data types        
    accumulatedNew = cell(sizeData,1);
    perStageNew = cell(sizeData,1);

%iterate over each sheet to parse data types into their respective arrays        
for i = 1:sizeData
%read in specific data of sheet i
    data = readcell(strcat(folder1,filename1),'Sheet', sheet_name{i});
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
    %parse out temperature data for each sheet and store in array
    isSearch = cellfun(@(x)isequal(x,'Temp'),data);
    [row,col] = find(isSearch);
    newTemps(i) = cell2mat(data(row,(col+1)));
end    
    

%DS2, old  

%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder2,filename2));
    %[~,sheet_name]=xlsfinfo(strcat(filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);
    data = readcell(strcat(folder2,filename2),'Sheet', sheet_name{1});
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
    data = readcell(strcat(folder2,filename2),'Sheet', sheet_name{i});
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
%for any stage-to-stage interval
%Arrhenius Stage-wise
%Stages
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;

%Stages to observe
    startStages = [1 8];
    endStages = 9:10;
    
%Plot Stages
    colors = [0 0 1; 1 0 0; 0 1 0];
        %col for fit, good data, dropped data
    
%remove T0 time point from analysis
    %perStageIntervalNew = cellfun(@(x) x(:,obsScores),perStageNew,'UniformOutput',false);
    %perStageIntervalOld = cellfun(@(x) x(:,obsScores),perStageOld,'UniformOutput',false);
    
    obsScores = (1:numScores);
    stageAbbsNoT0 = stageAbbs(obsScores);  
    stageNamesNoT0 = stageNames(obsScores);
    excludedCombinedTemps = [1 16 17 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];%[12 13]
    excludedTemps = [1 11 14:16 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];
        %for original Final data [1 16 17 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];
        %for combined data [1 11 14:16 sizeData-4 sizeData-3 sizeData-2 sizeData-1 sizeData];
    
%Initialize tables
    %Statistical Table
        compRegs = zeros(2,8); %2 cause comparing 2 data sets for each stage
        compRegs = array2table(compRegs);
        compRegs.Properties.VariableNames = {'numPoints','SX','SY',...
            'SXX','SYY','SXY','AveX','AveY'};
        compRegs.Properties.RowNames = ["New" "Old"];

%invert temps and temp range        
    invertedTemps = 1./Temps;
    invertedNewTemps = 1./newTemps;
    xmin = 1/307.15;
    xmax = 1/280.15;
    xrange = [xmin; xmax];

%convert raw data to log(inverse(t))    
    accumulatedNewInterval = cellfun(@(x) x(:,obsScores),accumulatedNew,'UniformOutput',false);
    accumulatedOldInterval = cellfun(@(x) x(:,obsScores),accumulatedOld,'UniformOutput',false);
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
   % standErr = nan(length(invertedTemps),length(lnPseudoKperStageNew{1}));

%plot old first then new over it
for startStage = startStages
    for endStage = endStages
        %initialize and store subplot data
        figure('Position',[100 100 600 450]); 
        if startStage == endStage
        else
        %calculate interval times by diff between acc(stages)
        dataRange = cellfun(@(x) x(:,[startStage endStage]),accumulatedOldInterval,'UniformOutput',false);
        diffDataRange = cellfun(@(x) diff(x,1,2),dataRange,'UniformOutput',false);
        %change any zeros to Nans
        convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'[^0123456789.]0','NaN')).',diffDataRange,'UniformOutput',false);
        convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'^0[^0123456789.]','NaN')).',convertZeros,'UniformOutput',false);        
        %convert times to Arrhenius ln(k)
        lnPseudoKperStageOld = cellfun(@(x) log(1./x),convertZeros,'UniformOutput',false);
         
    
        %Plot new N data (want new ploted first so old shows up overlapping)
            %calc total num of fittable embryos per stage for all temps, 2 or more
            %reps
                fittabeEmNum = 0;
                fittableMean = 0;
                overAllEms = 0;
                for temp = 1:length(lnPseudoKperStageOld)
                    if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                        if sum(~isnan(lnPseudoKperStageOld{temp})) > 0% 1
                            fittabeEmNum = fittabeEmNum + ...
                                sum(~isnan(lnPseudoKperStageOld{temp}));
                            fittableMean = fittableMean + 1;                    
                        end
                    end
                    overAllEms = overAllEms + sum(~isnan(lnPseudoKperStageOld{temp}));
                end 

            %initialize arrays to hold all fittable points    
                allX = nan(fittabeEmNum,1);
                allY = nan(fittabeEmNum,1);  
                emNum = 1;

            %iterate over temperaters to average embryos for each temp then plot
            for temp = 1:length(invertedTemps)
                hold on;

                %parse temperary 'data' of nans
                    data = lnPseudoKperStageOld{temp};
                    data = data(~isnan(data));

                %set/calculate errors
                    N = size(data,1); 
                        %number non-nan experiments in data
                    yStdErr = std(data)/sqrt(N-1); 
                        %standard Error
                    ci95 = tinv([0.025 0.975], N-1); 
                        %95 confidence interval
                    yci95 = yStdErr.*ci95;
                    yneg = yci95(1);
                    ypos = yci95(2);
                    errTemp = 1.96*0.5; 
                        %1.96 is multiplier associated with a two sided confidence interval
                    xneg = invertedTemps(temp)-1./(newTemps(temp)+errTemp);
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
                        %fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                            %ypos,xneg,xpos,'o','Color',colors(1,:),'LineWidth',2.5);
                        fittedOld = scatter(repmat(invertedTemps(temp),1,length(data)),data,[],colors(1,:),'o','LineWidth',2.5);
                    elseif length(data)==1
                        %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                         %   xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);
                        unfittedOld = scatter(repmat(invertedTemps(temp),1,length(data)),data,[],colors(2,:),'*','LineWidth',2.5);
                    else
                        %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                         %   xneg,xpos,'*','Color',colors(2,:),'LineWidth',2.5);
                        scatter(repmat(invertedTemps(temp),1,length(data)),data,[],colors(2,:),'*','LineWidth',2.5);
                    end

            end
            %combine allX/Y
                allData = horzcat(allX,allY);

            %store regression stats
                %individual points

                    compRegs.numPoints(1) = length(allData);
                    compRegs.SX(1) = real(sum(allData(:,1)));
                    compRegs.SY(1) = real(sum(allData(:,2)));
                    compRegs.SXX(1) = real(sum(allData(:,1).^2));
                    compRegs.SYY(1) = real(sum(allData(:,2).^2));
                    compRegs.SXY(1) = real(sum(allData(:,1).*allData(:,2)));
                    compRegs.AveX(1) = real(mean(allData(:,1)));
                    compRegs.AveY(1) = real(mean(allData(:,2)));

            %plot linear regression
                lineFit = fit(allData(:,1), allData(:,2), 'poly1');
                Ea1 = round((-lineFit.p1*R)*tokJConvert,0);

            %plot fit lines
                yLin = polyval([lineFit.p1 lineFit.p2],xrange);
                lm = plot(xrange,yLin, 'b', 'LineWidth',2.5);


        %Plot new data on top
        %calculate interval times by diff between acc(stages)
            dataRange = cellfun(@(x) x(:,[startStage endStage]),accumulatedNewInterval,'UniformOutput',false);
            diffDataRange = cellfun(@(x) diff(x,1,2),dataRange,'UniformOutput',false);
            %change any zeros to Nans
            convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'[^0123456789.]0','NaN')).',diffDataRange,'UniformOutput',false);
            convertZeros = cellfun(@(x) str2num(regexprep(num2str(x.'),'^0[^0123456789.]','NaN')).',convertZeros,'UniformOutput',false);        
            %convert times to Arrhenius ln(k)
            lnPseudoKperStageNew = cellfun(@(x) log(1./x),convertZeros,'UniformOutput',false);
        %Plot new N data (want new ploted first so old shows up overlapping)
            %calc total num of fittable embryos per stage for all temps, 2 or more
            %reps
                fittabeEmNum = 0;
                fittableMean = 0;
                overAllEms = 0;
                for temp = 1:length(lnPseudoKperStageNew)
                    if sum(temp ~= excludedCombinedTemps)/length(excludedCombinedTemps) == 1
                        if sum(~isnan(lnPseudoKperStageNew{temp})) > 0
                            fittabeEmNum = fittabeEmNum + ...
                                sum(~isnan(lnPseudoKperStageNew{temp}));
                            fittableMean = fittableMean + 1;                    
                        end
                    end
                    overAllEms = overAllEms + sum(~isnan(lnPseudoKperStageNew{temp}));
                end 

            %initialize arrays to hold all fittable points    
                allX = nan(fittabeEmNum,1);
                allY = nan(fittabeEmNum,1);  
                emNum = 1;

            %iterate over temperaters to average embryos for each temp then plot
            for temp = 1:length(invertedNewTemps)
                hold on;

                %parse temperary 'data' of nans
                    data = lnPseudoKperStageNew{temp};
                    data = data(~isnan(data));

                %set/calculate errors
                    N = size(data,1); 
                        %number non-nan experiments in data
                    yStdErr = std(data)/sqrt(N-1); 
                        %standard Error
                    ci95 = tinv([0.025 0.975], N-1); 
                        %95 confidence interval
                    yci95 = yStdErr.*ci95;
                    yneg = yci95(1);
                    ypos = yci95(2);
                    errTemp = 1.96*0.5; 
                        %1.96 is multiplier associated with a two sided confidence interval
                    xneg = invertedNewTemps(temp)-1./(newTemps(temp)+errTemp);
                    xpos = xneg;

                %build all data for specific temp for linear regression, while
                %plotting error bar
                    meanLnPseudoK = mean(data);
                    if isnan(meanLnPseudoK)==1
                    elseif sum(temp ~= excludedCombinedTemps)/length(excludedCombinedTemps) == 1&&length(data)>0%1
                        nextEmNum = emNum+size(data,1);
                        allY(emNum:nextEmNum-1,1) = data;
                        allX(emNum:nextEmNum-1,1) = repmat(invertedNewTemps(temp),...
                            size(data,1),1);
                        emNum = nextEmNum;
                    %fit to individuals
                        %fitted = errorbar(invertedTemps(temp),meanLnPseudoK,yneg,...
                         %   ypos,xneg,xpos,'o','Color',[0 1 1],'LineWidth',2.5);
                        fittedNew = scatter(repmat(invertedNewTemps(temp),1,length(data)),data,[],[0 1 1],'o','LineWidth',2.5);
                    elseif length(data)==1
                        %errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                         %   xneg,xpos,'*','Color',[1 0.5 0],'LineWidth',2.5); 
                        unfittedNew = scatter(repmat(invertedNewTemps(temp),1,length(data)),data,[],[1 0.5 0],'*','LineWidth',2.5); 
                    else
                       % errorbar(invertedTemps(temp),meanLnPseudoK,yneg,ypos,...
                       %     xneg,xpos,'*','Color',[1 0.5 0],'LineWidth',2.5);
                        scatter(repmat(invertedNewTemps(temp),1,length(data)),data,[],[1 0.5 0],'*','LineWidth',2.5);
                    end

            end
            %combine allX/Y
                allData = horzcat(allX,allY);

            %store regression stats
                %individual points

                    compRegs.numPoints(2) = length(allData);
                    compRegs.SX(2) = real(sum(allData(:,1)));
                    compRegs.SY(2) = real(sum(allData(:,2)));
                    compRegs.SXX(2) = real(sum(allData(:,1).^2));
                    compRegs.SYY(2) = real(sum(allData(:,2).^2));
                    compRegs.SXY(2) = real(sum(allData(:,1).*allData(:,2)));
                    compRegs.AveX(2) = real(mean(allData(:,1)));
                    compRegs.AveY(2) = real(mean(allData(:,2)));

            %plot linear regression
                lineFit = fit(allData(:,1), allData(:,2), 'poly1');
                Ea2 = round((-lineFit.p1*R)*tokJConvert,0);

            %plot fit lines
                yLin = polyval([lineFit.p1 lineFit.p2],xrange);
                lm = plot(xrange,yLin, 'c', 'LineWidth',2.5);



        %perform statistical analysis            
            pVal = ancova(compRegs);
            text(1./(273.15+30),-5.25,strcat("E_a = "),'Color','k','FontSize',19)
            text(1./(273.15+27),-5.1,strcat(num2str(Ea1)),'Color','b','FontSize',19)
            text(1./(273.15+25.5),-5.1,strcat(",",num2str(Ea2)),'Color','c','FontSize',19)
            text(1./(273.15+23.5),-5.1,strcat(" kJ/mol"),'Color','k','FontSize',19)
            text(1./(273.15+27),-6,strcat("pVal = ",num2str(round(pVal,3))),'Color','k','FontSize',19)


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
                
            ylabel('ln(1/time [s^{-1})]');
            xlabel('1/(T [^oC] + 273.15)');
            
            %dummies
            unfittedOld = scatter(0,0,[],colors(2,:),'*','LineWidth',2.5);
            unfittedNew = scatter(0,0,[],[1 0.5 0],'*','LineWidth',2.5); 


            legend([fittedOld unfittedOld fittedNew unfittedNew],{'fitted New', 'unfitted New', 'fitted Original', 'unfitted Original'},'Location','northeast')
            %{
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
            %}
            %set gca to be legible    
                set(gca,'linewidth', 3, 'FontSize', 23);
                    %28

            %place titles as text inside plots
                    text(mean(xrange),-11.9,{(strcat(stageAbbsNoT0(startStage),...
                        " to ",(stageAbbsNoT0(endStage))))}, 'FontSize', 24,...
                        'HorizontalAlignment','center');
                        %" to ",(stageAbbsNoT0(stage))))}, 'FontSize', 28,...
                        %'HorizontalAlignment','center');
        end
    end
end


%{
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
%}       