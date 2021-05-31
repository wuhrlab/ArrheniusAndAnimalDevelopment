%%
%For comparing the Eisen lab's data conclusions to our own

%%
%Load in the data
newTable = readtable('/Users/josephcrapse/Desktop/ArrheniusPaper/Revisions Round 2/EisenSourceData/TimeLapseInstructions_run_trial15.txt','DatetimeType','text');%'Format','%q%q%q%f%q%f%f%q%f%f%f');
newTable.Properties.VariableNames = {'date', 'position', 'firstFilename', 'orientation',...
    'fileNameStructure', 'mryTime', 'tfTime', 'strain', 'temp', 'scaling',...
    'zStack'};
    %key 
        %V1 == date shot
        %V2 == position on said date
        %V3 == some jpg filename
        %V4 == orientation
        %V5 == filename structure
        %V6 == membrane_reaches_yolk time
        %V7 == trachea_fill timw
        %V8 == species/strain
        %V9 == temperature recorded at
        %V10 == dilation (scaling)
        %V11 == zStack
%determine which data is actually D. melanogaster then copy over to a table
isMe = cellfun(@(x) strcmp(string(x),"Me"),newTable.strain);
meTable = newTable(isMe,:);
%remove all zStack duplicates (zStack != 0)
meTable = meTable(meTable.zStack == 0,:);

meTable.date = cellstr(meTable.date);
sizeData = size(meTable);

%adjust the dates so they match the date format from data below

%combine the data and position columns to match second batch of data below
%badDates = cellfun(@(x) contains(x,'00'), meTable.date);
%meTable.date(badDates) = cellfun(@(x) x(3:length(x)),meTable.date(badDates),'UniformOutput',false);
date_position = cellfun(@(x,y) strcat(x,'_',y), meTable.date,meTable.position,'UniformOutput',false);
%variable names matching the second batch of data later below
varNames = {'conditions','event','date_position','finalTimepoint','dilation'};
%build each stage's table with the necessary data and identifiers
eventMry = repmat("membrane_reaches_yolk",sizeData(1),1);
yolkT = table(meTable.temp, eventMry, date_position, meTable.mryTime, meTable.scaling);
eventTf = repmat("trachea_fill",sizeData(1),1);
trachT = table(meTable.temp, eventTf, date_position, meTable.tfTime, meTable.scaling);

%give new variable names, kinda unecessary, these get changed later
yolkT.Properties.VariableNames = varNames;
trachT.Properties.VariableNames = varNames;

%simplify the tables to a date+pos+temp identifier column and an absolute
%time column
date_pos_temp = cellfun(@(x,y) strcat(x,'#',y), yolkT.date_position,yolkT.conditions,'UniformOutput',false);
absoluteTime = yolkT.finalTimepoint.*yolkT.dilation;
yolkT = table(date_pos_temp,absoluteTime);

date_pos_temp = cellfun(@(x,y) strcat(x,'#',y), trachT.date_position,trachT.conditions,'UniformOutput',false);
absoluteTime = trachT.finalTimepoint.*trachT.dilation;
trachT = table(date_pos_temp,absoluteTime);

%load data for the second 'batch' of data (different format than the first)
mainFolder = '/Users/josephcrapse/Desktop/ArrheniusPaper/Revisions Round 2/EisenSourceData/';
stages = {'pole_bud_appears','pole_cell_invagination','amnioproctodeal_invagination','amnioserosa','clypeolabrum_retracts','heart-shaped_midgut'};
abvs = {'pba', 'pci', 'ai', 'a', 'cr','hsm'};
%cycle through each stage to build a table from the respective folder
for s = 1:length(stages)
    %construct the directory for the specific folder for the current stage
    direct = dir(char(strcat(mainFolder,stages(s))));
    directT = struct2table(direct);
    %slect files with Me only (these are D. melanogaster) and save
    %filenames
    isMe = cellfun(@(x) contains(string(x),"Me"),directT.name);
    fileNames = directT.name(isMe);
    if s == 5
        %specifically for 5 need to remove last file, 35 degrees, as it is
        %empty and an unused temperature
        fileNames = fileNames(1:length(fileNames)-1);
    end
    %loop through the files for each temperature at the current stage and
    %combine them into one table
    for i = 1:length(fileNames)
        if i == 1
            midT = readtable(char(strcat(mainFolder,stages(s),'/',fileNames(i))));
        else
            tempT = readtable(char(strcat(mainFolder,stages(s),'/',fileNames(i))));
            midT = vertcat(tempT,midT);
        end
    end
    %modifications to 'cleanup' or simplify, by removing unnecessary text
    %and combining date, pos, temp into one identifer
    midT.conditions = cellfun(@(x) (erase(x,'Me')),midT.conditions,'UniformOutput',false);
    date_pos_temp = cellfun(@(x,y) strcat(x,'#',y), midT.date_position,midT.conditions,'UniformOutput',false);
    %modifications to get specific time in minutes and move key variable to
    %left, then build table
    absoluteTime = midT.finalTimepoint.*midT.dilation;
    finalT = table(date_pos_temp, absoluteTime);%table(midT.date_position, midT.conditions, absoluteTime);%
    finalT.Properties.VariableNames = {'date_pos_temp', strcat(abvs{s},'Time')};
    
    %save the temparary table into specific tables for each stage
    switch s
        case 1
            poleAppT = finalT;
        case 2
            poleInvT = finalT;
        case 3
            amnInvT = finalT;
        case 4
            amnT = finalT;
        case 5
            clyRetT= finalT;
        case 6
            heartT = finalT;
    end       
end

%%
%now that tables are organized, reorder data in same manner as mine for
%likewaise analysis, also need to reduce to lowest common denominator for
%unique identifiers
%check highest unique denom
    %check if there are non-unique embryo markers (there was 1, so had to
    %drop both repeats)
    stages = {'pole_bud_appears','membrane_reaches_yolk','pole_cell_invagination','amnioproctodeal_invagination','amnioserosa','clypeolabrum_retracts','heart-shaped_midgut','trachea_fill'};
    abvs = {'pba', 'mry', 'pci', 'ai', 'a', 'cr','hsm','tf'};
    %initialize data table, from one arbitrary table   
    data = table(finalT.date_pos_temp);
    data.Properties.VariableNames = {'date_pos_temp'};
    %copy over each table to a temparary
    for i = 1:length(stages)
        switch i 
            case 1 
                temparayT = poleAppT;
            case 2
                temparayT = yolkT;
            case 3
                temparayT = poleInvT;
            case 4
                temparayT = amnInvT;
            case 5
                temparayT = amnT;
            case 6
                temparayT = clyRetT;
            case 7
                temparayT = heartT;
            case 8
                temparayT = trachT;
        end
        %Check each stage's table for unique embryo identifiers indx
        [~, uniques] = unique(temparayT.date_pos_temp);
        nonUnique = setdiff(1:length(temparayT.date_pos_temp),uniques);
        %save the name of non unique identifiers
        nonUnEm = temparayT.date_pos_temp(nonUnique);
        em2Delete = [];
        %find these identifiers in the data table, if they exist
        for num = 1:length(nonUnEm)
            tempArr = find(cell2mat(cellfun(@(x) strcmp(x, nonUnEm(num)),data.date_pos_temp,'UniformOutput',false)));
            em2Delete =  vertcat(em2Delete,tempArr);
        end
        %delete repeats
        data(em2Delete,:) = [];
    end
    
    %loop through all stages and modify the data table with a new variable
    %to be filled with that stages timings, date_pos_temp as identifier
    for s = 1:length(stages)
        data = addvars(data, nan(length(data.date_pos_temp),1),'NewVariableNames', strcat(abvs{s},'Time'));
        indx = nan(length(data.date_pos_temp),1);
        tempArray = nan(length(data.date_pos_temp),1);
            switch s
                case 1
                    %loop through all the different unique embryos decided
                    %on above
                    for i = 1:length(data.date_pos_temp)
                        %find the indx for the (i)th embryo
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),poleAppT.date_pos_temp,'UniformOutput',false)));
                        %copy over the time value if the embryo exists in
                        %the specific stage's table at the (i)th position
                        %else leave the (i)th position as a Nan
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = poleAppT.pbaTime(tempVal);
                        else
                        end
                    end
                    %copy over all the time information into the
                    %appropriate stage's position in the master table
                    data(:,s+1) = array2table(tempArray);
                case 2
                    for i = 1:length(data.date_pos_temp)
                        %cannot us contains, as it gives multiple hits for
                        %unique responses, need to edit the original
                        %datetime to remove 2 leading zeros
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),yolkT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = yolkT.absoluteTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
                case 3
                    for i = 1:length(data.date_pos_temp)
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),poleInvT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = poleInvT.pciTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
                case 4
                    for i = 1:length(data.date_pos_temp)
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),amnInvT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = amnInvT.aiTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
                case 5
                    for i = 1:length(data.date_pos_temp)
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),amnT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = amnT.aTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
                case 6
                    for i = 1:length(data.date_pos_temp)
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),clyRetT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = clyRetT.crTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
                case 7
                    for i = 1:length(data.date_pos_temp)
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),heartT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = heartT.hsmTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
                case 8
                    for i = 1:length(data.date_pos_temp)
                        tempVal = find(cell2mat(cellfun(@(x) strcmp(x, data.date_pos_temp{i}),trachT.date_pos_temp,'UniformOutput',false)));
                        if isempty(tempVal) == 0
                            indx(i,:) = tempVal;
                            tempArray(i,:) = trachT.absoluteTime(tempVal);
                        else
                        end
                    end
                    data(:,s+1) = array2table(tempArray);
            end
    end
    
%split identifer (dont need it anymore) into date_pos and temp
splitDim = split(data.date_pos_temp,'#');
date_pos = splitDim(:,1);
temps = str2double(splitDim(:,2));
data = addvars(data(:,2:9),date_pos,temps,'Before','pbaTime');

%normalize times to yolk stage (Acc) or previous stage (Stg)
%first initialize new tables
dataAcc = data(:,1:9);
dataStg = data(:,1:9);
stgDurs = strcat(abvs(1:7),'_2',abvs(2:8));
accDurs = [strcat(abvs(1),'_2',abvs(2)),strcat(abvs(2),'_2',abvs(3:8))];
for s = 1:length(stages)-1
    
    if s == 1
        %order is reversed to make it positive, since this is T-1 and T0
        dataAcc(:,s+2) = table(data.mryTime - data{:,s+2});
        %dataStg(:,s+2) = table(dataStg{:,s+2} - dataStg{:,s+1});
        %no real T0 for stage, so can calculate all the same
        %calc pba -> mry same as all others
        dataStg(:,s+2) = table(data{:,s+3} - data{:,s+2});
    else
        dataAcc(:,s+2) = table(data{:,s+3} - data.mryTime);
    	dataStg(:,s+2) = table(data{:,s+3} - data{:,s+2});
    end
    dataAcc.Properties.VariableNames(3:9) = accDurs;
    dataStg.Properties.VariableNames(3:9) = stgDurs;
end

dataAbs = data;
for s = 1:length(stages)
   dataAbs(:,s+2) = table(dataAbs{:,s+2} - data.mryTime);
end
%%
%write an excel file with sheets like our fly/frog data
unTemps = unique(temps);
for i = 1:length(unTemps)
    indx = find(data.temps == unTemps(i));
    writetable(data(indx,[1,3:10]),'/Users/josephcrapse/Desktop/KuntzEisenData.xlsx','Sheet',num2str(unTemps(i)))
end

%%
%preliminary analysis
%constants
    R = 8.3144598;
    tokJConvert = 1/1000;
    xrange = 1./([34.5 15.5]+273.15);
    
for type = 1:2   
    switch type
        case 1
            dataType = dataAcc;
            durs = accDurs;
            durs2 = replace(durs, '_2','-');
            figure('Name','Accumulated','Position',[10 10 1200 500])
        case 2
            dataType = dataStg;
            durs = stgDurs;
            durs2 = replace(durs, '_2','-');
            figure('Name','Stagewise','Position',[10 10 1200 500])
            
            stagRegM = zeros(length(stgDurs),8);
            stagRegT = array2table(stagRegM);
            stagRegT.Properties.VariableNames = {'numPoints','SX','SY',...
                'SXX','SYY','SXY','AveX','AveY'};
            stagRegT.Properties.RowNames = stgDurs;
            
            %Initialize fit comparison Table
                AICTable = table('Size',[1,4],'VariableTypes',{'double',...
                    'cell','cell','cell'});
                AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};
    end
    
    for s = 1:length(durs)
        %pullout fittable embryos (ie, no nans)
            fittedEmsIndx = ~isnan(dataType{:,s+2}) & dataType{:,s+2} > 0;
            fittedXOverall = 1./(dataType.temps(fittedEmsIndx)+273.15);
            fittedYOverall = log(1./dataType{fittedEmsIndx,s+2});
                %for some reason there are values at 0 or lower...some out
                %of order stages?!
            sp(s) = subplot(2,4,s);

        %plot embryos
            scatter(fittedXOverall,fittedYOverall,'LineWidth',2)
            hold on

        %fit overall
            xrangeTemps = xrange(1):0.00001:xrange(2);
            fitCurv = fit(fittedXOverall,fittedYOverall,'poly2');
            curv = polyfit(fittedXOverall,fittedYOverall,2);
            lm = polyfit(fittedXOverall,fittedYOverall,1);
            fitEval = polyval([fitCurv.p1 fitCurv.p2 fitCurv.p3],xrangeTemps);
            plot(xrangeTemps,fitEval,'b--','LineWidth',2);
            overallEa(s) = (-fitCurv.p1*R)*tokJConvert;
            ci = confint(fitCurv,0.68);
            overallCi(:,s) = (ci(:,1).*R).*tokJConvert;
            overallErrs(s) = round(abs(overallEa(s) + overallCi(1,s)),1);

        %fit our linear
            fittedEmsIndx = ~isnan(dataType{:,s+2}) & dataType.temps<=27.5 & dataType{:,s+2} > 0;
            fittedX = 1./(dataType.temps(fittedEmsIndx)+273.15);
            fittedY = log(1./dataType{fittedEmsIndx,s+2});
                %fitted embryos are different
                %for some reason there are values at 0 or lower...some out
                %of order stages?!
            fitLm = fit(fittedX,fittedY,'poly1'); 
            fitEval = polyval([fitLm.p1 fitLm.p2],xrange);
            plot(xrange,fitEval,'b','LineWidth',2);
            linearEa(s) = (-fitLm.p1*R)*tokJConvert;   
            ci = confint(fitLm,0.68);
            linearCi(:,s) = (ci(:,1).*R).*tokJConvert;
            linearErrs(s) = round(abs(linearEa(s) + linearCi(1,s)),1);
            
        %store data for stagewise pval analysis
            if type == 2
                stagRegT.numPoints(s) = sum(fittedEmsIndx);
                stagRegT.SX(s) = sum(fittedX);
                stagRegT.SY(s) = sum(fittedY);
                stagRegT.SXX(s) = sum(fittedX.^2);
                stagRegT.SYY(s) = sum(fittedY.^2);
                stagRegT.SXY(s) = sum(fittedX.*fittedY);
                stagRegT.AveX(s) = mean(fittedX);
                stagRegT.AveY(s) = mean(fittedY); 
            end
            
        %calculate 'BIC' for overall 
   
            AICTable.numPoints = length(fittedYOverall);
            AICTable.yi = num2cell(fittedYOverall,1);
            AICTable.xi = num2cell(fittedXOverall,1);
            
            BICcomp = fitBICCompCalc(AICTable,lm,curv,1);            

        %report E_as, and xline
            xline(1./(mean([27.5 30])+273.15),'r--','lineWidth',2);
            text(1./(26.75+273.15), -2.1, strcat(...
                '$ln(\frac{L_{Q}}{L_{L}}) = ',num2str(round(log(BICcomp),...
                2, 'significant')),'$'),'Color','k','FontSize',18,'Interpreter','latex');
            text(1./(28.5+273.15), -1.1, strcat("E_a = ",...
                num2str(round(linearEa(s),2,'significant'))," +/- ",...
                num2str(round(linearErrs(s),2,'significant'))),'Color','blue','FontSize',18);

        %graphics
            if s == 1 || s == 5
                ylabel('ln(1/time [min^{-1}])') 
            else
                set(gca,'yticklabel',[])
            end
            if s >=4
                xlabel('1/(T [^oC] + 273.15)')
            else
                set(gca,'xticklabel',[])
            end
            set(gca,'LineWidth',2,'FontSize',18)
            %title(durs2(s))
            text(mean(xrange),-7.75,durs2(s),'FontSize', 18,'HorizontalAlignment','center');
            xlim(xrange)
            xticks(1./[301.1500 294.1500]);
            xticklabels({'1/28','1/21'});
            switch type
                case 1
                    ylim([-8.5 -0.5])
                case 2
                    ylim([-8.5 -0.5])
            end         
    end
    
                
    %Adjust Subplot positions to tighten image
    %Shift left
        for stage = 1:length(durs)
            if stage ~= 1 && stage ~= 5
                prevPos = get(sp(stage-1),'Position');
                prevPos = [prevPos(1)+prevPos(3)+0.0075 prevPos(2) ...
                    prevPos(3) prevPos(4)];
                set(sp(stage),'Position', prevPos)
            end
        end
    %shift up
        for stage = 1:length(durs)
            if stage >= 5 && stage <=8
                prevPos = get(sp(1),'Position');
                curPos = get(sp(stage),'Position');
                curPos(2) = prevPos(2)-prevPos(4)-0.015;
                set(sp(stage),'Position', curPos)
            end
        end
    %
    %plot figure for comparing overall vs linear Eas
        figure('Position',[10 10 700 400])
        errorbar(overallEa,overallErrs,'o','LineWidth',1.5);
        hold on
        errorbar(linearEa,linearErrs,'o','LineWidth',1.5);
        xticklabels(durs2)
        ylim([10 110])
        xlim([0 8])
        ylabel('Activation energy (kJ)')
        xlabel({'Stage interval';'^*times not greater than 0 were ignored'})
        legend({'Overall','Linear'})
        set(gca,'LineWidth',1.5,'FontSize',16)
        switch type
            case 1
            	title('Compare E_as (Accumulated)')
                accLinearEa = linearEa;
                accLinearErrs = linearErrs;
        	case 2
            	title('Compare E_as (Stagewise)')
                stgLinearEa = linearEa;
                stgLinearErrs = linearErrs;
        end
    %}
end

figure('Position',[10 10 600 450])
        errorbar(stgLinearEa,stgLinearErrs,'bo','LineWidth',2);
        hold on
        %errorbar(accLinearEa,accLinearErrs,'bo','LineWidth',1.5);
        xticklabels(durs2)
        xtickangle(45)
        ylim([10 110])
        xlim([0 8])
        ylabel('Activation energy (kJ)')
        xlabel({'Developmental interval';'^*times less than 0 were ignored'})
        %legend({'Stagewise','Accumulated'})
        title('Apparent activation energies');
        set(gca,'LineWidth',2,'FontSize',22) 
    %draw braces to show significance        
        y1 = stgLinearEa(1)+stgLinearErrs(1)+1;
    	drawbrace([1,y1], [3,y1], 10,'LineWidth',2.5,'Color','k');
         	text(mean([1 3]),y1+8,'\fontsize{25}***','HorizontalAlignment','center');
            
        y1 = stgLinearEa(4)+stgLinearErrs(4)+1;
    	drawbrace([3,y1], [4,y1], 10,'LineWidth',2.5,'Color','k');
         	text(mean([3 4]),y1+8,'\fontsize{25}**','HorizontalAlignment','center');
                   
        
%%
%Calculate pvals for significance
figure('Position',[100 100 1150 900]);

numStages = size(dataStg);
numStages = numStages(2)-2;

pValM = nan(numStages);
pwrVals = nan(numStages);
pValMStr = cell(numStages);

durs2 = replace(durs, '_2','-');

for x = 1:numStages
    for y = 1:numStages
        if x == y
            pValM(y,x) = 1;
                %stages compared against same stage should = 1
            pwrVals(y,x) = 0;
            txt = sprintf('%.1e',1);
            pValMStr{y,x} = replace(txt,'+0','+');
        elseif x < y
            compStages = {char(stgDurs(x)),char(stgDurs(y))};
            compRegs = stagRegT(compStages,:);
        %run stage regression table through home-made ancova function
            pValM(y,x) = real(ancova(compRegs));
            pwrVals(y,x) = calcPower(compRegs);
            txt = sprintf('%.1e',pValM(y,x));
            if pValM(y,x) < 1
                pValMStr{y,x} = replace(txt,'-0','-');
            elseif pValM(y,x) >= 1
                pValMStr{y,x} = replace(txt,'+0','+');
            end
        else
            pValM(y,x) = NaN;
            pwrVals(y,x) = NaN;
            pValMStr{y,x} = '';
        end
    end  
end  

%plot 4 plots, one of p-values, then color coded for thhreshold sig at
%0.05,0.01,0.001
%pvals
    %subplot(2,2,1)
    %pValMStr = sprintf('%.1e',num2str(pValM))
    h1 = heatmap(durs2,...
        durs2,...
        round(pValM,4,'significant'),'CellLabelColor','none','CellLabelFormat','%.1e');

    h1.Title = {'\fontsize{44}Kuntz and Eisen ANCOVA p-values'};
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
        
    set(gca,'FontSize', 40);
    xlabel('\fontsize{44}\bf{Interval 1}')
    ylabel({'\fontsize{44}\bf{Interval 2}'})
    
    
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
        
        prevPos = get(cbax,'Position');
        prevPos = [prevPos(1)-0.01 prevPos(2) ...
            prevPos(3) prevPos(4)];
        set(cbax,'Position', prevPos)
 
        
    IPx = h1.InnerPosition(1);
    IPy = h1.InnerPosition(2);
    textAxis = axes('Position',[0 0 1 1],'visible','off');
    scaleX = h1.InnerPosition(3)/numStages;
    scaleY = h1.InnerPosition(4)/numStages;
    
for x = 1:numStages
    for y = x:numStages
        num = pValM(y,x);
        if num <= 1e-9
            num = 0e0;
        else
        end
        txt = sprintf('%0.0E',num);
        if txt == '0E+00'
            txt = replace(txt,'+0','+'); 
        elseif str2num(txt) < 1 
            txt = replace(txt,'-0','-');
        elseif str2num(txt) > 1 || str2num(txt) == 1 
            txt = replace(txt,'+0','+');
        end
        xShift = IPx - 0.5*scaleX;
        yShift = IPy + h1.InnerPosition(4) + 0.5*scaleY;
        text(xShift + scaleX*x,yShift - scaleY*y,txt,'Color','k','FontSize',40,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
end    
hold on
plot([IPx h1.InnerPosition(3)+IPx],[IPy IPy],'k','LineWidth',4)
plot([IPx IPx],[IPy h1.InnerPosition(4)+IPy],'k','LineWidth',4)
plot([IPx-0.01 h1.InnerPosition(3)+IPx],[h1.InnerPosition(4)+IPy h1.InnerPosition(4)+IPy],'w','LineWidth',4)
ylim([0 1])
xlim([0 1])

%%
%BIC calculations
numStgs = length(abvs);

BICMat = nan(numStgs);
CVMat = nan(numStgs);

%build AIC table
AICTable = table('Size',[1,4],'VariableTypes',{'double',...
                    'cell','cell','cell'});
AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};

for startStage = 1:numStgs
    for endStage = startStage:numStgs
    %calculate interval times by subtracting start from endstage
        durData = dataAbs{:,endStage+2} - dataAbs{:,startStage+2};
    %build locial array to mark fittable embryos
        fittedEmsIndx = ~isnan(durData) & durData > 0;
            %if something is nan or 0 thats "wrong"

    %initialize arrays to hold all fittable points X,Y, and means        
    	allX = 1./(dataAbs{fittedEmsIndx,2}+273.15);
    	allY = log(1./durData(fittedEmsIndx)); 

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
                
          	%points as replicates
          	
            	polyNMdl1 = polyfit(allX, allY, 1);
             	polyNMdl2 = polyfit(allX, allY, 2);
                
           	BICComps = fitBICCompCalc(AICTable,polyNMdl1,polyNMdl2,1);
   
        %Store BIC comp values and Eas (just in case)
            preDec = split(num2str(abs(log(BICComps))),'.');
            BICMat(endStage,startStage) = log(BICComps);
        %if length(preDec{1}) >= 2
         %   BICMat(endStage,startStage) = round(log(BICComps),0);
        %else
         %   BICMat(endStage,startStage) = round(log(BICComps),1);
        %end
        %calculate CVs for each temperature
        uniqueTemps = unique(dataAbs.temps);
        CVPerTemp = nan(1,length(uniqueTemps));
        for temp = 1:length(uniqueTemps)
            tempIndx = dataAbs.temps == uniqueTemps(temp);
            meanDataPerTemp = mean(durData(tempIndx&fittedEmsIndx));
            stdDataPerTemp = std(durData(tempIndx&fittedEmsIndx));
            CVPerTemp(temp) = stdDataPerTemp./meanDataPerTemp;
        end
            CVMat(endStage,startStage) = mean(CVPerTemp,'omitnan')*100;
        else
        end
    end
end

%plot the relevent BIC data range and the associated stage abbreviations
    fh = figure('position',[100 100 1000 800]);
    ylabels = abvs(2:numStgs);
    xlabels = abvs(1:numStgs-1);
    %exp space
    h = heatmap(xlabels, ylabels, round(BICMat(2:numStgs,1:numStgs-1),2,'significant'),'CellLabelColor','none');

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
        
    set(gca,'FontSize', 40);
    xlabel('\fontsize{44}\bf{Start score code}')
    ylabel({'\fontsize{44}\bf{End score code}'})
    
    %determine temperatures used and adjust title accordingly
        annotation('textbox',[0.960, 0.932, 0.00005, 0.00005], 'string', '+','FontSize',40,'LineWidth',0.01)
        annotation('textbox',[0.975, 0.172, 0.00005, 0.00005], 'string', '-','FontSize',40,'LineWidth',0.01)
            %reduce textbox size to near zero to get rid of it
        title('\fontsize{44}Kuntz and Eisen Fly ln(L_{Q}/L_{L})');
    
    
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

%{      
%plot the relevent CVs data range and the associated stage abbreviations
    fh = figure('position',[100 100 1000 800]);
    ylabels = abvs(2:numStgs);
    xlabels = abvs(1:numStgs-1);
    %exp space
        h = heatmap(xlabels, ylabels, round(CVMat(2:numStgs,1:numStgs-1),2,'significant'),'CellLabelColor','black');

%Aestetics
    
    caxis([0 30])
    h.MissingDataColor = [1 1 1];
    h.GridVisible = 'off';
    h.ColorbarVisible = 'on';
    %colormap
        topx = flip(0.0:0.001:1);
        topmap = horzcat(ones(length(topx),1),(1*(topx)).',(1*(topx)).');
        botx = flip(0.0:0.001:1);
        botmap = horzcat((1*flip(botx)).',ones(length(botx),1),(1*flip(botx)).');
        map = vertcat(botmap,topmap);
        colormap(h,map)
   
    set(gca,'FontSize', 32);
    xlabel('\fontsize{38}\bf{Start score code}')
    ylabel({'\fontsize{38}\bf{End score code}'})
    
    %determine temperatures used and adjust title accordingly
        annotation('textbox',[0.950, 0.942, 0.00005, 0.00005], 'string', '+','FontSize',32,'LineWidth',0.01)
        annotation('textbox',[0.964, 0.1435, 0.00005, 0.00005], 'string', '-','FontSize',32,'LineWidth',0.01)
            %reduce textbox size to near zero to get rid of it
        title('\fontsize{40}Eisen Fly CVs');

%}        
    IPx = h.InnerPosition(1);
    IPy = h.InnerPosition(2);
    textAxis = axes('Position',[0 0 1 1],'visible','off');
    scaleX = h.InnerPosition(3)/(numStgs-1);
    scaleY = h.InnerPosition(4)/(numStgs-1);
%
    for x = 1:numStgs-1
        for y = x:numStgs-1
            txt = num2str(round(BICMat(y+1,x),2,'significant'));
            if str2num(txt) < 10 && contains(txt,'.') == 0
                txt = strcat(txt,'.0');
            else
            end
            xShift = IPx - 0.5*scaleX;
            yShift = IPy + h.InnerPosition(4) + 0.5*scaleY;
            text(xShift + scaleX*x,yShift - scaleY*y,txt,'Color','k','FontSize',h.FontSize,'HorizontalAlignment','center','VerticalAlignment','middle')                
        end
    end

    
%} 
%%
%recreating Eisen fig 3 but with fits
figure('Position',[100 100 900 1000])
%comparison image
subplot(2,1,1)
I = imread('/Users/josephcrapse/Desktop/ArrheniusPaper/Revisions Round 2/EisenComparison/draftFigures/eisenFig3BC.png');
imshow(I)

subplot(2,1,2)
numStgs = length(abvs);
color = [21, 40, 209; 33, 154, 252; 0, 158, 42; 250, 246, 2; 255, 193, 23;...
    196, 22, 6; 194, 27, 250; 12, 0, 92]/255;

%initialize table to calc pvals
    stagRegM = zeros(length(abvs),8);
    stagRegT = array2table(stagRegM);
    stagRegT.Properties.VariableNames = {'numPoints','SX','SY',...
        'SXX','SYY','SXY','AveX','AveY'};
    stagRegT.Properties.RowNames = abvs;


for stage = 1:numStgs
    proportion = dataAbs{:,stage+2}./dataAbs.tfTime ;
    indx = ~isnan(proportion);
    indxLinear = ~isnan(proportion) & dataAbs.temps <= 27.5;
    fittedX = dataAbs.temps(indx);
    fittedY = proportion(indx);
    scatter(proportion(indx),dataAbs.temps(indx),[],color(stage,:),'s','filled')
        hold on
    %fit linearity in transposed space, then plot in same space as scatter
    lm = fit(dataAbs.temps(indx), proportion(indx),'poly1');
    lmEval = polyval([lm.p1 lm.p2],[15 35]);
    ov = plot(lmEval, [15 35], 'Color',color(stage,:),'LineWidth',2);
    
    lm = fit(dataAbs.temps(indxLinear), proportion(indxLinear),'poly1');
    lmEval = polyval([lm.p1 lm.p2],[15 35]);
    lin = plot(lmEval,[15 35],'--', 'Color',color(stage,:),'LineWidth',2);
    
    %table to calculate pvals later
        stagRegT.numPoints(stage) = sum(indx);
        stagRegT.SX(stage) = sum(fittedX);
        stagRegT.SY(stage) = sum(fittedY);
        stagRegT.SXX(stage) = sum(fittedX.^2);
        stagRegT.SYY(stage) = sum(fittedY.^2);
        stagRegT.SXY(stage) = sum(fittedX.*fittedY);
        stagRegT.AveX(stage) = mean(fittedX);
        stagRegT.AveY(stage) = mean(fittedY); 
    
end
grid on
grid minor
legend([ov lin],{'All Temps','<= 27.5 ^oC'})
title('Replotted data from Figure 3C from Kuntz and Eisen')
ylabel('Temperature ^oC')
xlabel('Proportion of Development')
set(gca,'LineWidth',2,'FontSize',14)


%Calculate pvals for significance
figure('Position',[100 100 1100 900]);

pValM = nan(numStgs);
pwrVals = nan(numStgs);

for x = 1:numStgs
    for y = x:numStgs
        if x == y
            pValM(y,x) = 1;
                %stages compared against same stage should = 1
            pwrVals(y,x) = 0;
        elseif x == 2 && y == numStgs
            pValM(y,x) = 1;
                %stages compared against same stage should = 1
            pwrVals(y,x) = 0;
        else
            compStages = {char(abvs(x)),char(abvs(y))};
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
    h1 = heatmap(abvs,...
        abvs,...
        round(pValM,4),'CellLabelColor','black','CellLabelFormat','%.1e');

    h1.Title = {'\fontsize{40}Eisen Fly ANCOVA p-values',...
        '\fontsize{40}developmental stage proportions'};
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
    xlabel('\fontsize{38}\bf{Score 1}')
    ylabel({'\fontsize{38}\bf{Score 2}'})
    
    
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