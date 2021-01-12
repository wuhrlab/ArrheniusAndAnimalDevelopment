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
%Global Search

%Run for various network sizes, initialize those sizes here
    netSizeS = 2;
%table for each overall run
    %columns = bestSolution, worstSolution, runTime
    %will have to parse solutions to make fit in table
    data = table('Size',[length(netSizeS) 11],'VariableTypes',...
        ["double","cell","cell","double", "cell","cell","cell","double",...
        "cell", "double", "double"]);
    data.Properties.VariableNames = {'RxnNum','BestSolX_Ea',...
        'BestSolX_A','BestSolFval','BestSolX0', 'WorstSolX_Ea',...
        'WorstSolX_A','WorstSolFval','WorstSolX0', 'RunTime','NumUniqueSols'};
    data.Properties.RowNames = cellstr(num2str(netSizeS));
%initialize constants   
    R = 8.3144598*(10^(-3));

%run global solver for each network size
for i = 1:length(netSizeS)
%reset problem conditions
    %number of attempts to find global min and evals
        numMsAttempts = 100;
        maxEvaluation = 225000;
        
%set up optimization solver
    %set Ea and A bound, constraint?
        netSize = netSizeS(i);
        Eaub = repmat(100,netSize,1);
        Aub = repmat(50,netSize,1);  
        Ealb = repmat(20,netSize,1);
        Alb = repmat(-10,netSize,1);
            %log space so can search uniformly, transform once in
            %base e space

%set initial conditions, do worst predicted based on math
        worstEa = 20 - log(netSize-1)*(R*295.15);
        x0Ea = [100; repmat(worstEa,netSize-1,1)];
        x0A = [37.2926197009076; repmat(4.79899926729709,netSize-1,1)];
  
%Set problem options and define the problem to be optimized 
        opts = optimoptions(@fmincon,'Algorithm','SQP',...
            'OptimalityTolerance',10^-17,'MaxFunctionEvaluations',...
            maxEvaluation);
    %Optimize curvature with calcCurvature
         problem = createOptimProblem('fmincon',...
            'objective',@(x) calcCurvature(x(:,1),x(:,2)),...
            'x0',[x0Ea,x0A],...
            'lb',[Ealb,Alb],...
            'ub',[Eaub,Aub],...
            'nonlcon',@(x) kConstraint(x(:,1),x(:,2)),...
            'options',opts);
        
%reset ms
    ms = MultiStart;

%run actual problem
    %returns the Fval (minIntegralEvaluation) and X (associated Eas and As) ...
        %...in a GlobalOptimSolution structure
    startTic = tic;
        %time the program run time for stats
    [x,f,exitflag,output,solutions] = run(ms, problem, numMsAttempts);
    tm = toc(startTic);
    
%Store statistics and results for each network size run to compare later
    data.RunTime(i) = tm;
    data.RxnNum(i) = netSize;
    data.BestSolX_Ea(i) = {solutions(1,1).X(:,1)};
        %Ea column 1, ln(A) column 2
    data.BestSolX_A(i) = {solutions(1,1).X(:,2)};
    data.BestSolFval(i) = solutions(1,1).Fval;
    data.BestSolX0(i) = {solutions(1,1).X0};
    lS = length(solutions);
    data.WorstSolX_Ea(i) = {solutions(1,lS).X(:,1)};
    data.WorstSolX_A(i) = {solutions(1,lS).X(:,2)};
    data.WorstSolFval(i) = solutions(1,lS).Fval;
    data.WorstSolX0(i) = {solutions(1,lS).X0};
    data.NumUniqueSols(i) = lS;
end

%%
%Response to reviewers
%Theoretical most nonlinear, with noise in A
%initialize experimental variables and function
    %subplot(2,2,3);
    R = 8.3144598*(10^(-3));
    netSize = 1000;

    worstEa = 20 - log(netSize-1)*(R*295.15);
    Ea = [100; repmat(worstEa,netSize-1,1)];
    lnA = [37.2926197009076; repmat(4.79899926729709,netSize-1,1)];

    %Set temperature range and initialize temps
        rangeT = [280.15 307.65];
        resolutionT = 30;
        temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';

    %plot "mean"
        %write function with NO noise
        	fcn1 = writelnKFcnVectorized(Ea,lnA);
        %calculate lnKs using temp range, resolution, and function(Ea,lnA)
            lnKs = double(fcn1(temps));
            
            magenta = plot(1./temps,lnKs,'m--','LineWidth',2.5);
            hold on
    %two kinds of noise, noise applied to ever A for ever embryo
    %or noise applied uniformly to all A within an embryo but diff
    %between diff embryos
    numRuns = 100;
    lnKsNoise = NaN(numRuns,1);
    for runs = 1:numRuns
        for i = 1:length(temps)
            %calculate random noise, gassuian distributed, in non-log space
                noise = normrnd(0,0.1*exp(lnA));
            %write function with noise
                fcn1 = writelnKFcnVectorized(Ea,log(exp(lnA)+noise));
            %calculate lnKs using temp range, resolution, and function(Ea,lnA)
                lnKsNoise(runs,i) = double(fcn1(temps(i)));
        end
    end
    %plot data with noise
        error = std(lnKsNoise,1).';
        blue = errorbar(1./temps,mean(lnKsNoise,1).',error,error,'bo',...
            'LineWidth',2.5);
    
    %Aestetics
    %control position so can focus on curvature
        xlim([min(1./temps) max(1./temps)])
        xticks(1./([30 20 10]+273.15))
        
    %Conditional titles and legends
        title({'Worst Case + Noise (In each Rxn)'},'FontSize',25)

    %calculate regressions
            tsK = [287.45 289.25 290.35 291.75 293.65 294.25 295.95 ...
                297.45 298.05 300.15];
            rtsK = [282.55 303.25 306.55];
            errTemp = 0.5; 
            linFit = polyfit(1./tsK',double(fcn1(tsK))',1);
                      
            xneg = 1./tsK-1./(tsK+errTemp);
            data = errorbar(1./tsK,double(fcn1(tsK)),xneg,'horizontal',...
                'oc','LineWidth',2.5);
            %data.Color = [127,200,175]/255;
            rxneg = 1./rtsK-1./(rtsK+errTemp);
            rdata = errorbar(1./rtsK,double(fcn1(rtsK)),rxneg,...
                'horizontal','or','LineWidth',2.5);
            
        %plot resultant lnKs from function(Ea,lnA)      
            fitEval = polyval(linFit,1./[280.15 307.15]);
            blue = plot(1./[280.15 307.15],fitEval,'b','LineWidth',2.5);
            ylim([min(fitEval)-.25 max(fitEval)])
            
            magenta = plot(1./temps,lnKs.','--r','LineWidth',2.5);

        %legend
            ylabs = cell2mat(cellfun(@str2num,yticklabels,'un',0));
            yticks(ylabs(ylabs == floor(ylabs)));
            yticklabels(ylabs(ylabs == floor(ylabs)));
            xlabel('1/(T [^oC] + 273.15)','FontSize',28)
            xticklabels(strcat('1/',strsplit(num2str(round((1./...
                xticks-273.15),1)))));
        %convert decimal xticks to 1/T in K
            set(gca,'linewidth', 3, 'FontSize', 28);
            ylabel('ln(k)','FontSize',28)
            legend([blue magenta], {'100x Simulated Data',...
                'Underlying Model'},'FontSize',24);            
 
%%   
%Figure 4B
%Complexity
%Theoretical most nonlinear plotted with random
%initialize experimental variables and function
    figure('Position',[50 400 800 600])

    R = 8.3144598*(10^(-3));
    netSize = 1000;
    
    %Set temperature range and initialize temps
    	rangeT = [280.15 307.65];
    	resolutionT = 60;
    	temps = (rangeT(1):range(rangeT)/(resolutionT-1):rangeT(2)).';
        above =  temps<(273.15+14);
        below =  temps<(273.15+30);
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

                tangT = find(round(temps,1)==295.1);
                tang = calcTangent(1./temps, worstlnKs, tangT);    
                black = plot(1./temps,tang,'k','LineWidth',2.5);
                    hold on
                %{
                linFitWorst = polyfit(1./fitTs',double(fcnWorst(fitTs))',1);
                fitEval = polyval(linFitWorst,1./[280.15 307.15]);
                black = plot(1./[280.15 307.15],fitEval,'k','LineWidth',2.5);
                    hold on
%}
                ylabel('ln(k) worst case','FontSize',28)
                ylim([min(tang) max(tang)]); 
                ylabs = cell2mat(cellfun(@str2num,yticklabels,'un',0));
                yticks(ylabs(ylabs == floor(ylabs)));
                yticklabels(ylabs(ylabs == floor(ylabs)));

    %random        
        Ea = (100-20).*rand(1000,1) + 20;
        lnKset = (0+12.3884).*rand(1000,1) - 12.3884;
            %control the ks to be "reasonable" and better control the A
        lnA = lnKset + Ea/(R*295.5);        
        %plot data with noise
    %plot "mean"
            %write function with NO noise
                fcnRnd = writelnKFcnVectorized(Ea,lnA);
            %calculate lnKs using temp range, resolution, and function(Ea,lnA)
                rndlnKs = double(fcnRnd(temps));

                yyaxis left 
                cyan = plot(1./temps,rndlnKs,'--','LineWidth',2.5);
                    hold on
            
                tangT = find(round(temps,1)==295.1);
                tang = calcTangent(1./temps, rndlnKs, tangT);    
                black = plot(1./temps,tang,'k','LineWidth',2.5);
                    hold on
%}
                 %{  
                linFitRnd = polyfit(1./fitTs',double(fcnRnd(fitTs))',1);
                fitEval = polyval(linFitRnd,1./[280.15 307.15]);
                black = plot(1./[280.15 307.15],fitEval,'k','LineWidth',2.5);
                    hold on
%}
                ylabel('ln(k) random network','FontSize',28)
                ylim([min(tang) max(tang)])

        %Aestetics
        %control position so can focus on curvature
            xticks(1./([30 20 10]+273.15))

        %Conditional titles and legends
            %title({'1000x Random Network vs Worst Case'},'FontSize',25)

        %legend
            ylabs = cell2mat(cellfun(@str2num,yticklabels,'un',0));
            yticks(ylabs(ylabs == floor(ylabs)));
            yticklabels(ylabs(ylabs == floor(ylabs)));
            xlabel('1/(T [^oC] + 273.15)','FontSize',28)
            xlim(1./[307.15 280.15])
            xticklabels(strcat('1/',strsplit(num2str(round((1./...
                xticks-273.15),1)))));
        %convert decimal xticks to 1/T in K
            set(gca,'linewidth', 3, 'FontSize', 28);
            legend([cyan magenta black], {'Random Model',...
                'Worst Case Model','Tangent at 295 ^oK'},'FontSize',24);            
        
%%
%Figure App2 A
    fh = figure('Position',[50 50 600 500]);
 
%Frog data import
    %Data Import
    %for all following analysis
    %Destination
        %filename = 'FrogScoresFinalCombinedCleaned.xlsx';
        filename = 'FrogScoresFinalCombinedCleanedGastRescore.xlsx';

        folder = '/Volumes/MyAppleFriend/Arrhenius Paper Work/Independent Work/Frog Scores/';

    %read in all sheet names
        [~,sheet_name]=xlsfinfo(strcat(folder,filename));

    %read in example sheet to build data variable arrays from 
        sizeData = length(sheet_name);
        data = readcell(strcat(folder,filename),'Sheet', sheet_name{1});

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


%calculate concavity in frog
%as well as collect error, temps, numSamples

        obsScores = (2:12);
        accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
        stageAbbsNoT0 = stageAbbs(obsScores);  
        stageNamesNoT0 = stageNames(obsScores);    

    %Initialize fit comparison Table
        AICTable = table('Size',[1,4],'VariableTypes',{'double',...
            'cell','cell','cell'});
        AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};

    %initialize matrix to hold errors for monte carlo
        monteTimeErrorsExtremes = cell(length(accumulatedInterval{1}));
        intervalAllXExtremes = cell(length(accumulatedInterval{1}));

    %initialize fit comparison arrays
        AICValues = zeros(numScores,4);

    %initialize Data for Arrhenius space x and y
        lnPseudoKperStage = cell(length(accumulatedInterval),1);
        invertedTemps = 1./Temps;

    %initialize presentation data matrecies    
        sampleSizeMat = cell(length(accumulatedInterval{1}));
        uniqueTemps = cell(length(accumulatedInterval{1}));
        isConcaveBioFrog = nan(numScores);
        isConvexBioFrog = nan(numScores);
        isLinearBioFrog = nan(numScores);

    %control what temperatures to exclude from the linear fit BIC comparison
    %calculation
        counter = 1;

    %calculate poly n 1 and 2 fits for every developmental interval then compare 
    %them using BIC
    excludedTemps = [76];
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
            fitTemps = zeros(length(lnPseudoKperStage),1);
            fittabeEmNum = 0;
            fitEmbs = zeros(length(lnPseudoKperStage),1);
            for temp = 1:length(lnPseudoKperStage)
            %if not one of the excluded temperature and more then one embryo
                if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                    if sum(~isnan(lnPseudoKperStage{temp})) > 1
                        fitTemps(temp) = 1;
                        fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}));
                        fitEmbs(temp) = sum(~isnan(lnPseudoKperStage{temp}));
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
                    yStd = std(data,'omitnan'); 
                    tempararyError(temp) = yStd;
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
                monteTimeErrorsExtremes{endStage,startStage} = tempararyError.'; 
                sampleSizeMat{endStage,startStage} = fitEmbs;
                %collect true/false if I should use a temp from Temps
                uniqueTemps{endStage,startStage} = fitTemps;
                
            %Build BIC comparision variables, comparing models, then run
            %through BIC program
            if startStage~=endStage
            %store data points for AIC
            %all points as replicates

            polyNMdl2 = polyfit(allX, allY, 2);

            %Check concavity        
                isConcaveBioFrog(endStage,startStage) = polyNMdl2(1) > 0; 
                isConvexBioFrog(endStage,startStage) = polyNMdl2(1) < 0; 
                isLinearBioFrog(endStage,startStage) = polyNMdl2(1) == 0; 
                
            else
            end
        end
    end

%MonteCarlo, frog
    %compute one underlying model for Temps
        R = 8.3144598*(10^(-3));
        netSize = 1000;
    
	%Set temperature range and initialize temps
        temps = Temps;

    %worst case
        worstEa = 20 - log(netSize-1)*(R*295.15);
        Ea = [100; repmat(worstEa,netSize-1,1)];
        lnA = [37.2926197009076; repmat(4.79899926729709,netSize-1,1)];

	%plot "mean"
        %write function with NO noise
            fcnWorst = writelnKFcnVectorized(Ea,lnA);
        %calculate lnKs using temp range, resolution, and function(Ea,lnA)
            worstlnKs = double(fcnWorst(temps));

    numRuns = 100;   
    isConcaveMonteFrog = zeros(numScores);
    isConvexMonteFrog = zeros(numScores);
    isLinearMonteFrog = zeros(numScores);

    for startStage = 1:(length(monteTimeErrorsExtremes)-1)
        for endStage = (startStage+1):length(monteTimeErrorsExtremes)
            for run = 1:numRuns
                allY = nan(sum(sampleSizeMat{endStage,startStage},'omitnan'),1);
                allX = nan(sum(sampleSizeMat{endStage,startStage},'omitnan'),1);

                %introduce error into Temp and worstlnKs
                dataPos = 1;
                for temp = find(uniqueTemps{endStage,startStage}).'
                    numEmbs = sampleSizeMat{endStage,startStage}(temp);               
                    %repeat worst case ln(k) x number of embryos at this temp
                    %introduce error in equal to error for this interval at
                    %this temp
                    tempY = repmat(worstlnKs(temp),numEmbs,1) + normrnd(0,...
                        monteTimeErrorsExtremes{endStage,startStage}(temp),...
                        numEmbs,1);

                    tempX = repmat(Temps(temp),numEmbs,1) + normrnd(0,0.5,numEmbs,1);

                    allY(dataPos:dataPos+numEmbs-1) = tempY;
                    allX(dataPos:dataPos+numEmbs-1) = 1./tempX;

                    dataPos = dataPos + numEmbs;
                end

                polyNMdl2 = polyfit(allX, allY, 2);            

                isConcaveMonteFrog(endStage,startStage) = isConcaveMonteFrog(...
                    endStage,startStage) + double(polyNMdl2(1) > 0); 
                isConvexMonteFrog(endStage,startStage) = isConvexMonteFrog(...
                    endStage,startStage) + double(polyNMdl2(1) < 0); 
                isLinearMonteFrog(endStage,startStage) = isLinearMonteFrog(...
                    endStage,startStage) + double(polyNMdl2(1) == 0); 
            end
        end
    end     
    
    
%Data Import
%import fly data
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

        
%calculate concavity in fly
%as well as collect error, temps, numSamples

        obsScores = (1:12);
        accumulatedInterval = cellfun(@(x) x(:,obsScores),accumulated,'UniformOutput',false);
        stageAbbsNoT0 = stageAbbs(obsScores);  
        stageNamesNoT0 = stageNames(obsScores);    

    %Initialize fit comparison Table
        AICTable = table('Size',[1,4],'VariableTypes',{'double',...
            'cell','cell','cell'});
        AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};

    %initialize matrix to hold errors for monte carlo
        monteTimeErrorsExtremes = cell(numScores);
        intervalAllXExtremes = cell(numScores);

    %initialize fit comparison arrays
        AICValues = zeros(numScores,4);

    %initialize Data for Arrhenius space x and y
        lnPseudoKperStage = cell(length(accumulatedInterval),1);
        invertedTemps = 1./Temps;

    %initialize presentation data matrecies    
        sampleSizeMat = cell(numScores);
        uniqueTemps = cell(numScores);
        isConcaveBioFly = nan(numScores);
        isConvexBioFly = nan(numScores);
        isLinearBioFly = nan(numScores);

    %control what temperatures to exclude from the linear fit BIC comparison
    %calculation
        counter = 1;

    %calculate poly n 1 and 2 fits for every developmental interval then compare 
    %them using BIC
    excludedTemps = [76];
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
            fitTemps = zeros(length(lnPseudoKperStage),1);
            fittabeEmNum = 0;
            fitEmbs = zeros(length(lnPseudoKperStage),1);
            for temp = 1:length(lnPseudoKperStage)
            %if not one of the excluded temperature and more then one embryo
                if sum(temp ~= excludedTemps)/length(excludedTemps) == 1
                    if sum(~isnan(lnPseudoKperStage{temp})) > 1
                        fitTemps(temp) = 1;
                        fittabeEmNum = fittabeEmNum + ...
                        sum(~isnan(lnPseudoKperStage{temp}));
                        fitEmbs(temp) = sum(~isnan(lnPseudoKperStage{temp}));
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
                    yStd = std(data,'omitnan'); 
                    tempararyError(temp) = yStd;
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
                monteTimeErrorsExtremes{endStage,startStage} = tempararyError.'; 
                sampleSizeMat{endStage,startStage} = fitEmbs;
                %collect true/false if I should use a temp from Temps
                uniqueTemps{endStage,startStage} = fitTemps;
                
            %Build BIC comparision variables, comparing models, then run
            %through BIC program
            if startStage~=endStage
            %store data points for AIC
            %all points as replicates

            polyNMdl2 = polyfit(allX, allY, 2);

            %Check concavity        
                isConcaveBioFly(endStage,startStage) = polyNMdl2(1) > 0; 
                isConvexBioFly(endStage,startStage) = polyNMdl2(1) < 0; 
                isLinearBioFly(endStage,startStage) = polyNMdl2(1) == 0; 
                
            else
            end
        end
    end

%MonteCarlo, fly
%worst network
    %compute one underlying model for Temps
        R = 8.3144598*(10^(-3));
        netSize = 1000;
    
	%Set temperature range and initialize temps
        temps = Temps;

    %worst case
        worstEa = 20 - log(netSize-1)*(R*295.15);
        Ea = [100; repmat(worstEa,netSize-1,1)];
        lnA = [37.2926197009076; repmat(4.79899926729709,netSize-1,1)];

	%plot "mean"
        %write function with NO noise
            fcnWorst = writelnKFcnVectorized(Ea,lnA);
        %calculate lnKs using temp range, resolution, and function(Ea,lnA)
            worstlnKs = double(fcnWorst(temps));
            
%random network
    %random        
        Ea = (100-20).*rand(1000,1) + 20;
        lnKset = (0+12.3884).*rand(1000,1) - 12.3884;
            %control the ks to be "reasonable" and better control the A
        lnA = lnKset + Ea/(R*295.5);        
        %plot data with noise
    %plot "mean"
            %write function with NO noise
                fcnRnd = writelnKFcnVectorized(Ea,lnA);
            %calculate lnKs using temp range, resolution, and function(Ea,lnA)
                rndlnKs = double(fcnRnd(temps));
            

    numRuns = 100;   
    isConcaveMonteFly = zeros(numScores);
    isConvexMonteFly = zeros(numScores);
    isLinearMonteFly = zeros(numScores);
%Initialize worst case fit comparison Table
    BICMatWorstCase = nan(numScores);
	AICTable = table('Size',[1,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
	AICTable.Properties.VariableNames = {'numPoints','yi','xi','mdl'};

%Initialize worst case fit comparison Table
    BICMatRndCase = nan(numScores);
	AICTableRnd = table('Size',[1,4],'VariableTypes',{'double',...
    	'cell','cell','cell'});
	AICTableRnd.Properties.VariableNames = {'numPoints','yi','xi','mdl'};    
for startStage = 1:(length(monteTimeErrorsExtremes)-1)
    for endStage = (startStage+1):length(monteTimeErrorsExtremes)
        BICComps = nan(numRuns,1);
        BICCompsRnd = nan(numRuns,1);
        for run = 1:numRuns
            allY = nan(sum(sampleSizeMat{endStage,startStage},'omitnan'),1);
            allX = nan(sum(sampleSizeMat{endStage,startStage},'omitnan'),1);
            
            %introduce error into Temp and worstlnKs
            dataPos = 1;
            for temp = find(uniqueTemps{endStage,startStage}).'
                numEmbs = sampleSizeMat{endStage,startStage}(temp);               
                %repeat worst case ln(k) x number of embryos at this temp
                %introduce error in equal to error for this interval at
                %this temp
                %worst case network
                tempY = repmat(worstlnKs(temp),numEmbs,1) + normrnd(0,...
                    monteTimeErrorsExtremes{endStage,startStage}(temp),...
                    numEmbs,1);
                
                tempX = repmat(Temps(temp),numEmbs,1) + normrnd(0,0.5,numEmbs,1);

                allY(dataPos:dataPos+numEmbs-1) = tempY;
                allX(dataPos:dataPos+numEmbs-1) = 1./tempX;
                
                %random network
                tempYRnd = repmat(rndlnKs(temp),numEmbs,1) + normrnd(0,...
                    monteTimeErrorsExtremes{endStage,startStage}(temp),...
                    numEmbs,1);
                
                tempXRnd = repmat(Temps(temp),numEmbs,1) + normrnd(0,0.5,numEmbs,1);
                
                allXRnd(dataPos:dataPos+numEmbs-1) = tempYRnd;
                allYRnd(dataPos:dataPos+numEmbs-1) = 1./tempXRnd;
                
                %increment counter
                dataPos = dataPos + numEmbs;
            end         
            
            %fit and store concavity
            polyNMdl2 = polyfit(allX, allY, 2);            

            isConcaveMonteFly(endStage,startStage) = isConcaveMonteFly(...
                endStage,startStage) + double(polyNMdl2(1) > 0); 
            isConvexMonteFly(endStage,startStage) = isConvexMonteFly(...
                endStage,startStage) + double(polyNMdl2(1) < 0); 
            isLinearMonteFly(endStage,startStage) = isLinearMonteFly(...
                endStage,startStage) + double(polyNMdl2(1) == 0); 
            
            %extra code for C, BIC, worst case net
                AICTable.numPoints = length(allY);
                AICTable.yi = num2cell(allY,1);
                AICTable.xi = num2cell(allX,1);
            
                polyNMdl1 = polyfit(allX, allY, 1);

                BICComps(run) = fitBICCompCalc(AICTable,polyNMdl1,polyNMdl2,1);  
                
            %extra code for B, BIC, random net
                AICTableRnd.numPoints = length(allY);
                AICTableRnd.yi = num2cell(allY,1);
                AICTableRnd.xi = num2cell(allX,1);
            
                polyNMdl1Rnd = polyfit(allXRnd, allYRnd, 1);
                polyNMdl2Rnd = polyfit(allXRnd, allYRnd, 2);
            
                BICCompsRnd(run) = fitBICCompCalc(AICTableRnd,polyNMdl1Rnd,polyNMdl2Rnd,1);
        end
        %copy over mean values
        BICMatWorstCase(endStage,startStage) = round(log(mean(BICComps)),0);
        BICMatRndCase(endStage,startStage) = round(log(mean(BICCompsRnd)),1);
        %BICMatWorstCase(endStage,startStage) = round(median(log(BICComps)),0);
        %BICMatRndCase(endStage,startStage) = round(median(log(BICCompsRnd)),1);
    end
end
       
%A
%Concavity
%First run Monte's above and Fly/Frog BIC section that records concavity

    flyConcave = sum(isConcaveBioFly,'All','omitnan');
    flyConvex = sum(isConvexBioFly,'All','omitnan');
    flyLinear = sum(isLinearBioFly,'All','omitnan');
    
    flyConcaveMonte = sum(isConcaveMonteFly,'All','omitnan');
    flyConvexMonte = sum(isConvexMonteFly,'All','omitnan');
    flyLinearMonte = sum(isLinearMonteFly,'All','omitnan');
    
    frogConcave = sum(isConcaveBioFrog,'All','omitnan');
    frogConvex = sum(isConvexBioFrog,'All','omitnan');
    frogLinear = sum(isLinearBioFrog,'All','omitnan');
    
    frogConcaveMonte = sum(isConcaveMonteFrog,'All','omitnan');
    frogConvexMonte = sum(isConvexMonteFrog,'All','omitnan');
    frogLinearMonte = sum(isLinearMonteFrog,'All','omitnan');
    %determine total number of measurable events to normalize by
        totals = [(flyConvex + flyConcave + flyLinear) (frogConvex + frogConcave + frogLinear)...
            (flyConvexMonte + flyConcaveMonte + flyLinearMonte) (frogConvexMonte + frogConcaveMonte + frogLinearMonte)];
    %group and normalize yeses and nos
        convex = [flyConvex frogConvex flyConvexMonte frogConvexMonte]./totals;
        concave = [flyConcave frogConcave flyConcaveMonte frogConcaveMonte]./totals;
        linear = [flyLinear frogLinear flyLinearMonte frogLinearMonte]./totals;
        
        bar([convex;concave;linear])
        ylabel('Frequency')
        title({'Concavity of Developmental Intervals','Entire Temperature Range'})
        legend({'Time-lapse Fly Data','Time-lapse Frog Data','Monte Carlo w/ Fly Errors','Monte Carlo w/ Frog Errors'})
        xticklabels({'Concave \newline   Down', 'Concave \newline     Up','Linear'})
        set(gca,'FontSize', 25,'LineWidth',3);
        
%B
%mean BIC over 100x runs per interval
%random network
%plot the relevent BIC data range and the associated stage abbreviations
    fh = figure('position',[100 100 1000 800]);
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    %exp space
        h = heatmap(xlabels, ylabels, round(BICMatRndCase(2:12,1:11),2),'CellLabelColor','black');

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

        title('\fontsize{40}ln(L_{Q}/L_{L}): Random Model');    
    
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
        
        
%C
%mean BIC over 100x runs per interval
%worst case network
%plot the relevent BIC data range and the associated stage abbreviations
    fh = figure('position',[100 100 1000 800]);
    ylabels = stageAbbsNoT0(2:12);
    xlabels = stageAbbsNoT0(1:11);
    %exp space
        h = heatmap(xlabels, ylabels, round(BICMatWorstCase(2:12,1:11),2),'CellLabelColor','black');

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

        title('\fontsize{40}ln(L_{Q}/L_{L}): Worst Case Model');    
    
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
        
