%%
%Data Import
%GAPDH Enzyme Assay
%Destination
    filename = 'earlyGAPDH1to3.xlsx';
    folder = '/Users/josephcrapse/Desktop/ArrheniusPaper/Enzymatic Assays/GAPDH/';
    
    
%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);

%take temps from sheet names
    splitCells = cellfun(@(x)split(x,'c'),sheet_name,'UniformOutput',false);
    cellTemps = cellfun(@(x) x(1),splitCells);
    temps = unique(str2double(cellTemps),'sorted');
    colors = rand(length(temps),3);

%group data by temperature
    allDataByTemp = cell(1,length(temps));
    for i = 1:length(temps)
        sheets = find(str2double(cellTemps)==temps(i));
        for j = 1:length(sheets)
            data = readcell(strcat(folder,filename),'Sheet',...
                sheet_name{sheets(j)});
            allDataByTemp{i}{j} = cell2mat(data((3:length(data)),:));
        end
    end
  %}  
%initialize ks
    k = cell(1,length(allDataByTemp));
    
    minTime = 11.5;
    maxTime = 14.5;

%unnormalized
for i = 1:length(allDataByTemp)
    dataSize = length(allDataByTemp{i});
    for j = 1:dataSize
    %calculate k
        %subsect data between time 5.5 to 14
        range = floor(((allDataByTemp{i}{j}(:,1)<maxTime)+(...
            allDataByTemp{i}{j}(:,1)>minTime))/2);
        subdata = allDataByTemp{i}{j}(range==1,:);
        %linear fit
        hold on;
        pfit = polyfit(subdata(:,1), subdata(:,2),1);
        %assign slope to k
        k{i}(j) = pfit(1);
    end
end
  
%controls
%
%initialize ks
    kCont = cell(1,length(allDataByTemp));
    
    minTime = 17.5;
    maxTime = 19.5;

%unnormalized
    contCases = find(ismember(temps,[10,15,20,25,30,35,40,45,50]));
        %fix the above later so its more general
for i = contCases
    dataSize = length(allDataByTemp{i});
    for j = 1:dataSize
    %calculate k
        %subsect data between time 5.5 to 14
        range = floor(((allDataByTemp{i}{j}(:,1)<maxTime)+(...
            allDataByTemp{i}{j}(:,1)>minTime))/2);
        subdata = allDataByTemp{i}{j}(range==1,:);
        %linear fit
        hold on;
        pfit = polyfit(subdata(:,1), subdata(:,2),1);
        %assign slope to k
        kCont{i}(j) = pfit(1);
    end
end

    meanLogks = cellfun(@(x) mean(log(x)),k);
    meanLogkConts = cellfun(@(x) mean(log(x)),kCont);
%plot with errorbars for exp  
figure('Position',[100 100 600 500])
    errs = cellfun(@(x)std(log(x),1,'omitnan')./sqrt(length(x)),kCont);
    ci95 = cellfun(@(x)tinv(0.025, length(x)-1),kCont);
    %
    yStdErr = errs;
    yci95 = yStdErr;
    errsAllTemps = yci95;
%calculate 95% CI for temeprature
    blue = errorbar(1./(temps+273.15),meanLogkConts,errsAllTemps,'bo',...
        'LineWidth',1.5);
    hold on
    
%plot with errorbars for cont  
    errs = cellfun(@(x)std(log(x),1,'omitnan')./sqrt(length(x)),kCont);
    ci95 = cellfun(@(x)tinv(0.025, length(x)-1),kCont);
    
    yStdErr = errs;
    yci95 = yStdErr;
    errsAllTemps = yci95;
%calculate 95% CI for temeprature
    hold on
    
    xlabel('1/(T [^oC] + 273.15)');
    xlim(1./[325 285]);
    xticks(1./[323.15, 313.15, 303.15 293.15 283.15]);
    xticklabels({'1/50','1/40','1/30','1/20','1/10'});
    ylabel('ln(k) (rate)')
    %title('Updated: Enzyme Kinetics of GAPDH')
    set(gca,'linewidth', 2.5, 'FontSize', 20);
    
    lm = polyfit(1./(temps(2:6)+273.15),meanLogkConts(2:6),1);
    lmVal = polyval(lm,1./(temps+273.15));
    linfit = plot(1./(temps+273.15),lmVal,'b--','LineWidth',1.5);
    
    lm = polyfit(1./(temps+273.15),meanLogkConts,2);
    lmVal = polyval(lm,1./(temps+273.15));
    qfit = plot(1./(temps+273.15),lmVal,'m--','LineWidth',1.5);
    
    legend([blue, linfit, qfit],{'0.5 mM GAP, 5 mM NAD+',...
        'Linear Fit: Temps = 15-35', 'Quadratic Fit: All Temps'},...
        'Location','northeast');
    %legend([red, blue],{'Exp','Control'},'Location','northeast');
    xlabel('1/(T [^oC] + 273)');
    xlim(1./(([54 4]+273.15)));
    xticks(1./[323.15, 313.15, 303.15 293.15 283.15]);
    xticklabels({'1/50','1/40','1/30','1/20','1/10'});
    ylabel('ln(k) (rate)')
    ylabs = yticks;
    yticks(ylabs(floor(ylabs)==ylabs));
    %title('Updated: Enzyme Kinetics of GAPDH')
    set(gca,'linewidth', 2.5, 'FontSize', 20);

    
%control comparison
    figure('Position',[100 100 600 500])

    tenMM = meanLogks;
    twentyMM = meanLogkConts;
    dColors = distinguishable_colors(length(tenMM));
    plot([-7 -3],[-7 -3],'k--','LineWidth',1.5)
    hold on
    for i = 1:length(tenMM)
        dataPoints(i) = scatter(tenMM(i),twentyMM(i),[],dColors(i,:),...
            'Filled','LineWidth',3);
        txt = strcat(string(temps(i))," ^oC");
        if i < length(tenMM)-3
            text(tenMM(i),twentyMM(i),txt,'HorizontalAlignment','center','VerticalAlignment','cap','FontSize',23)
        elseif sum(i == [8 9]) == 1
        	text(tenMM(i),twentyMM(i),txt,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',23)
        else
            text(tenMM(i),twentyMM(i),txt,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',23)
        end
    end
    xlim([-7 -3])
    ylim([-7 -3])
    xlabel('ln(k) at 0.25mM GAP, 0.25mM NAD+')
    ylabel('ln(k) at 0.5mM GAP, 5mM NAD+')
    title('GAPDH')
    for i = 1:length(temps)
        %legLabs(i) = strcat(string(temps(i))," ^oC")
    end
    %legend(dataPoints,legLabs,'Location','northwest','NumColumns',2);
    set(gca,'linewidth', 3, 'FontSize', 24);    
%display controls as fold change

%%
%Data Import
%B-Gal Assay (controls done seperately)
%Destination
%non-controls   
    filename = 'BGal_Reps_All.xlsx';
    folder = '/Users/josephcrapse/Desktop/ArrheniusPaper/Enzymatic Assays/BGal/';  
    
%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);

%take temps from sheet names
    splitCells = cellfun(@(x)split(x,'c'),sheet_name,'UniformOutput',false);
    cellTemps = cellfun(@(x) x(1),splitCells);
    temps = unique(str2double(cellTemps),'sorted');
    colors = horzcat(((1:length(cellTemps))/length(cellTemps)).',repmat(...
        0,length(cellTemps),1),flip(((1:length(cellTemps))/length(cellTemps)).'));

%group non-controls by data by temperature
    allDataByTemp = cell(1,length(temps));
    for i = 1:length(temps)
        sheets = find(str2double(cellTemps)==temps(i));
        for j = 1:length(sheets)
            data = readcell(strcat(folder,filename),'Sheet',...
                sheet_name{sheets(j)});
            allDataByTemp{i}{j} = cell2mat(data((3:length(data)),:));
        end
    end
  %}  
%initialize ks
    k = cell(1,length(allDataByTemp));
    
    minTime = 10.5;
    maxTime = 11.5;

%unnormalized
for i = 1:length(allDataByTemp)
    dataSize = length(allDataByTemp{i});
    for j = 1:dataSize
    %calculate k
        %subsect data between time 5.5 to 14
        range = floor(((allDataByTemp{i}{j}(:,1)<maxTime)+(allDataByTemp{i}{j}(:,1)>minTime))/2);
        subdata = allDataByTemp{i}{j}(range==1,:);
        %linear fit
        hold on;
        pfit = polyfit(subdata(:,1), subdata(:,2),1);
        %assign slope to k
        k{i}(j) = pfit(1);
    end
end
    
    filename = 'BGal_Conts_All.xlsx';
    folder = '/Users/josephcrapse/Desktop/ArrheniusPaper/Enzymatic Assays/BGal/';   
    
%read in all sheet names
    [~,sheet_name]=xlsfinfo(strcat(folder,filename));
    
%read in example sheet to build data variable arrays from 
    sizeData = length(sheet_name);

%take temps from sheet names
    splitCells = cellfun(@(x)split(x,'c'),sheet_name,'UniformOutput',false);
    cellTemps = cellfun(@(x) x(1),splitCells);
    contTemps = unique(str2double(cellTemps),'sorted');
    colors = horzcat(((1:length(cellTemps))/length(cellTemps)).',...
        repmat(0,length(cellTemps),1),flip(((1:length(cellTemps))/...
        length(cellTemps)).'));

%group non-controls by data by temperature
    allContDataByTemp = cell(1,length(contTemps));
    for i = 1:length(contTemps)
        sheets = find(str2double(cellTemps)==contTemps(i));
        for j = 1:length(sheets)
            data = readcell(strcat(folder,filename),'Sheet',...
                sheet_name{sheets(j)});
            allContDataByTemp{i}{j} = cell2mat(data((3:length(data)),:));
        end
    end


%initialize ks
    kCont = cell(1,length(allContDataByTemp));
    
    minTime = 10.5;
    maxTime = 11.5;

%unnormalized
        %fix the above later so its more general
for i = 1:length(allContDataByTemp)
    dataSize = length(allContDataByTemp{i});
    for j = 1:dataSize
    %calculate k
        %subsect data between time 5.5 to 14
        range = floor(((allContDataByTemp{i}{j}(:,1)<maxTime)+(...
            allContDataByTemp{i}{j}(:,1)>minTime))/2);
        subdata = allContDataByTemp{i}{j}(range==1,:);
        %linear fit
        hold on;
        pfit = polyfit(subdata(:,1), subdata(:,2),1);
        %assign slope to k
        kCont{i}(j) = pfit(1);
    end
end
    
%plot k and kCont against T    
    figure('Position',[100 100 800 600])

    meanLogks = cellfun(@(x) mean(log(x)),k);
    meanLogkConts = cellfun(@(x) mean(log(x)),kCont);
    
    exp = scatter(1./(temps+273.15), meanLogks,'ro','LineWidth',2);
    hold on
    cont = scatter(1./(contTemps+273.15), meanLogkConts,'b*','LineWidth',2);
    fit = polyfit(1./(temps(2:4)+273.15), meanLogks(2:4),1);
        %fit the "interior" temperatures
    fiteVal = polyval(fit,1./(temps(1:6)+273.15));
    plot(1./(temps(1:6)+273.15),fiteVal,'k--','LineWidth',1.5);
  
    legend([exp, cont],{'Exp','Control'},'Location','northeast');
    xlabel('1/(T [^oC] + 273)');
    xlim(1./[343.15 278]);
    xticks(1./[343.15, 333.15, 323.15, 313.15, 303.15 293.15 283.15]);
    ylabel('ln(k) (rate)')
    title('Enzyme Kinetics of \beta-gal')
    set(gca,'linewidth', 2.5, 'FontSize', 20);
%
%plot with errorbars for exp  
figure('Position',[100 100 600 500])
    errs = cellfun(@(x)std(log(x),1,'omitnan')./sqrt(length(x)),k);
    ci95 = cellfun(@(x)tinv(0.025, length(x)-1),k);
    yStdErr = errs;
    yci95 = yStdErr;
    errsAllTemps = yci95;
%calculate 95% CI for temeprature
    red = errorbar(1./(temps(1:8)+273.15),meanLogks(1:8),errsAllTemps(1:8),'b+','LineWidth',2);
    hold on
    
%plot with errorbars for cont  
    hold on
    fitKData = [];
    fitTData = [];
    fitTs = [15,25,35];
    for i = 1:length(allDataByTemp)-1
        red = scatter(repmat(1./(temps(i)+273.15),1,length(k{i})),log(k{i}),'bo','LineWidth',2);
        if sum(temps(i)==fitTs)>0
            fitKData = horzcat(fitKData,log(k{i}));
            fitTData = horzcat(fitTData,repmat(1./(temps(i)+273.15),1,length(k{i})));
        end
    end
        
        fit = polyfit(fitTData, fitKData,1);
        %fit the "interior" temperatures
        fiteVal = polyval(fit,1./(temps+273.15));
        black = plot(1./(temps+273.15),fiteVal,'b--','LineWidth',1.5);
        
    fitKData = [];
    fitTData = [];
    fitTs = [5,15,25,35,45,50,55,60];
    for i = 1:length(allDataByTemp)-1
        if sum(temps(i)==fitTs)>0
            fitKData = horzcat(fitKData,log(k{i}));
            fitTData = horzcat(fitTData,repmat(1./(temps(i)+273.15),1,length(k{i})));
        end
    end
        qfit = polyfit(fitTData, fitKData,2);
        qfiteVal = polyval(qfit,1./(temps+273.15));
        magenta = plot(1./(temps+273.15),qfiteVal,'m--','LineWidth',1.5);  
        
    legend([red, black, magenta],{'10mM ONPG','Linear Fit: Temps = 15-35','Quadratic Fit: All Temps'},'Location','northeast');

    xlabel('1/(T [^oC] + 273)');

    xlim(1./[343.15 273.15]);

    xticks(1./[323.15, 303.15 283.15]);
    xticklabels({'1/50','1/30','1/10'});
    ylabel('ln(k) (rate)')
    title('Enzyme Kinetics of \beta-gal')
    set(gca,'linewidth', 3, 'FontSize', 24);
 
%control plot
    figure('Position',[100 100 600 500])

    dColors = [0 0 1; 0 1 0; 1 0 0; 0 1 1];
    tenMM = meanLogks([1 3 6 8]);
    twentyMM = meanLogkConts([1 2 4 5]);
    plot([-5 0],[-5 0],'k--','LineWidth',1.5)
    hold on
    for i = 1:length(tenMM)
        dataPoints(i) = scatter(twentyMM(i),tenMM(i),[],dColors(i,:),'filled','LineWidth',3);
        txt = strcat(string(temps(i))," ^oC");
        text(twentyMM(i),tenMM(i),txt,'HorizontalAlignment','center','VerticalAlignment','cap','FontSize',23)
    end
    xlim([-5 0])
    ylim([-5 0])
    ylabel('ln(k) at 10 mM')
    xlabel('ln(k) at 20 mM')
    title('\beta-gal Kinetics: 10 mM vs 20 mM')
    title('\beta-gal')
    counter = 1;
    for i = [1 3 6 8]
        %legLabs(counter) = strcat(string(temps(i))," ^oC");
        counter = counter + 1;
    end
    %legend(dataPoints,legLabs,'Location','northwest')
    set(gca,'linewidth', 3, 'FontSize', 24);
