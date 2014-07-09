% The single week version shows only the results of a single week when the
% grid position remains the same. In addition, firing rates are also
% displayed.

% The RF version differs from the ripple version in that here we never show
% the energy plots. The LFP and Spike data are shown in separate plots.

% The rafiki version also shows the data from all electrodes. The layout is
% also slightly different. Effectively, we merge showRFSingleGridPosAbu.m
% with electrodeSelectionOneWeek.m, both in programs/abu

% Modified from showRFSingleGridPosRafiki

function displayReceptiveFields(monkeyName,expDates,protocolNames,folderSourceString,gridType,axisLims,aziCenters,eleCenters)

numberOfDays = length(expDates);
timePeriodPos=1; % 40 to 100 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts 
fontSizeSmall=10; fontSizeMedium=12; fontSizeLarge=14;
backgroundColor = 'w';
%axisLims = [2 7 -3.5 3.5];

% Colors
posColors = jet(numberOfDays); 
posColors(numberOfDays+1,:) = [0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots

hRMSPlot = subplot('Position',[0.05 0.75 0.125 0.14],'XTickLabel',[],'YTickLabel',[]);
hFRPlot = subplot('Position',[0.05 0.6 0.125 0.14],'XTickLabel',[],'YTickLabel',[]);
hGridPlot = subplot('Position',[0.05 0.45 0.125 0.14],'XTickLabel',[],'YTickLabel',[]);
hRFPlotLFPAllElectrodes  = subplot('Position',[0.05 0.05 0.425 0.395]);

hRFPlotLFP     = subplot('Position',[0.225  0.8 0.125 0.1]);
hMaxValPlotLFP = subplot('Position',[0.225 0.71 0.125 0.06]);
hMeansLFP      = subplot('Position',[0.225 0.45 0.125 0.24],'XTickLabel',[],'YTickLabel',[],'box','on');

hRFPlotFR   = subplot('Position',[0.35 0.8 0.125 0.1],'YTickLabel',[]);%,'XTickLabel',[]);
hMaxValPlotFR  = subplot('Position',[0.35 0.71 0.125 0.06],'YAxisLocation','right');
hMeansFR    = subplot('Position',[0.35 0.45 0.125 0.24],'XTickLabel',[],'YTickLabel',[],'box','on');

[~,aValsUnique,eValsUnique] = loadRFParams(monkeyName,expDates{1},protocolNames{1},folderSourceString,gridType,'LFP');
numAzi=length(aValsUnique);
numEle=length(eValsUnique);
%numAzi=11; numEle=11;
numRows=numEle;numCols=numAzi;
hERPPlots = getPlotHandles(numRows,numCols,[0.525 0.5 0.425 0.4],0,0);
hFRPlots   = getPlotHandles(numRows,numCols,[0.525 0.05 0.425 0.4],0,0);

hERPAndFRPlots{1} = hERPPlots;
hERPAndFRPlots{2} = hFRPlots;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[electrodeList,electrodeString]=getElectrodeList(monkeyName,expDates{1},protocolNames{1},folderSourceString,gridType);
hElectrodeNum = uicontrol('Unit','Normalized', 'Position',[0.05 0.9 0.075 0.1], ...
    'Style','popup','String',electrodeString,'FontSize',fontSizeSmall);

dayString='';
for e=1:length(expDates)
    dayString = cat(2,dayString,[num2str(e) ' ' expDates{e} protocolNames{e} '|']);
end
dayString = [dayString 'allDays'];
    
hDayNum = uicontrol('Unit','Normalized', 'Position',[0.05 0.875 0.075 0.1], ...
    'Style','popup','String',dayString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingTextWidth=0.25; timingBoxWidth=0.125; timingHeight=0.5;
hStimulusCenterPanel = uipanel('Title','Stimulus Center','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.15 0.925 0.25 0.075]);

centerRange = [-3.75 -1.75];
uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'Position',[0 timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Center','FontSize',fontSizeMedium);

hAzi = uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(centerRange(1)),'FontSize',fontSizeSmall);
hEle = uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(centerRange(2)),'FontSize',fontSizeSmall);

uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'Position',[0.5 0.5 0.25 0.5], ...
    'Style','text','String','Radius','FontSize',fontSizeMedium);

hRadius = uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(0.1),'FontSize',fontSizeSmall);

uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'Position',[0 0 0.5 0.5], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotStimulus_Callback});

uicontrol('Parent',hStimulusCenterPanel,'Unit','Normalized', ...
    'Position',[0.5 0 0.5 0.5], ...
    'Style','pushbutton','String','plot all','FontSize',fontSizeMedium, ...
    'Callback',{@plotAllStimuli_Callback});

%%%%%%%%%%%%%%%%%%%%%% All Electrodes Plotting option %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlottingPanelAllElectrodes = uipanel('Title','All electrodes','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.4 0.925 0.2 0.075]);

uicontrol('Parent',hPlottingPanelAllElectrodes,'Unit','Normalized', ...
    'Position',[0 0.5 0.5 0.5], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotAllChannels_Callback});

uicontrol('Parent',hPlottingPanelAllElectrodes,'Unit','Normalized', ...
    'Position',[0.5 0.5 0.5 0.5], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium,'Callback',{@holdOnAllElectrodes_Callback});

uicontrol('Parent',hPlottingPanelAllElectrodes,'Unit','Normalized', ...
    'Position',[0.5 0 0.5 0.5], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@claAllElectrodes_Callback});

uicontrol('Parent',hPlottingPanelAllElectrodes,'Unit','Normalized', ...
    'Position',[0 0 0.5 0.5], ...
    'Style','pushbutton','String','rescale','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleAllElectrodes_Callback});

%%%%%%%%%%%%%%%%%%%%%% Plotting option %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlottingPanel = uipanel('Title','Single electrode','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.8 0.925 0.2 0.075]);

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', ...
    'Position',[0 0.5 0.5 0.5], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotSingleChannel_Callback});

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', ...
    'Position',[0.5 0.5 0.5 0.5], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium,'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', ...
    'Position',[0.5 0 0.5 0.5], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', ...
    'Position',[0 0 0.5 0.5], ...
    'Style','pushbutton','String','rescale','FontSize',fontSizeMedium, ...
    'Callback',{@rescale_Callback});

%%%%%%%%%%%%%%%%%%%%%% Time, frequency and cLims %%%%%%%%%%%%%%%%%%%%%%%%%%
timingTextWidth=0.25; timingBoxWidth=0.125; timingHeight=0.5;
hTimingPanel = uipanel('Title','Ranges','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.6 0.925 0.2 0.075]);

% Time
timeRange = [-0.1 0.3];
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Time range','FontSize',fontSizeMedium);

hTimeMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(timeRange(1)),'FontSize',fontSizeSmall);
hTimeMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(timeRange(2)),'FontSize',fontSizeSmall);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 0 0.25 0.5], ...
    'Style','text','String','Thres','FontSize',fontSizeMedium);

hRMSThreshold = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.25 0 0.25 0.5], ...
    'Style','edit','String',num2str(60),'FontSize',fontSizeSmall);

% Frequency
% frequencyRange = [1 250];
% uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
%     'Position',[0 0 timingTextWidth timingHeight], ...
%     'Style','text','String','Frequency range','FontSize',fontSizeMedium);
% 
% hFreqMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
%     'BackgroundColor', backgroundColor, ...
%     'Position',[timingTextWidth 0 timingBoxWidth timingHeight], ...
%     'Style','edit','String',num2str(frequencyRange(1)),'FontSize',fontSizeSmall);
% hFreqMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
%     'BackgroundColor', backgroundColor, ...
%     'Position',[timingTextWidth+timingBoxWidth 0 timingBoxWidth timingHeight], ...
%     'Style','edit','String',num2str(frequencyRange(2)),'FontSize',fontSizeSmall);

% CLims for Energy
rmsRange = [0 200];
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0.5 timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','ERP','FontSize',fontSizeMedium);

hRMSMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.5+timingTextWidth timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(rmsRange(1)),'FontSize',fontSizeSmall);
hRMSMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.5+timingTextWidth+timingBoxWidth timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(rmsRange(2)),'FontSize',fontSizeSmall);

% dEnergy
frRange = [0 200];
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0.5 0 timingTextWidth timingHeight], ...
    'Style','text','String','Firing (s/s)','FontSize',fontSizeMedium);

hFRMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.5+timingTextWidth 0 timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(frRange(1)),'FontSize',fontSizeSmall);
hFRMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.5+timingTextWidth+timingBoxWidth 0 timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(frRange(2)),'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plotStimulus_Callback(~,~)
        
        paramsStimulus(1) = str2double(get(hAzi,'string'));
        paramsStimulus(2) = str2double(get(hEle,'string'));
        paramsStimulus(3) = str2double(get(hRadius,'string'));
        paramsStimulus(4) = paramsStimulus(3);
        paramsStimulus(5)=0;
        paramsStimulus(6)=1;
        
        [~,~,boundaryXStimulus,boundaryYStimulus] = gauss2D(paramsStimulus);
        hold(hRFPlotLFPAllElectrodes,'on');
        plot(hRFPlotLFPAllElectrodes,boundaryXStimulus,boundaryYStimulus,'k');
    end

    function plotAllStimuli_Callback(~,~)
        
        for i=1:length(aziCenters)
            
            paramsStimulus(1) = aziCenters(i);
            paramsStimulus(2) = eleCenters(i);
            paramsStimulus(3) = str2double(get(hRadius,'string'));
            paramsStimulus(4) = paramsStimulus(3);
            paramsStimulus(5)=0;
            paramsStimulus(6)=1;
            
            [~,~,boundaryXStimulus,boundaryYStimulus] = gauss2D(paramsStimulus);
            hold(hRFPlotLFPAllElectrodes,'on');
            plot(hRFPlotLFPAllElectrodes,boundaryXStimulus,boundaryYStimulus,'k');
        end
    end

    function plotAllChannels_Callback(~,~)
        
        dayNum = get(hDayNum,'val');
        rmsThreshold = str2double(get(hRMSThreshold,'string'));
        
        cla(hRMSPlot); cla(hFRPlot);cla(hGridPlot);
            
        if dayNum<=numberOfDays
            useTheseExpDates      = expDates(dayNum);
            useTheseProtocolNames = protocolNames(dayNum);
        else
            useTheseExpDates      = expDates;
            useTheseProtocolNames = protocolNames;
        end

        analyzeAllElectrodes(hRMSPlot,hFRPlot,hGridPlot,hRFPlotLFPAllElectrodes,monkeyName,useTheseExpDates,useTheseProtocolNames,folderSourceString,gridType,timePeriodPos,rmsThreshold,axisLims);
        cRMS = caxis(hRMSPlot); set(hRMSMin,'String',num2str(cRMS(1))); set(hRMSMax,'String',num2str(cRMS(2)));
        cFR = caxis(hFRPlot); set(hFRMin,'String',num2str(cFR(1))); set(hFRMax,'String',num2str(cFR(2)));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotSingleChannel_Callback(~,~)
        
        electrodeNum = electrodeList(get(hElectrodeNum,'val'));
        dayNum = get(hDayNum,'val');
        timeLims = [str2double(get(hTimeMin,'String')) str2double(get(hTimeMax,'String'))];

        cla(hRFPlotLFP); cla(hRFPlotFR);
        cla(hMeansLFP); cla(hMeansFR);
        cla(hMaxValPlotLFP); cla(hMaxValPlotFR);
        
        set(hRFPlotLFP,'NextPlot','add');
        set(hRFPlotFR,'NextPlot','add');
        
        for i=1:numberOfDays
            expDate = expDates{i};
            protocolName = protocolNames{i};
            
            % Get LFP Data
            allParamsLFP = loadRFParams(monkeyName,expDate,protocolName,folderSourceString,gridType,'LFP');
            paramsLFP = allParamsLFP{electrodeNum,timePeriodPos};
            [~,~,boundaryXLFP,boundaryYLFP] = gauss2D(paramsLFP);
            
            % Azi and Ele list
            aziListLFP(i) = paramsLFP(1); 
            eleListLFP(i) = paramsLFP(2);
            sigXLFP(i) = paramsLFP(3);
            sigYLFP(i) = paramsLFP(4);
            scalingFactorLFP(i) = paramsLFP(6);
            plot(hRFPlotLFP,aziListLFP(i),eleListLFP(i),'color',posColors(i,:),'Marker','*');
            plot(hRFPlotLFP,boundaryXLFP,boundaryYLFP,'color',posColors(i,:));
            
            [~,~,rfValsMaxLFP] = loadRMSAndMaxValues(monkeyName,expDate,protocolName,folderSourceString,gridType,'LFP');
            maxValLFP(i) = max(max(squeeze(rfValsMaxLFP(:,:,electrodeNum,timePeriodPos))));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Get Spike Data
            allParamsFR = loadRFParams(monkeyName,expDate,protocolName,folderSourceString,gridType,'Spikes');
            paramsFR = allParamsFR{electrodeNum,timePeriodPos};
            [~,~,boundaryXFR,boundaryYFR] = gauss2D(paramsFR);
            
            % Azi and Ele list
            aziListFR(i) = paramsFR(1); 
            eleListFR(i) = paramsFR(2);
            sigXFR(i) = paramsFR(3);
            sigYFR(i) = paramsFR(4);
            scalingFactorFR(i) = paramsFR(6);
            plot(hRFPlotFR,aziListFR(i),eleListFR(i),'color',posColors(i,:),'Marker','*');
            plot(hRFPlotFR,boundaryXFR,boundaryYFR,'color',posColors(i,:));
            
            % Get Spike Data
            [~,~,rfValsMaxFR] = loadRMSAndMaxValues(monkeyName,expDate,protocolName,folderSourceString,gridType,'Spikes');
            maxValFR(i) = max(max(squeeze(rfValsMaxFR(:,:,electrodeNum,timePeriodPos))));
        end
        
        axis(hRFPlotLFP,axisLims); axis(hRFPlotFR,axisLims);
        showMeans(hMeansLFP,aziListLFP,eleListLFP,sigXLFP,sigYLFP,posColors,scalingFactorLFP);
        showMeans(hMeansFR,aziListFR,eleListFR,sigXFR,sigYFR,posColors,scalingFactorFR);
        
        % Plot the mean ERP and firing rates
        hold(hMaxValPlotLFP,'on'); hold(hMaxValPlotFR,'on');
        for i=1:numberOfDays
            plot(hMaxValPlotLFP,i,maxValLFP(i),'color',posColors(i,:),'marker','o');
            plot(hMaxValPlotFR,i,maxValFR(i),'color',posColors(i,:),'marker','o');
        end
        axis(hMaxValPlotLFP,[0 numberOfDays+1 0 max(maxValLFP)]);
        axis(hMaxValPlotFR,[0 numberOfDays+1 0 max(max(maxValFR),1)]);
        ylabel(hMaxValPlotLFP,'\muV');
        ylabel(hMaxValPlotFR,'Spikes/s');
        
        %disp(maxValFR)
        % Show ERP and FRs        
        if dayNum<=numberOfDays
            expDatesToUse = expDates(dayNum);
            protocolNamesToUse = protocolNames(dayNum);
            maxValLFPToUse = maxValLFP(dayNum);
            maxValFRToUse = maxValFR(dayNum);
        else
            expDatesToUse = expDates;
            protocolNamesToUse = protocolNames;
            maxValLFPToUse = maxValLFP;
            maxValFRToUse = maxValFR;
        end
        
        normalizeData=0;
        if normalizeData
            maxValLFPToUse = ones(1,length(maxValLFPToUse));
            maxValFRToUse = ones(1,length(maxValFRToUse));
        end
        plotERPAndFRDataNew(hERPAndFRPlots,monkeyName,expDatesToUse,protocolNamesToUse,folderSourceString,gridType,electrodeNum,numRows,numCols,posColors(dayNum,:),timeLims,maxValLFPToUse,maxValFRToUse);
    end
    function cla_Callback(~,~)
        claGivenPlotHandle(hRFPlotLFP);claGivenPlotHandle(hRFPlotFR);
        claGivenPlotHandle(hMaxValPlotLFP);claGivenPlotHandle(hMaxValPlotFR);
        claGivenPlotHandle(hERPAndFRPlots{1});
        claGivenPlotHandle(hERPAndFRPlots{2});
    end
    function claAllElectrodes_Callback(~,~)
        cla(hRMSPlot);
        cla(hFRPlot);
        cla(hGridPlot);
        cla(hRFPlotLFPAllElectrodes);
    end
    function rescale_Callback(~,~)
        timeLims = [str2double(get(hTimeMin,'String')) str2double(get(hTimeMax,'String'))];
        %freqLims = [str2double(get(hFreqMin,'String')) str2double(get(hFreqMax,'String'))];
        
        rescaleAxes(hERPAndFRPlots{1},[timeLims getYLims(hERPAndFRPlots{1})]); changeTickLabels(hERPAndFRPlots{1});
        rescaleAxes(hERPAndFRPlots{2},[timeLims getYLims(hERPAndFRPlots{2})]); changeTickLabels(hERPAndFRPlots{2});
    end
    function rescaleAllElectrodes_Callback(~,~)
        rmsLims = [str2double(get(hRMSMin,'String')) str2double(get(hRMSMax,'String'))];
        frLims = [str2double(get(hFRMin,'String')) str2double(get(hFRMax,'String'))];
        
        caxis(hRMSPlot,rmsLims); caxis(hFRPlot,frLims);
    end
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
        holdOnGivenPlotHandle(hERPAndFRPlots{1},holdOnState);
        holdOnGivenPlotHandle(hERPAndFRPlots{2},holdOnState);
    end
    function holdOnAllElectrodes_Callback(source,~)
        holdOnState = get(source,'Value');
        holdOnGivenPlotHandle(hRFPlotLFPAllElectrodes,holdOnState);
    end
end
% Loading functions
function [numStimuli,rfValsRMS,rfValsMax,rfValsPowerOrMean] = loadRMSAndMaxValues(monkeyName,expDate,protocolName,folderSourceString,gridType,measure)
load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\RFMeasures\' measure '\rfValues.mat']);

if strcmp(measure,'LFP')
    rfValsPowerOrMean = rfValsPower;
elseif strcmp(measure,'Spikes')
    rfValsPowerOrMean = rfValsMean;
end

end
function [params,aValsUnique,eValsUnique] = loadRFParams(monkeyName,expDate,protocolName,folderSourceString,gridType,measure) %#ok<*STOUT>
poolingOption=2;
load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\RFMeasures\' measure '\rfParams' num2str(poolingOption) '.mat']); % poolingOption

if strcmp(measure,'LFP')
    params = paramsRMSScaled; %paramsMaxScaled,paramsPowerScaled,paramsRMS,paramsMax,paramsPower
elseif strcmp(measure,'Spikes')
    params = paramsRMSScaled;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rescaling
function holdOnGivenPlotHandle(plotHandles,holdOnState)

[numRows0,numCols0] = size(plotHandles);
if holdOnState
    for ii=1:numRows0
        for jj=1:numCols0
            set(plotHandles(ii,jj),'Nextplot','add');
        end
    end
else
    for ii=1:numRows0
        for jj=1:numCols0
            set(plotHandles(ii,jj),'Nextplot','replace');
        end
    end
end
end
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims = [yMin yMax];
end
function rescaleAxes(plotHandles,xyLims)

[numRows,numCols] = size(plotHandles);

for row=1:numRows
    for column=1:numCols 
        axis(plotHandles(row,column),xyLims);
    end
end
end
function changeTickLabels(plotHandles)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:numRows
    for j=1:numCols
        if (i==numRows && rem(j,2)==1)
            if j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        else
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end
function claGivenPlotHandle(plotHandles)
[numRows0,numCols0] = size(plotHandles);
for ii=1:numRows0
    for jj=1:numCols0
        cla(plotHandles(ii,jj));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showMeans(hMeans,aziList,eleList,sigx,sigy,posColors,scalingList)
cla(hMeans);

aziList = round(100*aziList)/100;
eleList = round(100*eleList)/100;
sigXList = sigx;
sigYList = sigy;
rfSizeList = round(100*sqrt((sigXList.^2 + sigYList.^2)/2))/100;

N = length(aziList);
text('Parent',hMeans,'String','Azi','Position',[0.2 1-1/(N+3)], 'HorizontalAlignment','center','Color','k');
text('Parent',hMeans,'String','Ele','Position',[0.5 1-1/(N+3)], 'HorizontalAlignment','center','Color','k');
text('Parent',hMeans,'String','RF' ,'Position',[0.8 1-1/(N+3)], 'HorizontalAlignment','center','Color','k');

for i=1:N
    text('Parent',hMeans,'String',num2str(aziList(i)),'Position',[0.2 1-(i+1)/(N+3)],'HorizontalAlignment','center','Color',posColors(i,:));
    text('Parent',hMeans,'String',num2str(eleList(i)),'Position',[0.5 1-(i+1)/(N+3)],'HorizontalAlignment','center','Color',posColors(i,:));
    text('Parent',hMeans,'String',num2str(rfSizeList(i)),'Position',[0.8 1-(i+1)/(N+3)],'HorizontalAlignment','center','Color',posColors(i,:));
end

% Remove bad indices
badFits = find(scalingList<0);
badRFs =  [find(rfSizeList<0.1) find(rfSizeList>=1)];

badIndices = unique([badFits badRFs]);

for i=1:length(badIndices)
    text('Parent',hMeans,'String','X','Position',[0.95 1-(badIndices(i)+1)/(N+3)],'HorizontalAlignment','center','Color','r');
end

aziList(badIndices)=[];
eleList(badIndices)=[];
rfSizeList(badIndices)=[];

grandMeanAzi = round(100*mean(aziList))/100;
grandMeanEle = round(100*mean(eleList))/100;

text('Parent',hMeans,'String',num2str(grandMeanAzi),'Position',[0.2 1-(N+2.5)/(N+3)],'HorizontalAlignment','center','Color','k');
text('Parent',hMeans,'String',num2str(grandMeanEle),'Position',[0.5 1-(N+2.5)/(N+3)],'HorizontalAlignment','center','Color','k');
text('Parent',hMeans,'String',num2str(round(100*mean(rfSizeList))/100),'Position',[0.8 1-(N+2.5)/(N+3)],'HorizontalAlignment','center','Color','k');
end
function plotERPAndFRDataNew(hERPAndFRPlots,monkeyName,expDatesToUse,protocolNamesToUse,folderSourceString,gridType,electrodeNum,numRows,numCols,plotColor,timeLims,maxValLFP,maxValFR)

numDays = length(expDatesToUse);
spikeRateLimit=20;
spikeColor = plotColor; lfpColor = plotColor;

hERPPlots = hERPAndFRPlots{1};
hFRPlots  = hERPAndFRPlots{2};

% Initialization
load([folderSourceString 'data\' monkeyName '\' gridType '\' expDatesToUse{1} '\' protocolNamesToUse{1} '\extractedData\parameterCombinations']);
load([folderSourceString 'data\' monkeyName '\' gridType '\' expDatesToUse{1} '\' protocolNamesToUse{1} '\segmentedData\LFP\lfpInfo']);

aValsUniqueFirstDay = aValsUnique;
eValsUniqueFirstDay = eValsUnique;

if length(eValsUnique) ~= numRows
    error('numRows not equal to the number of elevation entries');
end

if length(aValsUnique) ~= numCols
    error('numCols not equal to the number of azimuth entries');
end
    
numGoodFRDays=1;
for i=1:numDays
    disp(['working on day ' num2str(i) ' of ' num2str(numDays)]);
    expDate = expDatesToUse{i};
    protocolName = protocolNamesToUse{i};
    
    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\extractedData\parameterCombinations']);
    %%%%%%%%%% check if the entries are consistent  %%%%%%%%%%%%%%%%%%%%%%%
    if aValsUnique ~= aValsUniqueFirstDay
        error(['azimuth entries are different for ' expDate protocolName  ' and ' expDatesToUse{1} protocolNamesToUse{1}]);
    end
    
    if eValsUnique ~= eValsUniqueFirstDay
        error(['elevation entries are different for ' expDate protocolName  ' and ' expDatesToUse{1} protocolNamesToUse{1}]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\RFMeasures\LFP\meanLFPDataChan' num2str(electrodeNum)]);
    numStimuliLFP=numStimuli; clear numStimuli
    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\RFMeasures\Spikes\meanSpikeDataChan' num2str(electrodeNum) '_SID0']);
    numStimuliFR=numStimuli;
    % LFP data
    for e=1:numRows
        for a=1:numCols
            normalizedLFPData(i,e,a,:) = meanLFPData(e,a,:)/maxValLFP(i); %#ok<*NODEF>
            allNumStimuliLFP(i,e,a) = numStimuliLFP(e,a);
        end
    end
    
    if maxValFR(i) > spikeRateLimit
        for e=1:numRows
            for a=1:numCols
                normalizedSpikeData(numGoodFRDays,e,a,:) = meanSpikeData(e,a,:)/maxValFR(i);
                allNumStimuliFR(numGoodFRDays,e,a) = numStimuliFR(e,a);
            end
        end
        numGoodFRDays=numGoodFRDays+1;
    end
end

% Mean LFP Data
meanNormalizedLFPData = squeeze(mean(normalizedLFPData,1));
totalStimsLFP =   squeeze(sum(allNumStimuliLFP,1));
disp('LFP data: ')
disp(totalStimsLFP);

if numGoodFRDays>1
    disp(['Max spikeRate > ' num2str(spikeRateLimit) ' on ' num2str(numGoodFRDays-1) ' days']);
    meanNormalizedSpikeData = squeeze(mean(normalizedSpikeData,1));
    totalStimsFR =   squeeze(sum(allNumStimuliFR,1));
    disp(totalStimsFR);
    showFR=1;
else
    disp(['None of the days showed max spiking activity > ' num2str(spikeRateLimit)]);
    showFR=0;
end

for a=1:numCols
    colPos = a;    
    for e=1:numRows
        rowPos = numRows-e+1;
        plot(hERPPlots(rowPos,colPos),timeVals,squeeze(meanNormalizedLFPData(e,a,:)),'color',lfpColor);
        
        if showFR
            plot(hFRPlots(rowPos,colPos),timeValsFR,squeeze(meanNormalizedSpikeData(e,a,:)),'color',spikeColor);
        end
        
        if colPos==1
            text('Parent',hERPPlots(rowPos,colPos),'Position',[0.1 0.1],'Units','Normalized', ...
                    'String', num2str(eValsUniqueFirstDay(e)));
            text('Parent',hFRPlots(rowPos,colPos),'Position',[0.1 0.1],'Units','Normalized', ...
                    'String', num2str(eValsUniqueFirstDay(e)));
        end
        
        if rowPos==1
            text('Parent',hERPPlots(rowPos,colPos),'Position',[0.5 1.1],'Units','Normalized', ...
                    'String', num2str(aValsUniqueFirstDay(a)));
            text('Parent',hFRPlots(rowPos,colPos),'Position',[0.5 1.1],'Units','Normalized', ...
                    'String', num2str(aValsUniqueFirstDay(a)));    
        end
    end
end

rescaleAxes(hERPPlots,[timeLims getYLims(hERPPlots)]); changeTickLabels(hERPPlots);
rescaleAxes(hFRPlots,[timeLims getYLims(hFRPlots)]); changeTickLabels(hFRPlots);
end
function analyzeAllElectrodes(hRMSPlot,hFRPlot,hGridPlot,hRFPlotLFPAllElectrodes,monkeyName,expDates,protocolNames,folderSourceString,gridType,timePeriodPos,rmsThreshold,axisLims)

numDays = length(expDates);
badElectrodes=[];

electrodeList=getElectrodeList(monkeyName,expDates{1},protocolNames{1},folderSourceString,gridType);
numElectrodes = length(electrodeList);
for i=1:numDays    
    % Impedances
    impedanceFile = [folderSourceString 'data\' monkeyName '\' gridType '\' expDates{i} '\impedanceValues.mat'];
    if exist(impedanceFile,'file')
        load(impedanceFile);
        badElectrodes = unique([badElectrodes find(impedanceValues>2500)]);
    else
        disp([impedanceFile ' does not exist! Bad Electrodes not found..']);
    end

    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDates{i} '\' protocolNames{i} '\RFMeasures\LFP\rfValues.mat']);
    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDates{i} '\' protocolNames{i} '\RFMeasures\LFP\rfParams2.mat']);
    
    rfVals(i,:,:,:,:) = rfValsRMS; %#ok<*AGROW>
    params = paramsRMSScaled;

    for jj=1:numElectrodes
        j=electrodeList(jj);
        azi(i,j) = params{j,timePeriodPos}(1);
        ele(i,j) = params{j,timePeriodPos}(2);
%         rfSizeAzi(i,j) = params{j,timePeriodPos}(3);
%         rfSizeEle(i,j) = params{j,timePeriodPos}(4);
%         rfSize(i,j) = sqrt((params{j,timePeriodPos}(3)^2 + params{j,timePeriodPos}(4)^2)/2);
%         scalingFactor(i,j) = params{j,timePeriodPos}(6);
%         theta(i,j)    = params{j,timePeriodPos}(5);
    end
    
    % Get the same for spikes
    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDates{i} '\' protocolNames{i} '\RFMeasures\Spikes\rfValues.mat']);
    load([folderSourceString 'data\' monkeyName '\' gridType '\' expDates{i} '\' protocolNames{i} '\RFMeasures\Spikes\rfParams2.mat']);
    rfValsSpikes(i,:,:,:,:) = rfValsRMS;
%    paramsSpikes = paramsRMSScaled;
% 
%     for jj=1:numElectrodes
%         j=electrodeList(jj);
%         aziSpikes(i,j) = paramsSpikes{j,timePeriodPos}(1);
%         eleSpikes(i,j) = paramsSpikes{j,timePeriodPos}(2);
%         rfSizeAziSpikes(i,j) = paramsSpikes{j,timePeriodPos}(3);
%         rfSizeEleSpikes(i,j) = paramsSpikes{j,timePeriodPos}(4);
%         rfSizeSpikes(i,j) = sqrt((paramsSpikes{j,timePeriodPos}(3)^2 + paramsSpikes{j,timePeriodPos}(4)^2)/2);
%         scalingFactorSpikes(i,j) = paramsSpikes{j,timePeriodPos}(6);
%         
%         thetaSpikes(i,j)    = paramsSpikes{j,timePeriodPos}(5);
%     end
end

meanRFVals = squeeze(mean(rfVals,1));
meanRFValsSpikes = squeeze(mean(rfValsSpikes,1));

% numChans=size(meanRFVals,3);
% numTimePeriods = size(meanRFVals,4);

for ii=1:numElectrodes
    i=electrodeList(ii);
    [row,column] = electrodePositionOnGrid(i,gridType);
    maxVals(i) = max(max(squeeze(meanRFVals(:,:,i,timePeriodPos))));
    maxValsByElectrode(row,column) = maxVals(i);
    
    maxValsSpikes(i) = max(max(squeeze(meanRFValsSpikes(:,:,i,timePeriodPos))));
    maxValsByElectrodeSpikes(row,column) = maxValsSpikes(i);
end

% Plot
axes(hRMSPlot);imagesc(maxValsByElectrode); set(hRMSPlot,'XTickLabel',[],'YTickLabel',[]); ylabel(hRMSPlot,'LFP RMS'); %#ok<*MAXES>
axes(hFRPlot);imagesc(maxValsByElectrodeSpikes); set(hFRPlot,'XTickLabel',[],'YTickLabel',[]); ylabel(hFRPlot,'Firing rate');
goodElectrodes=setdiff(find(maxVals>rmsThreshold),badElectrodes);

% Plot grid in color
[~,electrodeColorNames,electrodeArray] = showElectrodeLocationsInColor([],hGridPlot,1,0,1,gridType);
showElectrodeLocations([],setdiff(electrodeList,goodElectrodes),'w',hGridPlot,1,1,gridType);
showElectrodeLocations([],badElectrodes,'k',hGridPlot,1,1,gridType);

% Receptive fields
hold(hRFPlotLFPAllElectrodes,'on');
axes(hRFPlotLFPAllElectrodes);
del=0.02;
for i=1:length(goodElectrodes)
    
    thisAzi = azi(:,goodElectrodes(i));
    thisEle = ele(:,goodElectrodes(i));
    
    [row,col] = find(goodElectrodes(i) == electrodeArray);
    colorList = electrodeColorNames{row,col};
    
    for j=1:numDays
        plot(hRFPlotLFPAllElectrodes,thisAzi(j),thisEle(j),'color',colorList,'Marker','+','markersize',12);     
    end
    
    if numDays==1
        if thisAzi<axisLims(1) || thisAzi > axisLims(2) || thisEle<axisLims(3) || thisEle > axisLims(4)
            disp(['Electrode ' num2str(goodElectrodes(i)) ', center: (' num2str(thisAzi) ',' num2str(thisEle) ') is out of range']);
        else
            text(azi(:,goodElectrodes(i))+del,ele(:,goodElectrodes(i)),num2str(goodElectrodes(i)));
        end
    end
end
hold(hRFPlotLFPAllElectrodes,'off');
axis(hRFPlotLFPAllElectrodes,axisLims);
end
function [electrodeList,electrodeString]=getElectrodeList(monkeyName,expDate,protocolName,folderSourceString,gridType)

folderName = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];
folderSegment = [folderName 'segmentedData\'];
load([folderSegment 'LFP\lfpInfo']);
electrodeList = analogChannelsStored;

electrodeString='';
for i=1:length(electrodeList);
    electrodeString = cat(2,electrodeString,['elec' num2str(electrodeList(i)) '|']);
end
electrodeString=electrodeString(1:end-1);
end