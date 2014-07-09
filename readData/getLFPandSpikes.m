% This program extracts the LFP and Spikes from the .nev file and stores the data in
% the specified folder

% Inputs
% fileName: name of the .nev file
% analogElectrodesToStore: List of analog electrodes to store.
% folderIn: folder where filename is located. Default: C:\Supratim\rawData\
% folderOut: folder where the output data will be stored. Default: folderIn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program needs the following matlab files
% makeDirectory.m
% appendIfNotPresent.m
% removeIfPresent.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray,
% September 2008

% June 2014: In previous versions we only recorded data from the electrodes
% but not input/auxilliary channels. The code is modified to allow reading
% both. "analogChannel" now is sepatared into either an "electrode" or an
% "analogInput".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function getLFPandSpikes(fileName,analogElectrodesToStore,folderIn,folderOut,...
    goodStimTimes,timeStartFromBaseLine,deltaT,Fs,hFile,neuralChannelsToStore,getLFP,getSpikes)

% Initialize
if ~exist('folderIn','var')     folderIn = 'C:\Supratim\rawData\';      end
if ~exist('folderOut','var')    folderOut = folderIn;                   end
if ~exist('hFile','var')        hFile = [];                             end
if ~exist('getLFP','var')       getLFP=1;                               end
if ~exist('getSpikes','var')    getSpikes=1;                            end

% Append appropriate tags
fileName = appendIfNotPresent(fileName,'.nev');
folderIn = appendIfNotPresent(folderIn,'\');
folderOut = appendIfNotPresent(folderOut,'\');

if isempty(hFile)
    % Load the appropriate DLL
    dllName = 'N:\programs\requiredResources\nsNEVLibrary64.dll';
    [nsresult] = ns_SetLibrary(dllName);
    if (nsresult ~= 0)
        error('DLL was not found!');
    end
end

% Load data file and display some info about the file open data file

if ~isempty(hFile)
    % Get file information
    [nsresult, fileInfo] = ns_GetFileInfo(hFile);
    % Gives you entityCount, timeStampResolution and timeSpan
    if (nsresult ~= 0)
        error('Data file information did not load!');
    end
else
    disp('Getting file info...')
    [nsresult, hFile] = ns_OpenFile([folderIn fileName]);
    if (nsresult ~= 0)
        error('Data file did not open!');
    else
        disp('Done getting file info...');
    end

    % Get file information
    [nsresult, fileInfo] = ns_GetFileInfo(hFile);
    % Gives you entityCount, timeStampResolution and timeSpan
    if (nsresult ~= 0)
        error('Data file information did not load!');
    end
end

% Build catalogue of entities
[~, entityInfo] = ns_GetEntityInfo(hFile, 1:fileInfo.EntityCount);

% List of EntityIDs needed to retrieve the information and data
%eventList = find([entityInfo.EntityType] == 1);
analogList = find([entityInfo.EntityType] == 2);
segmentList = find([entityInfo.EntityType] == 3);
neuralList = find([entityInfo.EntityType] == 4);

% How many of a particular entity do we have
% cEvent = length(eventList);
cAnalog = length(analogList);
cSegment = length(segmentList);
cNeural = length(neuralList);

% From the entityInfo.Labels, identify the channels that have been recorded
% eventLabels = strvcat(entityInfo(eventList).EntityLabel);
analogLabels = char(entityInfo(analogList).EntityLabel);
segmentLabels = char(entityInfo(segmentList).EntityLabel);
neuralLabels = char(entityInfo(neuralList).EntityLabel);

%%%%%%% Separate AnalogChannels into Electrode and AnalogInput data %%%%%%%%%%%%%%%%
electrodeCount = 0;
ainpCount = 0;

for i=1:cAnalog
    if strcmp(analogLabels(i,1:4),'chan')
        electrodeCount = electrodeCount+1;
        electrodeNums(electrodeCount) = str2num(analogLabels(i,5:end)); %#ok<*AGROW,*ST2NM>
        electrodeListIDs(electrodeCount,:) = analogList(i);
        
    elseif strcmp(analogLabels(i,1:4),'ainp')
        ainpCount = ainpCount+1;
        analogInputNums(ainpCount) = str2num(analogLabels(i,5:end));
        analogInputListIDs(ainpCount) = analogList(i);
    end
end

segmentChannelNums = str2num(segmentLabels(:,5:end)); %#ok<ST2NM>
neuralChannelNums  = str2num(neuralLabels(:,5:end)); %#ok<ST2NM>

% Display these numbers
disp(['Total number of Analog channels recorded: ' num2str(cAnalog) ', electrodes: ' num2str(electrodeCount) ', Inp: ' num2str(ainpCount)]);
disp(['Total number of Segments recorded: ' num2str(cSegment)]);
disp(['Total number of Neurons recorded: ' num2str(cNeural)]);

% Get the desired number of samples
numberOfItems = [entityInfo.ItemCount];
analysisOnsetTimes = goodStimTimes + timeStartFromBaseLine;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% LFP Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if getLFP && (cAnalog>0)
    
    % Set appropriate time Range
    numSamples = deltaT*Fs;
    timeVals = timeStartFromBaseLine+ (1/Fs:1/Fs:deltaT); %#ok<NASGU>

    % The optional input analogChannelToStore contains the list of channels that
    % should be recorded
    
    if electrodeCount ~= 0
        if isempty(analogElectrodesToStore)
            disp('Analog electrode list not given, taking all available electrodes.');
            electrodeListIDsStored = electrodeListIDs;
            electrodesStored = electrodeNums;
        else
            disp(['LFP from the specified ' num2str(length(analogElectrodesToStore)) ' electrodes and ' num2str(ainpCount) ' Ainp channels will be stored.']);

            count=1;
            for i=1:length(analogElectrodesToStore)
                % Check that the selected channel actually exist
                findChannelLoc = find(electrodeNums == analogElectrodesToStore(i));
                if isempty(findChannelLoc)
                    disp(['Analog Channel ' num2str(analogElectrodesToStore(i)) ' does not exist.']);
                else
                    electrodeListIDsStored(count) = electrodeListIDs(findChannelLoc); %#ok<AGROW>
                    electrodesStored(count) = electrodeNums(findChannelLoc); %#ok<AGROW>
                    count=count+1;
                end
            end
        end
        cElectrodeListIDsStored = length(electrodeListIDsStored);
        disp([num2str(cElectrodeListIDsStored) ' electrodes obtained']);
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare folders
    isWin=1; % variable for makeDirectory
    makeDirectory(folderOut,isWin); % main directory to store both LFP and spikes

    % Make Diectory for storing LFP data
    outputFolder = [folderOut 'LFP\'];
    makeDirectory(outputFolder,isWin);

    % Now segment and store data in the outputFolder directory
    totalStim = length(analysisOnsetTimes);
    goodStimPos = zeros(1,totalStim);
    for i=1:totalStim
        [~,goodStimPos(i)] = ns_GetIndexByTime(hFile,analogList(1),analysisOnsetTimes(i),-1);
    end

    %%%%%%%%%%%%%%%%%%%%%%% Get data from electrodes %%%%%%%%%%%%%%%%%%%%%%
    if electrodeCount ~= 0
        if cElectrodeListIDsStored > 0
            for i=1:cElectrodeListIDsStored
                disp(['elec' num2str(electrodesStored(i))]);

                clear analogInfo
                [~, analogInfo] = ns_GetAnalogInfo(hFile, electrodeListIDsStored(i)); %#ok<NASGU>

                clear analogData
                analogData = zeros(totalStim,numSamples);
                for j=1:totalStim
                    [~, ~, analogData(j,:)] = ns_GetAnalogData(hFile, ...
                        electrodeListIDsStored(i), goodStimPos(j)+1, numSamples);
                end
                save([outputFolder 'elec' num2str(electrodesStored(i))],'analogData','analogInfo');
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%% Get analog input data %%%%%%%%%%%%%%%%%%%%%%%%%
    if ainpCount>0
        for i=1:ainpCount
            disp(['ainp' num2str(analogInputNums(i))]);
            
            clear analogInfo
            [~, analogInfo] = ns_GetAnalogInfo(hFile, analogInputListIDs(i)); %#ok<NASGU>
            
            clear analogData
            analogData = zeros(totalStim,numSamples);
            for j=1:totalStim
                [~, ~, analogData(j,:)] = ns_GetAnalogData(hFile, ...
                    analogInputListIDs(i), goodStimPos(j)+1, numSamples);
            end
            save([outputFolder 'ainp' num2str(analogInputNums(i))],'analogData','analogInfo');
        end
    end

    % Write LFP information. For backward compatibility, we also save
    % analogChannelsStored which is the list of electrode data
    if electrodeCount == 0
        electrodesStored = 0;
    end
    analogChannelsStored = electrodesStored; %#ok<NASGU>
    save([outputFolder 'lfpInfo.mat'],'analogChannelsStored','electrodesStored','analogInputNums','goodStimPos','timeVals');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spike Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if electrodeCount ~= 0
    if getSpikes

        % The optional input digitalChannelsToStore contains the list of channels that
        % should be recorded
        clear neuralListIDs neuralChannelsStored
        if ~exist('neuralChannelsToStore','var')
            disp('Neural electrode list not given, taking all available electrodes.');
            neuralChannelsToStore = unique(neuralChannelNums);
        end

        disp(['Spikes from the specified ' num2str(length(neuralChannelsToStore)) ' electrodes will be stored.']);
        count=1;
        for i=1:length(neuralChannelsToStore)
            % Check that the selected channel actually exist
            findChannelLocs = find(neuralChannelNums == neuralChannelsToStore(i));
            if isempty(findChannelLocs)
                disp(['Neural Channel ' num2str(neuralChannelsToStore(i)) ' does not exist.']);
            else
                for j=1:length(findChannelLocs)
                    neuralListIDs(count) = neuralList(findChannelLocs(j)); %#ok<AGROW>
                    neuralChannelsStored(count) = neuralChannelNums(findChannelLocs(j)); %#ok<AGROW>
                    count=count+1;
                end
            end
        end
        
        if count == 1 % count is 1 if no neural channel is stored in the loop above
            cNeuralListIDs = 0;
            neuralChannelsStored = [];
        else
            cNeuralListIDs = length(neuralListIDs);
        end
        disp([num2str(cNeuralListIDs) ' neurons obtained']);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare folders
        isWin=1; % variable for makeDirectory
        makeDirectory(folderOut,isWin); % main directory to store both LFP and spikes

        % Make Diectory for storing Spike data
        outputFolder = [folderOut 'Spikes\'];
        makeDirectory(outputFolder,isWin);

        % Now segment and store data in the outputFolder directory
        totalStim = length(analysisOnsetTimes);

        SourceUnitID = zeros(1,cNeuralListIDs);
        for i=1:cNeuralListIDs
            clear neuralInfo
            [~, neuralInfo] = ns_GetNeuralInfo(hFile, neuralListIDs(i));
            SourceUnitID(i) = neuralInfo.SourceUnitID;
            disp([neuralChannelsStored(i) SourceUnitID(i)]);

            clear neuralData
            [~,neuralData]  = ns_GetNeuralData(hFile, neuralListIDs(i),1,numberOfItems(neuralListIDs(i)));

            clear spikeData
            spikeData = cell(1,totalStim);
            for j=1:totalStim
                spikeData{j} = neuralData(intersect(find(neuralData>=analysisOnsetTimes(j))...
                    ,find(neuralData<analysisOnsetTimes(j)+deltaT)));
                if ~isempty(spikeData{j})
                    spikeData{j} = spikeData{j} - goodStimTimes(j);
                end
            end
            save([outputFolder 'elec' num2str(neuralChannelsStored(i)) ... 
                '_SID' num2str(SourceUnitID(i))],'spikeData','neuralInfo');
        end

        % Write Spike information
        save([outputFolder 'spikeInfo.mat'],'neuralChannelsStored','SourceUnitID');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Segment Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if getSpikes % Get segments if spikes are being recorded
        segmentChannelsToStore = neuralChannelsToStore;

        % The optional input segmentChannelsToStore contains the list of channels that
        % should be recorded
        clear segmentListIDs segmentChannelsStored
        if ~exist('segmentChannelsToStore','var')
            disp('Segment list not given, taking all available electrodes.');
            segmentChannelsToStore = unique(segmentChannelNums);
        end

        disp(['Segments from the specified ' num2str(length(segmentChannelsToStore)) ' electrodes will be stored.']);
        count=1;
        for i=1:length(segmentChannelsToStore)
            % Check that the selected channel actually exist
            findChannelLocs = find(segmentChannelNums == segmentChannelsToStore(i));
            if isempty(findChannelLocs)
                disp(['Segements from Channel ' num2str(segmentChannelsToStore(i)) ' does not exist.']);
            else
                for j=1:length(findChannelLocs)
                    segmentListIDs(count) = segmentList(findChannelLocs(j)); %#ok<AGROW>
                    segmentChannelsStored(count) = segmentChannelNums(findChannelLocs(j)); %#ok<AGROW>
                    count=count+1;
                end
            end
        end

        cSegmentListIDs = length(segmentListIDs);
        disp([num2str(cSegmentListIDs) ' segments obtained']);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare folders
        isWin=1; % variable for makeDirectory
        makeDirectory(folderOut,isWin); % main directory to store both LFP and spikes

        % Make Diectory for storing Spike data
        outputFolder = [folderOut 'Segments\'];
        makeDirectory(outputFolder,isWin);

        % Now segment and store data in the outputFolder directory

        for i=1:cSegmentListIDs
            clear segmentInfo
            [~, segmentInfo] = ns_GetSegmentInfo(hFile, segmentListIDs(i)); %#ok<NASGU>
            %disp([segmentChannelsStored(i)]);

            clear timeStamp segmentData sampleCount unitID
            [~,timeStamp,segmentData, sampleCount, unitID]  = ns_GetSegmentData(hFile, segmentListIDs(i),1:min(1000000,numberOfItems(segmentListIDs(i)))); %#ok<NASGU,ASGLU>

            save([outputFolder 'elec' num2str(segmentChannelsStored(i))],'segmentData','segmentInfo','timeStamp','unitID','sampleCount');
        end

        % Write Spike information
        numItems = numberOfItems(segmentListIDs); %#ok<NASGU>
        save([outputFolder 'segmentInfo.mat'],'segmentChannelsStored','numItems');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close file
ns_CloseFile(hFile);
clear mexprog;
end