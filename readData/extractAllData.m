% This function is used to extract the following from raw data files
% 1. Digital data
% 2. LFP data
% 3. Spike data
% 4. Segments

% Each data file is characterized by four parameters - monkeyName, expDate,
% protocolName and gridType.

% We assume that the raw data is initially stored in
% folderSourceString\data\rawData\{monkeyName}{expDate}\

function extractAllData(monkeyName,expDate,protocolName,folderSourceString,gridType,timeStartFromBaseLine,deltaT,electrodesToStore)

if ~exist('folderSourceString','var')   folderSourceString ='F:\';      end
if ~exist('timeStartFromBaseLine','var') timeStartFromBaseLine= -0.55;  end
if ~exist('deltaT','var')                deltaT = 1.024;                end

folderSourceString = appendIfNotPresent(folderSourceString,'\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileName = [monkeyName expDate protocolName '.nev'];
folderName0 = [folderSourceString 'data\' monkeyName '\'];
makeDirectory(folderName0);
folderName0 = [folderName0 gridType '\'];
makeDirectory(folderName0);
folderName1 = [folderName0 expDate '\'];
makeDirectory(folderName1);
folderName = [folderName1 protocolName '\'];
makeDirectory(folderName);

folderIn = [folderSourceString 'data\rawData\' monkeyName expDate '\'];
folderExtract = [folderName 'extractedData\'];
folderSegment = [folderName 'segmentedData'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the NEV file

% Load the appropriate DLL
dllName = 'N:\programs\requiredResources\nsNEVLibrary64.dll';
[nsresult] = ns_SetLibrary(dllName);
if (nsresult ~= 0)      error('DLL was not found!');                    end

% Load data file and display some info about the file open data file
[nsresult, hFile] = ns_OpenFile([folderIn fileName]);
if (nsresult ~= 0)      error('Data file did not open!');               end

% Get file information
[nsresult, fileInfo] = ns_GetFileInfo(hFile);
% Gives you entityCount, timeStampResolution and timeSpan
if (nsresult ~= 0)      error('Data file information did not load!');   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Digital Codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, entityInfo] = ns_GetEntityInfo(hFile,1:fileInfo.EntityCount);
eventList = find([entityInfo.EntityType] == 1);

digitalEventName = 'digin';
for i =1:length(eventList)
    if strcmp(digitalEventName,[entityInfo(eventList(i)).EntityLabel]);
        digID=i;
        break;
    end
end

if ~exist('digID','var')    
    error(['No event named ' digitalEventName]);
else
    numDigitalEvents=entityInfo(eventList(digID)).ItemCount;
end
disp(['Number of digital events: ' num2str(numDigitalEvents)]);

% Get the digital events
[~,digitalTimeStamps,digitalEvents] = ns_GetEventData(hFile,eventList(digID),1:numDigitalEvents);

% Blackrock collects data at 30,000 samples per second. Sometimes a digital
% code transition takes longer than this, and is counted twice. First we
% find out the double counts.

deltaLimit = 1.5/30000; 
dt = diff(digitalTimeStamps);
badDTPos = find(dt<=deltaLimit);

if ~isempty(badDTPos)
    disp([num2str(length(badDTPos)) ' of ' num2str(length(digitalTimeStamps)) ' (' num2str(100*length(badDTPos)/length(digitalTimeStamps),2) '%) are repeats and will be discarded']);
    digitalTimeStamps(badDTPos)=[];
    digitalEvents(badDTPos)=[];
end

goodStimTimes = extractDigitalDataGRF(digitalEvents,digitalTimeStamps,folderExtract,1);
save([folderExtract 'NEVFileInfo.mat'], 'fileInfo', 'entityInfo');

%%%%%%%%%%%%%%%%%%%%% Save Analog and Spike Data %%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=2000;
analogChannelsToStore = electrodesToStore;
neuralChannelsToStore = analogChannelsToStore;
getLFP=1;getSpikes=1;
getLFPandSpikes(fileName,analogChannelsToStore,folderIn,folderSegment, ...
    goodStimTimes,timeStartFromBaseLine,deltaT,Fs,hFile,neuralChannelsToStore,getLFP,getSpikes);
end