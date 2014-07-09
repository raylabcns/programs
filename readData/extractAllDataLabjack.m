% Labjack Data is stored using the LJStreamUD application which comes with
% the installation package. It runs only on windows. This application
% stores analog data, but can also record a 16 digit digital code as an
% analog data. However, to read these codes properly, we need to sample the
% data very fast (otherwise multiple codes get sent between the sampling
% interval). This program assumes that the data is collected using
% LJStreamUD and the sampling rate is 20 KHz. Further, digital code is
% stored on channel 0 and the analog data (eg: output of the photometer) is
% stored in channel 1 onwards.

% We assume that the raw data is initially stored in
% folderSourceString\data\rawData\{monkeyName}{expDate}\

function extractAllDataLabjack(monkeyName,expDate,protocolName,folderSourceString,gridType,timeStartFromBaseLine,deltaT,electrodesToStore)

if ~exist('folderSourceString','var')   folderSourceString ='F:\';      end
if ~exist('timeStartFromBaseLine','var') timeStartFromBaseLine= -0.55;  end
if ~exist('deltaT','var')                deltaT = 1.024;                end

folderSourceString = appendIfNotPresent(folderSourceString,'\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderName0 = [folderSourceString 'data\' monkeyName '\'];
makeDirectory(folderName0);
folderName0 = [folderName0 gridType '\'];
makeDirectory(folderName0);
folderName1 = [folderName0 expDate '\'];
makeDirectory(folderName1);
folderName = [folderName1 protocolName '\'];
makeDirectory(folderName);

folderExtract = [folderName 'extractedData\'];
folderSegment = [folderName 'segmentedData\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Labjack Data
[data,t] = readLabjackData(monkeyName,expDate,protocolName,folderSourceString,length(electrodesToStore));

%%%%%%%%%%%%%%%%%%%%%%%%%% Digital Codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = data{1};
dd = abs(diff(d));
dPos = find(dd>0)+1;

digitalEventsAll = d(dPos);
digitalTimeStampsAll = t(dPos);

% Labjack collects data at 20,000 samples per second. Sometimes a digital
% code transition takes longer than this, and is counted twice. First we
% find out the double counts.

Fs = unique(round(1./diff(t)));

if Fs~= 20000
    error('Labjack must be sampled at 20 kHz');
end

deltaLimit = 1.5/Fs; 
dt = diff(digitalTimeStampsAll);
badDTPos = find(dt<=deltaLimit);

if ~isempty(badDTPos)
    disp([num2str(length(badDTPos)) ' of ' num2str(length(digitalTimeStampsAll)) ' (' num2str(100*length(badDTPos)/length(digitalTimeStampsAll),2) '%) are repeats and will be discarded']);
    digitalTimeStampsAll(badDTPos)=[];
    digitalEventsAll(badDTPos)=[];
end

goodStimTimes = extractDigitalDataGRF(digitalEventsAll,digitalTimeStampsAll,folderExtract,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Save LFP Data %%%%%%%%%%%%%%%%%%%%%%%%
saveAnalogData(data,t,electrodesToStore,folderSegment,goodStimTimes,timeStartFromBaseLine,deltaT,Fs);
end

% Read Data from Labjack
function [analogData,timeVals] = readLabjackData(monkeyName,expDate,protocolName,folderSourceString,numChannels)

folderIn = [folderSourceString 'data\rawData\' monkeyName expDate '\'];

timeVals = [];
for i=1:numChannels
    analogData{i} = [];
end

fileExistFlag = 1;
fileCounter = 0;

while(fileExistFlag)
    fileName = [folderIn monkeyName expDate protocolName '_' num2str(fileCounter) '.dat'];
    
    if exist(fileName,'file')
        disp(['Reading from ' fileName]);
        [v,t]=readLabjackDataSingleFile(fileName,numChannels);
    
        timeVals = cat(1,timeVals,t);
        for j=1:numChannels
            analogData{j} = cat(1,analogData{j},v{j});
        end
        fileCounter = fileCounter+1;
    else
        fileExistFlag = 0;
    end 
end
end
function [v,t,dataDetails]=readLabjackDataSingleFile(fileName,numChannels)

fid = fopen(fileName,'r');

% Get Parameters of interest

dataDetails.dateStr = fscanf(fid,'%s',1);
dataDetails.timeStr = fscanf(fid,'%s',1);

for i=1:numChannels
    dataDetails.channelDetails{i} = fscanf(fid,'%s',5);
end

dataDetails.paramsList = fscanf(fid,'%s',5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = fscanf(fid,'%f%f%f%f%f');

t = data(1:numChannels*2+1:end);
for i=1:numChannels
    v{i} = data(i+1:numChannels*2+1:end); %#ok<*AGROW>
end

fclose(fid);
end
function saveAnalogData(data,t,analogChannelsStored,folderSegment,goodStimTimes,timeStartFromBaseLine,deltaT,Fs)

% Make Diectory for storing LFP data
makeDirectory(folderSegment);
outputFolder = [folderSegment 'LFP\'];
makeDirectory(outputFolder);

analysisOnsetTimes = goodStimTimes + timeStartFromBaseLine;
numSamples = deltaT*Fs;
timeVals = timeStartFromBaseLine+ (1/Fs:1/Fs:deltaT); %#ok<NASGU>

% Now segment and store data in the outputFolder directory
totalStim = length(analysisOnsetTimes);

for i=1:length(analogChannelsStored)
    dataThisChannel = data{i};
    
    clear analogData
    analogData = zeros(totalStim,numSamples);
    
    for j=1:totalStim
        pos = find(t<=analysisOnsetTimes(j), 1, 'last' );
        analogData(j,:) = dataThisChannel(pos+1:pos+numSamples);
    end
    save([outputFolder 'elec' num2str(analogChannelsStored(i))],'analogData');
end

% Write LFP information
save([outputFolder 'lfpInfo.mat'],'analogChannelsStored','timeVals');

end