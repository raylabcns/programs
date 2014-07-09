% This program contains the following functions from the old format
% 1. saveLLDataSRC and saveLLDataGRF
% 2. getStimResultsLLSRC and getStimResultsLLGRF
% 3. saveEyeAndBehaviorDataSRC and saveEyeAndBehaviorDataGRF
% 4. getEyePositionAndBehavioralDataSRC and getEyePositionAndBehavioralDataGRF

% 5. saveEyeDataInDegSRC, saveEyeDataInDeg (GRF suffix added),
% saveEyeDataStimPosGRF and getEyeDataStimPosGRF

function saveLLData(monkeyName,expDate,protocolName,folderSourceString,gridType,type)

FsEye=200;
eyeRangeMS{1} = [-320 320]; timePeriodMS{1} = [-500   500]; maxStimPos{1} = 20;
eyeRangeMS{2} = [-480 800]; timePeriodMS{2} = [-1000 1000]; maxStimPos{2} = 10;

folderName    = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];
folderExtract = [folderName 'extractedData'];
makeDirectory(folderExtract);

if strncmpi(protocolName,'SRC',3) % SRC
    [LL,targetInfo,psyInfo,reactInfo] = getStimResultsLLSRC(monkeyName,expDate,protocolName,folderSourceString); %#ok<*ASGLU,*NASGU>
    save([folderExtract '\LL.mat'],'LL','targetInfo','psyInfo','reactInfo');
    
    [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataSRC(monkeyName,expDate,protocolName,folderSourceString,eyeRangeMS{type},FsEye);
    save([folderExtract '\BehaviorData.mat'],'allTrials','goodTrials','stimData');
    save([folderExtract '\EyeData.mat'],'eyeData','eyeRangeMS');
    
    saveEyeDataInDegSRC(monkeyName,expDate,protocolName,folderSourceString,gridType);
else
    LL = getStimResultsLLGRF(monkeyName,expDate,protocolName,folderSourceString);
    save([folderExtract '\LL.mat'],'LL');

    [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataGRF(monkeyName,expDate,protocolName,folderSourceString,eyeRangeMS{type},FsEye);
    save([folderExtract '\BehaviorData.mat'],'allTrials','goodTrials','stimData');
    save([folderExtract '\EyeData.mat'],'eyeData','eyeRangeMS');
    
    saveEyeDataInDegGRF(monkeyName,expDate,protocolName,folderSourceString,gridType,timePeriodMS{type},FsEye,maxStimPos{type});
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LL,targetInfo,psyInfo,reactInfo] = getStimResultsLLSRC(monkeyName,expDate,protocolName,folderSourceString)

frameISIinMS=10;

monkeyName = removeIfPresent(monkeyName,'\');
expDate    = removeIfPresent(expDate,'\');
protocolName = removeIfPresent(protocolName,'\');
folderSourceString = appendIfNotPresent(folderSourceString,'\');

datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];

% Get Lablib data
header = readLLFile('i',datFileName);

if isfield(header,'eccentricityDeg')
    LL.azimuthDeg   = round(100*header.eccentricityDeg.data*cos(header.polarAngleDeg.data/180*pi))/100;
    LL.elevationDeg = round(100*header.eccentricityDeg.data*sin(header.polarAngleDeg.data/180*pi))/100;
else
    LL.azimuthDeg   = header.azimuthDeg.data;
    LL.elevationDeg = header.elevationDeg.data;
end

LL.sigmaDeg = header.sigmaDeg.data;
LL.spatialFreqCPD = header.spatialFreqCPD.data;
LL.orientationDeg = header.stimOrientationDeg.data;

% Stimulus properties
numTrials = header.numberOfTrials;
targetContrastIndex = [];
targetTemporalFreqIndex = [];

catchTrial = [];
instructTrial = [];
attendLoc = [];

eotCode=[];
startTime=[];

countTI=1;
countPI=1;
countRI=1;

for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    thisEOT = [trials.trialEnd.data];
    eotCode = cat(2,eotCode,thisEOT);
    
    if isfield(trials,'trialStart')
        startTime = cat(2,startTime,[trials.trialStart.timeMS]);
    end
    
    %     if isfield(trials,'stimulusOn')
    %         allStimulusIndex = [allStimulusIndex [trials.stimulusOn.data]'];
    %     end
    %
    %     if isfield(trials,'stimulusOnTime')
    %         allStimulusOnTimes = [allStimulusOnTimes [trials.stimulusOnTime.timeMS]'];
    %     end
    
    if isfield(trials,'trial')
        thisCatchTrial    = [trials.trial.data.catchTrial];
        thisInstructTrial = [trials.trial.data.instructTrial];
        thisAttendLoc     = [trials.trial.data.attendLoc];
        thisTargetContrastIndex     = [trials.trial.data.targetContrastIndex];
        thisTargetTemporalFreqIndex = [trials.trial.data.targetTemporalFreqIndex];
        
        catchTrial    = cat(2,catchTrial,thisCatchTrial);
        instructTrial = cat(2,instructTrial,thisInstructTrial);
        attendLoc     = cat(2,attendLoc,thisAttendLoc);
        
        targetContrastIndex    = cat(2,targetContrastIndex,thisTargetContrastIndex);
        targetTemporalFreqIndex = cat(2,targetTemporalFreqIndex,thisTargetTemporalFreqIndex);
        
        
        % Adding the getTrialInfoLLFile (from the MAC) details here
        
        
        % Nothing to update on catch trials or instruction trials
        if (thisCatchTrial) || (thisInstructTrial) %#ok<*BDSCI,BDLGI>
            % Do nothing
        else
            
            % All the information needed to plot the targetIndex versus
            % performace plot are stored in the variable targetInfo
            
            targetInfo(countTI).trialNum          = i; %#ok<*AGROW>
            targetInfo(countTI).contrastIndex     = thisTargetContrastIndex;
            targetInfo(countTI).temporalFreqIndex = thisTargetTemporalFreqIndex;
            targetInfo(countTI).targetIndex       = trials.trial.data.targetIndex;
            targetInfo(countTI).attendLoc         = thisAttendLoc;
            targetInfo(countTI).eotCode           = thisEOT;
            countTI=countTI+1;
            
            % In addition, if the EOTCODE is correct, wrong or failed, make
            % another variable 'psyInfo' that stores information to plot the
            % psychometric functions
            
            if thisEOT <=2 % Correct, wrong or failed
                psyInfo(countPI).trialNum          = i;
                psyInfo(countPI).contrastIndex     = thisTargetContrastIndex;
                psyInfo(countPI).temporalFreqIndex = thisTargetTemporalFreqIndex;
                psyInfo(countPI).attendLoc         = thisAttendLoc;
                psyInfo(countPI).eotCode           = thisEOT;
                countPI=countPI+1;
            end
            
            % Finally, get the reaction times for correct trials and store them
            % in 'reactInfo'
            
            if thisEOT==0
                reactInfo(countRI).trialNum          = i;
                reactInfo(countRI).contrastIndex     = thisTargetContrastIndex;
                reactInfo(countRI).temporalFreqIndex = thisTargetTemporalFreqIndex;
                reactInfo(countRI).attendLoc         = thisAttendLoc;
                reactInfo(countRI).reactTime         = trials.saccade.timeMS-trials.visualStimsOn.timeMS ...
                                                       -trials.stimulusOn.data(trials.trial.data.targetIndex+1)*frameISIinMS;
                countRI=countRI+1;
            end
        end
    end
end

LL.eotCode = eotCode;
LL.instructTrial = instructTrial;
LL.catchTrial = catchTrial;
LL.attendLoc = attendLoc;
LL.targetContrastIndex = targetContrastIndex;
LL.targetTemporalFreqIndex = targetTemporalFreqIndex;

LL.startTime = startTime/1000; % in seconds
end
function LL = getStimResultsLLGRF(monkeyName,expDate,protocolName,folderSourceString)

monkeyName = removeIfPresent(monkeyName,'\');
expDate    = removeIfPresent(expDate,'\');
protocolName = removeIfPresent(protocolName,'\');
folderSourceString = appendIfNotPresent(folderSourceString,'\');

datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];

% Get Lablib data
header = readLLFile('i',datFileName);

% Stimulus properties
numTrials = header.numberOfTrials;

allStimulusIndex = [];
allStimulusOnTimes = [];
gaborIndex = [];
stimType = [];
azimuthDeg = [];
elevationDeg = [];
sigmaDeg = [];
radiusDeg = [];
spatialFreqCPD =[];
orientationDeg = [];
contrastPC = [];
temporalFreqHz = [];

eotCode=[];
myEotCode=[];
startTime=[];

for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        eotCode = [eotCode [trials.trialEnd.data]];
    end
    
    if isfield(trials,'myTrialEnd')
        myEotCode = [myEotCode [trials.myTrialEnd.data]];
    end
    
    if isfield(trials,'trialStart')
        startTime = [startTime [trials.trialStart.timeMS]];
    end
    
    if isfield(trials,'stimulusOn')
        allStimulusIndex = [allStimulusIndex [trials.stimulusOn.data]'];
    end
    
    if isfield(trials,'stimulusOnTime')
        allStimulusOnTimes = [allStimulusOnTimes [trials.stimulusOnTime.timeMS]'];
    end
    
    if isfield(trials,'stimDesc')
        gaborIndex = [gaborIndex [trials.stimDesc.data.gaborIndex]];
        stimType = [stimType [trials.stimDesc.data.stimType]];
        azimuthDeg = [azimuthDeg [trials.stimDesc.data.azimuthDeg]];
        elevationDeg = [elevationDeg [trials.stimDesc.data.elevationDeg]];
        sigmaDeg = [sigmaDeg [trials.stimDesc.data.sigmaDeg]];
        if isfield(trials.stimDesc.data,'radiusDeg')
            radiusExists=1;
            radiusDeg = [radiusDeg [trials.stimDesc.data.radiusDeg]];
        else
            radiusExists=0;
        end
        spatialFreqCPD = [spatialFreqCPD [trials.stimDesc.data.spatialFreqCPD]];
        orientationDeg = [orientationDeg [trials.stimDesc.data.directionDeg]];
        contrastPC = [contrastPC [trials.stimDesc.data.contrastPC]];
        temporalFreqHz = [temporalFreqHz [trials.stimDesc.data.temporalFreqHz]];
    end
end

% Sort stim properties by stimType
for i=1:3
    gaborIndexFromStimulusOn{i} = find(allStimulusIndex==i-1);
    gaborIndexFromStimDesc{i} = find(gaborIndex==i-1);
end

if isequal(gaborIndexFromStimDesc,gaborIndexFromStimulusOn)
    for i=1:3
        LL.(['time' num2str(i-1)]) = allStimulusOnTimes(gaborIndexFromStimulusOn{i});
        LL.(['stimType' num2str(i-1)]) = stimType(gaborIndexFromStimulusOn{i});
        LL.(['azimuthDeg' num2str(i-1)]) = azimuthDeg(gaborIndexFromStimulusOn{i});
        LL.(['elevationDeg' num2str(i-1)]) = elevationDeg(gaborIndexFromStimulusOn{i});
        LL.(['sigmaDeg' num2str(i-1)]) = sigmaDeg(gaborIndexFromStimulusOn{i});
        if radiusExists
            LL.(['radiusDeg' num2str(i-1)]) = radiusDeg(gaborIndexFromStimulusOn{i});
        end
        LL.(['spatialFreqCPD' num2str(i-1)]) = spatialFreqCPD(gaborIndexFromStimulusOn{i});
        LL.(['orientationDeg' num2str(i-1)]) = orientationDeg(gaborIndexFromStimulusOn{i});
        LL.(['contrastPC' num2str(i-1)]) = contrastPC(gaborIndexFromStimulusOn{i});
        LL.(['temporalFreqHz' num2str(i-1)]) = temporalFreqHz(gaborIndexFromStimulusOn{i});
    end
else
    error('Gabor indices from stimuluOn and stimDesc do not match!!');
end

LL.eotCode = eotCode;
LL.myEotCode = myEotCode;
LL.startTime = startTime/1000; % in seconds
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataSRC(monkeyName,expDate,protocolName,folderSourceString,eyeRangeMS,Fs)

if ~exist('eyeRangeMS','var')           eyeRangeMS = [-480 800];        end    % ms
if ~exist('Fs','var')                   Fs = 200;                       end % Eye position sampled at 200 Hz.

eyeRangePos = eyeRangeMS*Fs/1000;

monkeyName = removeIfPresent(monkeyName,'\');
expDate    = removeIfPresent(expDate,'\');
protocolName = removeIfPresent(protocolName,'\');
folderSourceString = appendIfNotPresent(folderSourceString,'\');

datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];

% Get Lablib data
header = readLLFile('i',datFileName);

% Stimulus properties
numTrials = header.numberOfTrials;
stimNumber=1;
correctIndex=1;
trialEndIndex=1;

for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        allTrials.trialEnded(i) = 1;
        allTrials.catchTrials(trialEndIndex) = trials.trial.data.catchTrial;
        allTrials.instructTrials(trialEndIndex) = trials.trial.data.instructTrial;
        allTrials.trialCertify(trialEndIndex) = trials.trialCertify.data;
        allTrials.targetPosAllTrials(trialEndIndex) = trials.trial.data.targetIndex+1;
        allTrials.eotCodes(trialEndIndex) = trials.trialEnd.data;
        
        allTrials.fixWindowSize(trialEndIndex) = trials.fixWindowData.data.windowDeg.size.width;
        %allTrials.respWindowSize(trialEndIndex) = trials.respWindowData.data.windowDeg.size.width; % instead of responseWindowData
        allTrials.certifiedNonInstruction(trialEndIndex) = (allTrials.instructTrials(trialEndIndex)==0)*(allTrials.trialCertify(trialEndIndex)==0);

        if (allTrials.eotCodes(trialEndIndex)==0) &&  (allTrials.certifiedNonInstruction(trialEndIndex)==1) ...
                && (allTrials.catchTrials(trialEndIndex)==0) % Work on only Correct Trials, which are not instruction, catch or uncertified trials
            
            % Get Eye Data
            eyeX = trials.eyeXData.data;
            eyeY = trials.eyeYData.data;
            % eyeStartTime = trials.eyeXData.timeMS(1);  % This is wrong.
            % The eye data is synchronized with trialStartTime.
            eyeStartTime = trials.trialStart.timeMS;
            eyeAllTimes = eyeStartTime + (0:(length(eyeX)-1))*(1000/Fs);
            
            stimOnTimes  = [trials.stimulusOn.timeMS];
            numStimuli = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
            
            goodTrials.targetPos(correctIndex) = numStimuli;
            goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS;
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            clear stimType
            stimType = ([trials.stimDesc.data.type0]) .* ([trials.stimDesc.data.type1]);
            goodTrials.stimType{correctIndex} = stimType;
            
            for j=1:numStimuli

                if (stimType(j)==9) || (stimType(j)==1)  % Frontpad or Valid
                    
                    stimTime = stimOnTimes(j);
                    stp=find(eyeAllTimes>=stimTime, 1 );
                    
                    stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    stimData.stimPos(stimNumber) = j;
                    stimData.stimType(stimNumber) = stimType(j);
                    
                    if (stimType(j)==9) % First stimulus may not have sufficient baseline
                        eyeData(stimNumber).eyePosDataX = eyeX(stp:stp+eyeRangePos(2)-1);
                        eyeData(stimNumber).eyePosDataY = eyeY(stp:stp+eyeRangePos(2)-1);
                    else
                        eyeData(stimNumber).eyePosDataX = eyeX(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);
                        eyeData(stimNumber).eyePosDataY = eyeY(stp+eyeRangePos(1):stp+eyeRangePos(2)-1);  
                    end
                    
                    eyeData(stimNumber).eyeCal = trials.eyeCalibrationData.data.cal;
                    stimNumber=stimNumber+1;
                end
            end
            
            correctIndex=correctIndex+1;
        end
        trialEndIndex=trialEndIndex+1;
    end
end
end
function [allTrials,goodTrials,stimData,eyeData,eyeRangeMS] = getEyePositionAndBehavioralDataGRF(monkeyName,expDate,protocolName,folderSourceString,eyeRangeMS,Fs)

if ~exist('eyeRangeMS','var')           eyeRangeMS = [-480 800];        end    % ms
if ~exist('Fs','var')                   Fs = 200;                       end % Eye position sampled at 200 Hz.

eyeRangePos = eyeRangeMS*Fs/1000;

if ~exist('folderSourceString','var')   folderSourceString ='F:\';       end

monkeyName = removeIfPresent(monkeyName,'\');
expDate    = removeIfPresent(expDate,'\');
protocolName = removeIfPresent(protocolName,'\');
folderSourceString = appendIfNotPresent(folderSourceString,'\');

datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];

% Get Lablib data
header = readLLFile('i',datFileName);

% Stimulus properties
numTrials = header.numberOfTrials;
stimNumber=1;
correctIndex=1;
trialEndIndex=1;

for i=1:numTrials
    disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trialEnd')
        allTrials.trialEnded(i) = 1;
        allTrials.catchTrials(trialEndIndex) = trials.trial.data.catchTrial;
        allTrials.instructTrials(trialEndIndex) = trials.trial.data.instructTrial;
        allTrials.trialCertify(trialEndIndex) = trials.trialCertify.data;
        allTrials.targetPosAllTrials(trialEndIndex) = trials.trial.data.targetIndex+1;
        allTrials.eotCodes(trialEndIndex) = trials.trialEnd.data;
        
        allTrials.fixWindowSize(trialEndIndex) = trials.fixWindowData.data.windowDeg.size.width;
        allTrials.respWindowSize(trialEndIndex) = trials.responseWindowData.data.windowDeg.size.width;
        allTrials.certifiedNonInstruction(trialEndIndex) = (allTrials.instructTrials(trialEndIndex)==0)*(allTrials.trialCertify(trialEndIndex)==0);

        if (allTrials.eotCodes(trialEndIndex)==0) &&  (allTrials.certifiedNonInstruction(trialEndIndex)==1)
                %&& (allTrials.catchTrials(trialEndIndex)==0) % Work on only Correct Trials, which are not instruction or uncertified trials. Include catch trials
            
            isCatchTrial = (allTrials.catchTrials(trialEndIndex)==1);
            
            % Get Eye Data
            eyeX = trials.eyeXData.data;
            eyeY = trials.eyeYData.data;
            % eyeStartTime = trials.eyeXData.timeMS(1);  % This is wrong.
            % The eye data is synchronized with trialStartTime.
            eyeStartTime = trials.trialStart.timeMS;
            eyeAllTimes = eyeStartTime + (0:(length(eyeX)-1))*(1000/Fs);
            
            stimOnTimes  = [trials.stimulusOnTime.timeMS];
            numStimuli = allTrials.targetPosAllTrials(trialEndIndex); %=length(stimOnTimes)/3;
            
            goodTrials.targetPos(correctIndex) = numStimuli;
            goodTrials.targetTime(correctIndex) = stimOnTimes(end);
            goodTrials.fixateMS(correctIndex) = trials.fixate.timeMS;
            goodTrials.fixonMS(correctIndex) = trials.fixOn.timeMS;
            goodTrials.stimOnTimes{correctIndex} = stimOnTimes;
            
            % Find position of Gabor1
            gaborPos = find([trials.stimDesc.data.gaborIndex]==1); % could be 4 gabors for GRF protocol
            
            if isCatchTrial
                stimEndIndex = numStimuli;    % Take the last stimulus because it is still a valid stimulus
            else
                stimEndIndex = numStimuli-1;  % Don't take the last one because it is the target
            end
                
            if stimEndIndex>0  % At least one stimulus
                for j=1:stimEndIndex
                    stimTime = stimOnTimes(gaborPos(j));
                    stp=find(eyeAllTimes>=stimTime,1);
                    
                    stimData.stimOnsetTimeFromFixate(stimNumber) = stimTime-trials.fixate.timeMS;
                    stimData.stimPos(stimNumber) = j;
                    
                    startingPos = max(1,stp+eyeRangePos(1)); % First stimulus may not have sufficient baseline           
                    endingPos   = min(stp+eyeRangePos(2)-1,length(eyeX)); % Last stimulus may not have suffiecient length, e.g. for a catch trial
     
                    eyeData(stimNumber).eyePosDataX = eyeX(startingPos:endingPos);
                    eyeData(stimNumber).eyePosDataY = eyeY(startingPos:endingPos);
                    
                    eyeData(stimNumber).eyeCal = trials.eyeCalibrationData.data.cal;
                    stimNumber=stimNumber+1;
                end
            end
            correctIndex=correctIndex+1;
        end
        trialEndIndex=trialEndIndex+1;
    end
end
end
% Additional analysis
function saveEyeDataInDegSRC(monkeyName,expDate,protocolName,folderSourceString,gridType)
% The difference between saveEyeDataInDeg and saveEyeDataInDegSRC is that
% the frontPad stimuli are also saved in SRC.

folderName    = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];
folderExtract = [folderName 'extractedData\'];

clear eyeData 
load([folderExtract 'EyeData.mat']);

[eyeDataDegX,eyeDataDegY] = convertEyeDataToDeg(eyeData,1);

for i=1:length(eyeDataDegX)
    lengthEyeSignal = size(eyeDataDegX{i},1);
    eyeSpeedX{i} = [eyeDataDegX{i}(2:lengthEyeSignal)-eyeDataDegX{i}(1:lengthEyeSignal-1);0];
    eyeSpeedY{i} = [eyeDataDegY{i}(2:lengthEyeSignal)-eyeDataDegY{i}(1:lengthEyeSignal-1);0];
end

folderSave = [folderName 'segmentedData\eyeData\'];
makeDirectory(folderSave);
save([folderSave 'eyeDataDeg.mat'],'eyeDataDegX','eyeDataDegY');
save([folderSave 'eyeSpeed.mat'],'eyeSpeedX','eyeSpeedY');

end
function saveEyeDataInDegGRF(monkeyName,expDate,protocolName,folderSourceString,gridType,timePeriodMS,FsEye,maxStimPos)

folderName    = [folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\'];
folderExtract = [folderName 'extractedData\'];

clear eyeData 
load([folderExtract 'EyeData.mat']);

clear goodStimNums
load([folderExtract 'goodStimNums.mat']);

if exist([folderExtract 'validStimAfterTarget.mat'],'file')
    load([folderExtract 'validStimAfterTarget.mat']);
    disp(['Removing ' num2str(length(validStimuliAfterTarget)) ' stimuli from goodStimNums']);
    goodStimNums(validStimuliAfterTarget)=[];
end

clear stimResults
load([folderExtract 'stimResults.mat']);

goodStimPos = stimResults.stimPosition(goodStimNums);

% all stimPostions greater than 1
useTheseStims = find(goodStimPos>1);

[eyeDataDegX,eyeDataDegY] = convertEyeDataToDeg(eyeData(useTheseStims),1);
folderSave = [folderName 'segmentedData\eyeData\'];
makeDirectory(folderSave);
save([folderSave 'eyeDataDeg.mat'],'eyeDataDegX','eyeDataDegY');

% lengthEyeSignal = size(eyeDataDegX,2);
% eyeSpeedX = [eyeDataDegX(:,2:lengthEyeSignal)-eyeDataDegX(:,1:lengthEyeSignal-1) zeros(size(eyeDataDegX,1),1)];
% eyeSpeedY = [eyeDataDegY(:,2:lengthEyeSignal)-eyeDataDegY(:,1:lengthEyeSignal-1) zeros(size(eyeDataDegY,1),1)];
% save([folderSave 'eyeSpeed.mat'],'eyeSpeedX','eyeSpeedY');

% More data saved for GRF protocol
[eyeXAllPos,eyeYAllPos,xs] = getEyeDataStimPosGRF(monkeyName,expDate,protocolName,folderSourceString,timePeriodMS,FsEye,maxStimPos);
save([folderSave 'EyeDataStimPos.mat'],'eyeXAllPos','eyeYAllPos','xs','timePeriodMS');
end
function [eyeXAllPos,eyeYAllPos,xs,timePeriodMS] = getEyeDataStimPosGRF(monkeyName,expDate,protocolName,folderSourceString,timePeriodMS,Fs,maxStimPos)

intervalTimeMS=1000/Fs;

datFileName = [folderSourceString 'data\rawData\' monkeyName expDate '\' monkeyName expDate protocolName '.dat'];

% Get Lablib data
header = readLLFile('i',datFileName);

for j=1:maxStimPos                                  
    eyeXAllPos{j}=[];
    eyeYAllPos{j}=[];
    xs{j}=[];
end
    
for i=1:header.numberOfTrials
    trial = readLLFile('t',i);
    
    % Work on only Correct Trials, which are not instruction, catch or
    % uncertified trials
    if (trial.trialEnd.data == 0) && (trial.trial.data.instructTrial==0) && ...
            (trial.trialCertify.data==0) && (trial.trial.data.catchTrial==0)
        
        % get eye data
        eX=trial.eyeXData.data';
        eY=trial.eyeYData.data';
        cal=trial.eyeCalibrationData.data.cal;
        numUsefulStim = trial.trial.data.targetIndex; % these are the useful stimuli, excluding target
        
        stimOnTimes = trial.stimulusOnTime.timeMS;
        gaborPos = find([trial.stimDesc.data.gaborIndex]==1);
        if numUsefulStim>1
            for j=2:numUsefulStim
                stp = ceil((stimOnTimes(gaborPos(j)) - trial.trialStart.timeMS)/intervalTimeMS);
                
                eXshort = eX( stp + ((j-1)*timePeriodMS(1))/intervalTimeMS : stp + timePeriodMS(2)/intervalTimeMS);
                eYshort = eY( stp + ((j-1)*timePeriodMS(1))/intervalTimeMS : stp + timePeriodMS(2)/intervalTimeMS);
                
                eyeX = cal.m11*eXshort + cal.m21 * eYshort + cal.tX;
                eyeY = cal.m12*eXshort + cal.m22 * eYshort + cal.tY;
                
                eyeXAllPos{j} = cat(1,eyeXAllPos{j},eyeX);
                eyeYAllPos{j} = cat(1,eyeYAllPos{j},eyeY);
                
                if isempty(xs{j})
                    xs{j} = 0:intervalTimeMS:timePeriodMS(2) - (j-1)*timePeriodMS(1);
                end
            end
        end
    end
end
end