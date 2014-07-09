function electrodeNums = getElectrodesInRange(monkeyName,expDate,protocolName,folderSourceString,dRange,isWin,stimulusCenter,useTheseElectrodes,gridType)

if ~exist('isWin','var')                isWin=1;                        end
if ~exist('stimulusCenter','var')       stimulusCenter=[];              end
if ~exist('useTheseElectrodes','var')   useTheseElectrodes=[];          end
if ~exist('gridType','var')             gridType = 'Microelectrode';   end

if isWin==1
    folderseparator='\';
    folderCluster = 'N:\Programs\programsMatlabAddPath\generalPrograms\receptiveFieldEstimation\';
else
    folderseparator='/';
    folderCluster = '/Volumes/maunsell/';
end

% get Stimulus position information
if isempty(stimulusCenter)
    load([folderSourceString monkeyName folderseparator expDate folderseparator protocolName folderseparator 'extractedData' folderseparator 'parameterCombinations.mat']);
else
    aValsUnique=stimulusCenter(1);
    eValsUnique=stimulusCenter(2);
end

% get RF information
load([folderCluster monkeyName gridType 'RFData']);

count=1;
if isempty(useTheseElectrodes)
    useTheseElectrodes=highRMSElectrodes;
end

for i=1:length(useTheseElectrodes)
    azi = rfStats(useTheseElectrodes(i)).meanAzi;
    ele = rfStats(useTheseElectrodes(i)).meanEle;
    
    d = sqrt(sum((azi-aValsUnique)^2+(ele-eValsUnique)^2));
    
    if d>=dRange(1) && d<dRange(2)
        electrodeNums(count) = useTheseElectrodes(i); %#ok<AGROW>
        count=count+1;
    end
end

if count==1 % Nothing found
    electrodeNums=[];
end
end