function params = getRFcenter(aVals,eVals,rfVals,poolingOption,numStimuli)

if ~exist('poolingOption','var')        poolingOption=2;                end
if ~exist('numStimuli','var')           numStimuli=[];                  end

opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
    'Display','off','LargeScale','off','MaxFunEvals',500);

if poolingOption==1
    addAdditive=0;
elseif poolingOption==2
    addAdditive=0;
    rfVals = rfVals-mean(mean(rfVals));
elseif poolingOption==3
    addAdditive=1;
end

[maxElevationPos,maxAzimuthPos] = find(rfVals==max(max(rfVals)));

if addAdditive
    startPt = [aVals(maxAzimuthPos(1)) eVals(maxElevationPos(1)) 0.5 0.5 0 max(max(rfVals)) mean(mean(rfVals))];
else
    startPt = [aVals(maxAzimuthPos(1)) eVals(maxElevationPos(1)) 0.5 0.5 0 max(max(rfVals))];
end

maxRFVals = max(max(rfVals));

if maxRFVals==0
    %disp('All zeros');
    params = zeros(1,length(startPt));
else
    if isempty(numStimuli)
        params = fminsearch(@(params) gauss2D(params,aVals,eVals,rfVals),startPt,opts);
    else
        params = fminsearch(@(params) gauss2D(params,aVals,eVals,rfVals,numStimuli),startPt,opts);
    end
end
end