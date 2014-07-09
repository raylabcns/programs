% gaborStim.azimuthDeg=0;
% gaborStim.elevationDeg=0;
% gaborStim.contrastPC=100;
% gaborStim.sigmaDeg=0.3;
% gaborStim.radiusDeg=3*gaborStim.sigmaDeg;
% gaborStim.spatialFreqCPD=1;
% gaborStim.orientationDeg=0;

function showGaborStimulus(gaborStim,aPoints,ePoints,plothandle)

if ~exist('plotHandle','var')       plotHandle=gca;         end
axes(plotHandle)
gaborPatch = makeGaborStimulus(gaborStim,aPoints,ePoints);
sc(gaborPatch);