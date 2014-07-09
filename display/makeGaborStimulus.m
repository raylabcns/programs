% makeGaborStimulus(gaborStim,aVals,eVals)
% G(x,y) = exp(-((x'^2) + (y'^2))/2) sin(2*pi*f*x')
% where 
% x' = xcos(theta)-ysin(theta)
% y' = xsin(theta)+ycos(theta)

function [gaborPatch,aperature] = makeGaborStimulus(gaborStim,aVals,eVals,showGabor)

if ~exist('showGabor','var')        showGabor=0;            end

azi = gaborStim.azimuthDeg;
ele = gaborStim.elevationDeg;
sf  = gaborStim.spatialFreqCPD;
ori = gaborStim.orientationDeg;
C   = gaborStim.contrastPC/2;

if length(gaborStim.radiusDeg)==1    % Gabor
    radMax = gaborStim.radiusDeg;
    radMin = 0;
else
    radMax = gaborStim.radiusDeg(2);
    radMin = gaborStim.radiusDeg(1);
end


% Make a 2D grating
theta = pi*ori/180; % converting to radians

for e=1:length(eVals)  % 
    for a=1:length(aVals)
        xg = aVals(a)*cos(theta) - eVals(e)*sin(theta);
        grating(e,a) = C*sin(2*pi*sf*xg);
        
        distance = sqrt((aVals(a)-azi)^2+(eVals(e)-ele)^2);
        
        if (distance<=radMax) && (distance>=radMin)
            aperature(e,a) = 1;
        else
            aperature(e,a) = 0;
        end
    end
end

% Gaussian
params(1) = azi;
params(2) = ele;
params(3) = gaborStim.sigmaDeg;
params(4) = params(3);
params(5) = 0;
params(6) = 1;
[diffOut,GaussianEnvelope,boundaryX,boundaryY] = gauss2D(params,aVals,eVals,[]);

% set everything outside radius to zero
gaborPatch = 50+GaussianEnvelope.*aperature.*grating;

if showGabor
    colormap('gray');
    subplot(221);
    pcolor(aVals,eVals,grating); shading interp; colorbar;
    subplot(222);
    pcolor(aVals,eVals,GaussianEnvelope); shading interp; colorbar;
    subplot(223);
    pcolor(aVals,eVals,aperature); shading interp; colorbar;
    subplot(224);
    pcolor(aVals,eVals,gaborPatch); shading interp; colorbar
    hold on;
    plot(boundaryX,boundaryY,'r');
end