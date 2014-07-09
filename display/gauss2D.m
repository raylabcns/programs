% gauss2D generates a 2D gaussian with 7 parameters (x0,y0,a,b,theta,C,D)
% x0,y0 - center of the gaussian
% a, b  - standard deviation of the Gaussian along magor and minor axes
% theta - angle of rotation with respect to the coordinate system
% C     - scaling factor
% D     - additive constant

% xg = (x-x0)cos(theta) - (y-y0)sin(theta)
% yg = (x-x0)sin(theta) + (y-y0)cos(theta)
% W(x,y) = C*exp(-1/2(xg^2/a^2+yg^2/b^2))+D

% aVals, eVals - the x and y positions where the 2D Gaussian should be
% computed.

% The additive constant D is ignored if length(params) == 6

% Added optional parameter numStimuli
% This scales the computation of diffOut

function [diffOut,outVals,boundaryX,boundaryY] = gauss2D(params,aVals,eVals,rfVals,numStimuli)

if ~exist('rfVals','var')           rfVals=[];                      end
if ~exist('aVals','var')            aVals=[];                       end
if ~exist('numStimuli','var')       numStimuli=[];                  end

x0 = params(1); y0 = params(2);
sx  = params(3); sy  = params(4);
theta  = params(5);
C      = params(6);

if length(params)>6
    D = params(7);
else
    D = 0;
end

if isempty(aVals)
    outVals=[];
else
    for e=1:length(eVals)
        for a=1:length(aVals)
            xg = (aVals(a)-x0)*cos(theta) - (eVals(e)-y0)*sin(theta);
            yg = (aVals(a)-x0)*sin(theta) + (eVals(e)-y0)*cos(theta);

            outVals(e,a) = C*exp(-(1/2)*(xg^2/sx^2+yg^2/sy^2))+D;
        end
    end
    
end

% Generate the "Boundary", which are the (i,j) values at which
% outVals(i,j)=exp(-1/2);


xVals = -sx:sx/100:sx;
yVals = sy*sqrt(1-(xVals/sx).^2);

boundaryX1 = x0 + xVals*cos(theta) + yVals*sin(theta);
boundaryY1 = y0 - xVals*sin(theta) + yVals*cos(theta);

xVals = sx:-sx/100:-sx;
yVals = -sy*sqrt(1-(xVals/sx).^2);

boundaryX2 = x0 + xVals*cos(theta) + yVals*sin(theta);
boundaryY2 = y0 - xVals*sin(theta) + yVals*cos(theta);

boundaryX = [boundaryX1 boundaryX2];
boundaryY = [boundaryY1 boundaryY2];


if isempty(rfVals)
    diffOut=[];
else
    if isempty(numStimuli)
        diffOut = sum(sum((outVals-rfVals).^2));
    else
        diffOut = sum(sum((numStimuli .* (outVals-rfVals).^2)));
    end
end
end