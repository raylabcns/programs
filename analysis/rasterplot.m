% rasterplot(X) displays the raster. 

% The cell array X contains the spike times.
% trialNums: optional variable that denotes the trial number of X -
% default: 1:length(X)
% colorName - color of the raster
% d

function rasterplot(X,trialNums,colorName,d)

if ~exist('d','var')                d=1;                                end
if ~exist('trialNums','var')        trialNums=1:length(X);              end
if ~exist('colorName','var')        colorName='k';                      end

for i = 1:length(X)
    if ~isempty(X{i})
        rplot(X{i},trialNums(i),d,colorName); hold on;
    end
end
end

function hh = rplot(x, y, l,u,symbol,colorname)

if nargin < 6
    colorname = '';
end

%RPLOT = Modification of ERRORBAR for the needs of RASTER
%
%   Plots vertical lines for rasters 
%   INPUT: x - vector with absciae (times of the events)
%          y - scalar defining the position of the raster trait on the ordinate  
%          l - l(i)=const=length of the vertical bar 
%          u,symbol - not necessary arguments
%   EXAMPLE: y=0.65; e=0.01; rplot(x,y,e);
%
%   L. Shure 5-17-88, 10-1-91 B.A. Jones 4-5-93
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.18 $             $Date: 2001/04/15 12:03:51 $
%   $Modified  by Peter Denchev   $Date: 2002/08/07 
%   $Modified  by Supratim Ray: color option added  $Date: 2005/08/08
szx=max(size(x));
y=y(1);
y=ones(1,szx)*y;
l=l(1)/2;
l=ones(1,szx)*l;

if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
    if nargin > 2,
        if ~isstr(l),  
            l = l(:);
        end
        if nargin > 3
            if ~isstr(u)
                u = u(:);
            end
        end
    end
else
  [npt,n] = size(x);
end

if nargin == 3
    if ~isstr(l)  
        u = l;
        symbol = '-';
    else
        symbol = l;
        l = y;
        u = y;
        y = x;
        [m,n] = size(y);
        x(:) = (1:npt)'*ones(1,n);;
    end
end

if nargin == 4
    if isstr(u),    
        symbol = u;
        u = l;
    else
        symbol = '-';
    end
end


if nargin == 2
    l = y;
    u = y;
    y = x;
    [m,n] = size(y);
    x(:) = (1:npt)'*ones(1,n);;
    symbol = '-';
end

u = abs(u);
l = abs(l);
    
if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
  szx=size(x)
  szy=size(y)
  szl=size(l)
  szu=size(u)
  error('The sizes of X, Y, L and U must be the same.');
end

tee = 0; 
xl = x - tee;
xr = x + tee;
ytop = y + u;
ybot = y - l;
n = size(y,2);

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
xb = zeros(npt*9,n);
xb(1:9:end,:) = x;
xb(2:9:end,:) = x;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xl;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = ytop;
yb(5:9:end,:) = ytop;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ybot;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;

[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end

% adding the color option
if isempty(colorname)
    % do nothing
else
    clear col;
    col = colorname;
end

symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol); hold on

if ~hold_state, hold off; end

if nargout>0, hh = h; end
end