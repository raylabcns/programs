% H =getPSTH(X,d,Tmin,Tmax,smoothSigmaMs)
% This function returns the psth of the cell array X calculated with bin width of d ms. Default = 1 ms
% We assume that each cell in X has a sequence of spike times in sec. 
% Tmin and Tmax denote the range for which the PSTH is computed (in
% seconds)

% Supratim Ray 05/21/05
% getPSTH is used while loading Chronux because it also has a psth
% function which overrides this one

% Including the smoothing function.
% 11/6/14: Name changed from psth_SR to getPSTH

function [H,timeVals] = getPSTH(X,d,Tmin,Tmax,smoothSigmaMs)

if ~exist('smoothSigmaMs','var')      smoothSigmaMs=[];                     end

if nargin==1
    d = 0.001; % 1ms
else
    d = d/1000; %converting in s
end

spk = [];
numTrials = length(X);

for i=1:numTrials
    [A B] = size(X{i});
    if A == 1
        spk = cat(2,spk,X{i}); % row vector
    elseif B == 1
        spk = cat(2,spk,X{i}');
    end
end

N = (Tmax-Tmin)/d;  % number of bins

for i=1:N
    Htmp(i) = length(find(spk < Tmin+i*d)); %#ok<*AGROW>
    timeVals(i) = Tmin+i*d-d/2;
end

Htmp0 = length(find(spk<Tmin));
H = cat(2,Htmp(1)-Htmp0,diff(Htmp));

% Scale
% H contains the number of spikes in M trials and d second bin
H = (H/numTrials) /d; % spikes per trial per second

if ~isempty(smoothSigmaMs)
    H = gaussSmooth(H,(smoothSigmaMs/1000)/d);
end

end

function smoothData = gaussSmooth(data,sigma,windowLen)

if ~exist('sigma','var')            sigma = 1;                          end
if ~exist('windowLen','var')        windowLen = 11;                     end

window = exp(-0.5*(((1:windowLen) - (windowLen+1)/2)/sigma).^2);
window = window/sum(window);
smoothData = convn(data,window,'same');

end