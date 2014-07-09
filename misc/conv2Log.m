function Z = conv2Log(Y,Thres)

if nargin==1
    Thres = 10^(-16);
end

%%% Replace zeros with Thres or min(Y), whichever is smaller;
%minY = min(min(Y(find(Y>0))),Thres);
%Y(find(Y==0)) = minY;
%Z=log10(Y);

%%%% Replace anything less than Thres with Thres
Y(find(Y<Thres)) = Thres;
Z=log10(Y);

return