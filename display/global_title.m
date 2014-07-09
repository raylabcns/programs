function h_out = global_title(X,Y,str,backGround,FONTSIZE,horizAlign,fig)
% global_title(str)
%
% INPUT: X,Y        - position of the title [0 1]
%        str        - Title_string
%        backGround - backGroundColor
%        FONTSIZE   - fontSize of the title (default 12)
%        horizAlign - HorizontalAlignment (default - center)
%        fig        - the handle of the figure of interest (default -> current figure)

% Peter Denchev - 2004/02/16
% modified by Arun Sripati

if(~exist('backGround')|isempty(backGround)) backGround = [0.75 0.75 0.75]; end; 
if(~exist('FONTSIZE')) FONTSIZE = 12; end; 
if(~exist('horizAlign')) horizAlign = 'center'; end; 
if(~exist('fig')) fig = gcf; end; 

cax    = gca; set(cax,'Units','normalized');
pos    = get(cax,'position');

xa0    = pos(1); ya0    = pos(2);
deltaX = pos(3); deltaY = pos(4);

XLim   = get(cax,'XLim'); YLim   = get(cax,'YLim');
Xscale = (XLim(2)-XLim(1))/deltaX; Yscale = (YLim(2)-YLim(1))/deltaY;

x_text = XLim(1) + (X-xa0)*Xscale; y_text = YLim(1) + (Y-ya0)*Yscale;

figure(fig);
h_out = text(x_text,y_text,str,'fontSize',FONTSIZE,'BackgroundColor',backGround,'HorizontalAlignment',horizAlign);

return