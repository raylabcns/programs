function plotGaborBoundaryOnSCImage(aPoints,ePoints,boundaryX,boundaryY,colorName,plotHandle)

newBoundaryX = 1+length(aPoints)*(boundaryX-aPoints(1))/(aPoints(end)-aPoints(1));
newBoundaryY = length(ePoints) - length(ePoints)*(boundaryY-ePoints(1))/(ePoints(end)-ePoints(1)) + 1;
set(plotHandle,'NextPlot','add','box','on');
plot(plotHandle,newBoundaryX,newBoundaryY,'color',colorName,'linewidth',2);
%plot(plotHandle,newBoundaryX,newBoundaryY,'color','r','linewidth',2);