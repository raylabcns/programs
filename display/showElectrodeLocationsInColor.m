function [plotHandle,colorNames,electrodeArray,colorRows,colorCols] = showElectrodeLocationsInColor(gridPosition,plotHandle,holdOnState,hideElectrodeNums,startRow,gridType)

if ~exist('hideElectrodeNums','var')    hideElectrodeNums=0;            end
if ~exist('startRow','var')             startRow=5;                     end
if ~exist('gridType','var')             gridType = 'Microelectrode';    end

if strcmpi(gridType,'ECoG')
    numRows=8;numCols=10;
else 
    numRows=10;numCols=10;
end
[~,~,electrodeArray] = electrodePositionOnGrid(1,gridType);

if ~exist('plotHandle','var') || isempty(plotHandle)
    plotHandle = subplot('Position',gridPosition,'XTickLabel',[],'YTickLabel',[],'box','on');
end
if ~exist('holdOnState','var')
    holdOnState = 1;
end

if ~holdOnState
    cla(plotHandle);
end

axes(plotHandle) %#ok<MAXES>
dX = 1/numCols;
dY = 1/numRows;

lineXRow = zeros(2,numRows);lineYRow = zeros(2,numRows);
for i=1:numRows
    lineXRow(:,i) = [0 1]; lineYRow(:,i) = [i*dY i*dY];
end
lineXCol = zeros(2,numCols);lineYCol = zeros(2,numCols);
for i=1:numCols
    lineXCol(:,i) = [i*dX i*dX]; lineYCol(:,i) = [0 1];
end
line(lineXRow,lineYRow,'color','k'); hold on
line(lineXCol,lineYCol,'color','k'); 
hold off;

% Color the patches - choose the subgrid to color
colorRows = startRow:numRows;
colorCols = 1:numCols;

% The four corners have the following colors

% yellow (1 1 0)                      % Green (0 1 0)

% red (1 0 0)                         % Blue (0 0 1)

numColorRows = length(colorRows);
numColorCols = length(colorCols);

colorNames = cell(numColorRows,numColorCols);
for j=1:numColorRows
    for i=1:numColorCols
        
        row = colorRows(j);
        col = colorCols(i);
        
        electrodeNum = electrodeArray(row,col);
        
        if electrodeNum>0
            % red value does not depend on row, goes from 1 to 0 with col
            if numColorCols==1
                rValue = 1-i/numColorCols;
            else
                rValue = 1-(i-1)/(numColorCols-1);
            end
            
            % green value goes from 1 to 0 with row, does not depend on col
            if numColorRows==1
                gValue = 1-j/numColorRows;
            else
                gValue = 1-(j-1)/(numColorRows-1);
            end
            
            % blue value is 0 on the three corners and 1 at the bottom right
            % corner. We approximate this by the following
            if numColorRows==1
                bNum1 = j; bDen1 = 1;
            else
                bNum1 = j-1; bDen1 = numColorRows-1;
            end
            
            if numColorCols==1
                bNum2 = i; bDen2 = 1;
            else
                bNum2 = i-1; bDen2 = numColorCols-1;
            end
            
            bValue = sqrt((bNum1*bNum2)/(bDen1*bDen2));
            colorName = [rValue gValue bValue];
        else
            colorName = [1 1 1];
        end
            
        highlightRow=row;
        highlightCol=col;

        % Create patch
        patchX = (highlightCol-1)*dX;
        patchY = (numRows-highlightRow)*dY;
        patchLocX = [patchX patchX patchX+dX patchX+dX];
        patchLocY = [patchY patchY+dY patchY+dY patchY];
        
        patch(patchLocX,patchLocY,colorName);
        colorNames{row,col} = colorName;
    end
end

if ~hideElectrodeNums
    % Write electrode numbers
    for i=1:numRows
        textY = (numRows-i)*dY + dY/2;
        for j=1:numCols
            textX = (j-1)*dX + dX/2;
            if electrodeArray(i,j)>0
                text(textX,textY,num2str(electrodeArray(i,j)),'HorizontalAlignment','center');
            end
        end
    end
end
end