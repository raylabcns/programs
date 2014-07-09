function [eyeDataDegX,eyeDataDegY] = convertEyeDataToDeg(eyeData,returnCellArray)

if ~exist('returnCellArray','var')      returnCellArray=0;              end

N = length(eyeData);
if returnCellArray
    eyeDataDegX = cell(1,N);
    eyeDataDegY = cell(1,N);
else
    eyeDataDegX = zeros(N,length(eyeData(1).eyePosDataX));
    eyeDataDegY = zeros(N,length(eyeData(1).eyePosDataY));
end

for i=1:N   
    clear cal eyeX eyeY
    cal = eyeData(i).eyeCal;
    eyeX = eyeData(i).eyePosDataX;
    eyeY = eyeData(i).eyePosDataY; 
    
    if returnCellArray
        eyeDataDegX{i} = cal.m11*eyeX + cal.m21*eyeY + cal.tX;
        eyeDataDegY{i} = cal.m12*eyeX + cal.m22*eyeY + cal.tY;
    else
        eyeDataDegX(i,:) = cal.m11*eyeX + cal.m21*eyeY + cal.tX;
        eyeDataDegY(i,:) = cal.m12*eyeX + cal.m22*eyeY + cal.tY;
    end
end