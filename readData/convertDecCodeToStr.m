function codeInStrList = convertDecCodeToStr(codeInDecList,useSingelITC18Flag)

if ~exist('useSingelITC18Flag','var')       useSingelITC18Flag=1;       end

for i=1:length(codeInDecList)
    
    if useSingelITC18Flag       %8-14 represent the first letter, 1-7 represent the second.
        d=dec2bin(codeInDecList(i),16);
        codeInStrList(i,:) = [char(bin2dec(d(2:8))) char(bin2dec(d(9:15)))]; %#ok<*AGROW>
    else                        %8 bits are used for a hex letter
        y=dec2hex(codeInDecList(i));
        codeInStrList(i,:) = [char(hex2dec(y(1:2))) char(hex2dec(y(3:4)))];
    end
end
end