function decList = convertStrCodeToDec(codeInStrList,useSingelITC18Flag)

if ~exist('useSingelITC18Flag','var')       useSingelITC18Flag=1;       end

for i=1:size(codeInStrList,1)
    
    if useSingelITC18Flag       %8-14 represent the first letter, 1-7 represent the second.
        decList(i)=256*double(codeInStrList(i,1))+2*double(codeInStrList(i,2)); %#ok<*AGROW>
    else                        %8 bits are used for a hex letter
        decList(i)=256*double(codeInStrList(i,1))+double(codeInStrList(i,2));
    end
end
end