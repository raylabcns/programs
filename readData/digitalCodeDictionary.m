% digitalCodeDictionary(str);

% Digital Code dictionary
% This program maintains a list of all the codes that are put in the digital data
% stream by various Lablib programs.

function codeExistsFlag = digitalCodeDictionary(str)

if ~exist('str','var')    str='';     end

codeList = [...
    'AL';
    'AZ';
    'BR'; 
    'CO';
    'CT';
    'EC';
    'EL';
    'FI';
    'FO'; 
    'IT'; 
    'M0'; 
    'M1'; 
    'ON'; 
    'OF'; 
    'OR'; 
    'PA'; 
    'RA'; 
    'SA'; 
    'SF'; 
    'SI'; 
    'ST'; 
    'TC'; 
    'TE'; 
    'TF'; 
    'TG'; 
    'TS'; 
    'T0'; 
    'T1'; 
];

if isempty(str)
    disp(codeList);
    codeExistsFlag=1;
else
    codeExistsFlag=0;
    for i=1:size(codeList,1)
        if isequal(str,codeList(i,:))
            codeExistsFlag=1;
            break;
        end
    end
end