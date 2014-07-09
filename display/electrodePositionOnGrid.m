function [row,column,electrodeArray] = electrodePositionOnGrid(electrodeNum,gridType)

if strcmp(gridType,'Microelectrode')

    electrodeArray = ...
       [00 02 01 03 04 06 08 10 14 00;
        65 66 33 34 07 09 11 12 16 18;
        67 68 35 36 05 17 13 23 20 22;
        69 70 37 38 48 15 19 25 27 24;
        71 72 39 40 42 50 54 21 29 26;
        73 74 41 43 44 46 52 62 31 28;
        75 76 45 47 51 56 58 60 64 30;
        77 78 82 49 53 55 57 59 61 32;
        79 80 84 86 87 89 91 94 63 95;
        00 81 83 85 88 90 92 93 96 00];
      
    if electrodeNum<1 || electrodeNum>96
        disp('Electrode Number out of range');
        row=1;column=1;
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end

elseif strcmp(gridType,'ECoG')
    
    goodElectrodeList = [1:34 65:96];
    electrodeArray = ...
       [00 89 81 73 65 25 17 09 01 00;
        00 90 82 74 66 26 18 10 02 00;
        00 91 83 75 67 27 19 11 03 00;
        34 92 84 76 68 28 20 12 04 33;
        00 93 85 77 69 29 21 13 05 00;
        00 94 86 78 70 30 22 14 06 00;
        00 95 87 79 71 31 23 15 07 00;
        00 96 88 80 72 32 24 16 08 00];
    
    if isempty(setdiff(goodElectrodeList,electrodeNum))
        disp('Electrode Number out of range');
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end
end
end
    