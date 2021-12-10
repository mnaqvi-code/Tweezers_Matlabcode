function [ output_args ] = loadRuler( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global data pltVr datVr

[tFileName,tPathName,tFilterIndex] = uigetfile({'*.txt','Custom Ruler file (*.txt)'},'Pick data file',[pwd '\CustomRulers\'],'MultiSelect', 'off');

if isequal(tFileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
    if ~isempty(data)
        
    else
        data = 0;
    end
    return
else
    datVr.CustRul = dlmread(fullfile(tPathName, tFileName),'\t');
    
end

end

