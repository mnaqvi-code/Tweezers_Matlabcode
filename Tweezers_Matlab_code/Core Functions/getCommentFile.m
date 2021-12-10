function [ cFileName,cPathName,cFilterIndex,cData ] = getCommentFile( tPathName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global data;

[cFileName,cPathName,cFilterIndex] = uigetfile({'*.txt','Comment data file (*.dat)'},'Pick comment file',tPathName,'MultiSelect', 'off');

if isequal(cFileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
    
    cFileName = [];
    cPathName = [];
    cFilterIndex = [];
    cData = [];    
    return
else
    cData = [];
    h = waitbar(0,'Loading Comment File');
    waitbar(0,h);
    disp([getenv('USERNAME'), ' selected ', fullfile(cPathName, cFileName)])
    %cData{1}.data = importdata(fullfile(cPathName, cFileName));
    %keyboard
    cData.pathname = fullfile(cPathName, cFileName);
    waitbar(1,h);
    close(h)
end

end

