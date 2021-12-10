function [ output_args ] = saveData( varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global data pltVr datVr;

tempSplit = strsplit(data.pathname,'.');

[filename, pathname, filterindex] = uiputfile('*.mat','Save processed tweezers data files',[tempSplit{1}, '.mat']);
if isequal(filename,0) || isequal(pathname,0)
    disp('User selected Cancel')
else
    
    data.pathname = fullfile(pathname,filename);
    save([data.pathname(1:end-4) '.mat'],'data', 'pltVr', 'datVr');
    
    disp(['User saved ',fullfile(pathname,filename)])
end

end

