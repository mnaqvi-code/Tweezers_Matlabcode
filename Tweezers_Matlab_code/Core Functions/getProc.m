function [ output_args ] = getProc( varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% if nargin == 1;
%     data = varargin{1};
% end
global data pltVr datVr;
if ~isempty(data)
    [tFileName,tPathName,tFilterIndex] = uigetfile({'*.mat','Processed tweezers data files (*.mat)'},'Pick data files',strrep(data.pathname, '\', filesep),'MultiSelect', 'off');
else
    [tFileName,tPathName,tFilterIndex] = uigetfile({'*.mat','Processed tweezers data files (*.mat)'},'Pick data files',pwd,'MultiSelect', 'off');
end


if isequal(tFileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
    if ~isempty(data)
        
    else
        data = 0;
    end
    return
else
    data = [];
    if ischar(tFileName)
        h = waitbar(0,'Importing Data');
        waitbar(0,h);
        disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName)])
        load(fullfile(tPathName, tFileName),'data', 'pltVr', 'datVr');
        if iscell(data)
            data = data{1};
        end
        
        if ~isfield(data,'shiftX')
            data.shiftX = 0;
            data.shiftY = 0;
        end
        waitbar(1,h);
        close(h)
    else
        tData = cell(1,size(tFileName,2));
        h = waitbar(0,'Importing Data');
        for ii = 1:size(tFileName,2)
            waitbar((ii-1)/size(tFileName,2),h);
            disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName{ii})])
            datatemp = load(fullfile(tPathName, tFileName{ii}),'data');
            data = datatemp.data;
            if iscell(data)
                data = data{1};
            end
            if ~isfield(data,'shiftX')
                data.shiftX = 0;
                data.shiftY = 0;
            end
        end
        waitbar(1,h);
        close(h)
        clear ii
    end
end
end

