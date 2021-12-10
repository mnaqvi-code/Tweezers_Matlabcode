function [ data ] = getRaw ( varargin )
%% SELECT DATA

% if nargin == 1;    
%     data = varargin{1};
% end
global data; %#ok<REDEF>
if ~isempty(data)
    [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw tweezers data files (*.dat)'},'Pick data files',strrep(data{1}.pathname, '\', filesep),'MultiSelect', 'off');
else
    if exist('accessHistory.txt', 'file') == 2
        fid = fopen( 'accessHistory.txt', 'r' );
        hisRead = textscan(fid, '%s %s %s %s %s %s %s\n', 'delimiter', '|','collectoutput',true);
        fclose(fid);
        if exist(hisRead{1}{end,5},'dir') == 7
            [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw tweezers data files (*.dat)'},'Pick data files',hisRead{1}{end,5},'MultiSelect', 'off');
        else
            [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw tweezers data files (*.dat)'},'Pick data files',pwd,'MultiSelect', 'off');
        end
    else
        [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw tweezers data files (*.dat)'},'Pick data files',pwd,'MultiSelect', 'off');
    end
    
end

if isequal(tFileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
    if exist('data')
        
    else
        data = 0;
    end
    return
else
    if exist('accessHistory.txt', 'file') == 2
        fid = fopen( 'accessHistory.txt', 'a+' );
        fprintf( fid, '%s| %s| %s| %s| %s| %s| %s\n', getenv('USERNAME'), 'selected', tFileName, 'in', tPathName, 'at', datestr(datetime('now')));
        fclose(fid);
    else
        fid = fopen( 'accessHistory.txt', 'a+' );
        fprintf( fid, '%s| %s| %s| %s| %s| %s| %s', getenv('USERNAME'), 'selected', tFileName, 'in', tPathName, 'at', datestr(datetime('now')));
        fclose(fid);
    end
    
    if ischar(tFileName)
        tData = cell(1,1);
        h = waitbar(0,'Importing Data');
        waitbar(0,h);
        disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName)])
        tData{1} = importdata(fullfile(tPathName, tFileName));
        tData{1}.pathName = fullfile(tPathName, tFileName);
        waitbar(1,h);
        close(h)
    else
        tData = cell(1,size(tFileName,2));
        h = waitbar(0,'Importing Data');
        for ii = 1:size(tFileName,2)
            waitbar((ii-1)/size(tFileName,2),h);
            disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName{ii})])
            tData{ii} = importdata(fullfile(tPathName, tFileName{ii}));
            tData{ii}.pathName = fullfile(tPathName, tFileName{ii});
        end
        waitbar(1,h);
        close(h)
        clear ii
    end
end

%
%% SELECT COMMENT FILE

[cFileName,cPathName,cFilterIndex] = uigetfile({'*.txt','Comment data file (*.dat)'},'Pick comment file',tPathName,'MultiSelect', 'off');

if isequal(cFileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
    return
else
    cData = cell(1,1);
    h = waitbar(0,'Loading Comment File');
    waitbar(0,h);
    disp([getenv('USERNAME'), ' selected ', fullfile(cPathName, cFileName)])
    %cData{1}.data = importdata(fullfile(cPathName, cFileName));
    %keyboard
    cData{1}.pathName = fullfile(cPathName, cFileName);
    waitbar(1,h);
    close(h)
end
%% FORMAT ALL DATA
clear data
h = waitbar(0,'Formatting Data');
waitbar(0,h);
[file,time,comment] = textread(cData{1}.pathName,'%s %d %s','delimiter','\t');	%reads comment file to the variables file, time, comment using format mask
filenumber = strtok(tFileName,'.');
dayfiles = strmatch(date,file);
trapcommentsindex = strmatch('**k',comment);	%gets comments' indices with values of trap constants

% Get the data from the datafiles and put them in arrays

if not(iscell(tFileName))
    tFileName = {tFileName};  %fileName becomes cell array
end
if not(iscell(filenumber))
    filenumber = {filenumber};  %filenumber becomes cell array
end

for ii = 1:size(tFileName,2)
    waitbar((ii-1)/size(tFileName,2),h);
    if isempty(strmatch(filenumber(ii),file))
        commentindex{ii} = 1;
    else
        commentindex{ii} = strmatch(filenumber(ii),file);   %index of comments written during the saving of i-file !!gives a mistake if there are no commenta written during saving this file
    end
    
    data{ii}.commentindex = commentindex{ii};   %creation of cell array with comments related to i-file
    data{ii}.commenttime = time(commentindex{ii});
    data{ii}.commenttext = comment(commentindex{ii});
    data{ii}.filename = tFileName{ii};
    data{ii}.pathname = tPathName;
    data{ii}.date = date;
    data{ii}.time = tData{ii}.data(:,1);
    for jj = 2:size(tData{ii}.textdata,2)
        
        data{ii}.(genvarname(tData{ii}.textdata{jj})) = tData{ii}.data(:,jj);
        
    end
    % get trap constants from comment file and put them in struct.
    kindex = trapcommentsindex(find(trapcommentsindex < commentindex{ii}(1),1,'last'));     % looks for last trap constant before this file
    kindex = [kindex;trapcommentsindex(find(trapcommentsindex >= commentindex{ii}(1),3,'first'))];      % looks for three first trap constants during/after this file
    % output things to get right trapping parameters
    constantsstring = comment(kindex);
    H1 = figure(200);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    pause(0.1);
    
    colors      =   linspecer(8);
    set(gcf,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name','Choose Calibration Values');
    set(gcf,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
        'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
    constants = 0;
    for jj = 1:1:size(constantsstring,1)
        
        uicontrol('Style','Pushbutton', ...
            'Position',[0.1 0.7-(jj-1)*(0.2) 0.8 0.1], ...
            'String',constantsstring{jj},...
            'FontSize',0.5,...
            'BackgroundColor',colors(2,:),...
            'ForegroundColor',colors(5,:),...
            'Callback',{@delete_graphic,jj,constantsstring});
    end
    waitfor(H1);
    figure(1);
    constants = getappdata(0, 'constants');
    data{ii}.kx = constants(1);
    data{ii}.ky = constants(2);
    data{ii}.dvdx = constants(3);
    data{ii}.dvdy = constants(4);
    clear kindex constantsstring nconstants constants
    
    % calculate forces from trapping parameters and voltages
    data{ii}.fx = data{ii}.vx*data{ii}.kx/data{ii}.dvdx;
    data{ii}.fx0 = 0;
    data{ii}.fy = data{ii}.vy*data{ii}.ky/data{ii}.dvdy;
    data{ii}.fy0 = 0;
    data{ii}.ext = -(data{ii}.xpz+data{ii}.vx/data{ii}.dvdx);
    data{ii}.ext0 = 0;
    data{ii}.shiftX = 0;
    data{ii}.shiftY = 0;
    
end
waitbar(1,h);
close(h)

clear cData cFileName cFilterIndex cPathName comment commentindex dayfiles file filenumber h ii jj tData tFileName tFilterIndex tPathName time trapcommentsindex
end

function delete_graphic(hObject,eventdata,A,constantsstring)
    constants = sscanf(constantsstring{A}(4:end), '%f %f %f %f');
    setappdata(0, 'constants', constants);
    close('Choose Calibration Values');
end