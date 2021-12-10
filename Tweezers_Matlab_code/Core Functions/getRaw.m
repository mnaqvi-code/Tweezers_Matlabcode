function [ data ] = getRaw ( varargin )
%% SELECT DATA

% if nargin == 1;    
%     data = varargin{1};
% end
global data binaryD; %#ok<REDEF>

if exist('data') %#ok<*EXIST>
    
else
    data = []; %#ok<*NASGU>
end

if ~isempty(data) %#ok<*NODEF>
    [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw tweezers data files (*.dat)'},'Pick data files',strrep(data.pathname, '\', filesep),'MultiSelect', 'off'); %#ok<*ASGLU>
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
    if exist('data') %#ok<*EXIST>
        
    else
        data = 0; %#ok<*NASGU>
    end
    return
else
    if exist('accessHistory.txt', 'file') == 2
        fid = fopen( 'accessHistory.txt', 'a+' );
        fprintf( fid, '%s| %s| %s| %s| %s| %s| %s\n', getenv('USERNAME'), 'selected', tFileName, 'in', tPathName, 'at', datestr(clock, 0));%datestr(datetime('now'))
        fclose(fid);
    else
        fid = fopen( 'accessHistory.txt', 'a+' );
        fprintf( fid, '%s| %s| %s| %s| %s| %s| %s', getenv('USERNAME'), 'selected', tFileName, 'in', tPathName, 'at', datestr(clock, 0));
        fclose(fid);
    end
    
    if ischar(tFileName)
        tData = [];
        h = waitbar(0,'Importing Data');
        waitbar(0,h);
        disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName)])
        if binaryD == 0
            tData = importdata(fullfile(tPathName, tFileName));
            tData.pathname = fullfile(tPathName, tFileName);
            [pathstr,name,ext]=fileparts(fullfile(tPathName, tFileName));
            FileNameShort=sprintf('%s%s',name,ext);
            FileNameNoExt=name;
            FileFolder=pathstr;
            tData.FileName=FileNameShort;
            tData.FileFolder=FileFolder;
            
            waitbar(1,h);
            close(h)
        else
            clear g     
            g = getBinary(fullfile(tPathName, tFileName));
            if isfield(data, 'x')
                data = rmfield(data,'x');
            end
            
            [pathstr,name,ext]=fileparts(fullfile(tPathName, tFileName));
            FileNameShort=sprintf('%s%s',name,ext);
            FileNameNoExt=name;
            FileFolder=pathstr;
            data.FileName=FileNameShort;
            data.FileFolder=FileFolder;            
            data.pathname = fullfile(tPathName, tFileName);
            
            data.filename = tFileName; %#ok<*STRNU>
            if g.SampleFrequency > 1000
                data.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency/1000) ' kHz']};
            else
                data.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency) 'Hz']};
            end
            
            %f.DmaData.Length/13
            data.textdata{3,1} = {'Timestamp'};
            data.data = g.Data;
            if isfield(g, 'beadData')
                data.beadData = g.beadData;                
            end
            data.data(:,1) = data.data(:,1) - data.data(1,1);
    %             for ii = 2:12
    %                 data.data(:,ii) = data.data(:,ii+1);
    %             end
            data.textdata{3,2} = {'FrameNumber'};
            data.textdata{3,3} = {'Psd1_VxDiff'};
            data.textdata{3,4} = {'Psd1_VxSum'};
            data.textdata{3,5} = {'Psd1_VyDiff'};
            data.textdata{3,6} = {'Psd1_VySum'};
            data.textdata{3,7} = {'Psd2_VxDiff'};
            data.textdata{3,8} = {'Psd2_VxSum'};
            data.textdata{3,9} = {'Psd2_VyDiff'};
            data.textdata{3,10} = {'Psd2_VySum'};
            data.textdata{3,11} = {'MirrorX'};
            data.textdata{3,12} = {'MirrorY'};            
            data.textdata{3,13} = {'Status'};
            data.textdata{3,14} = {'Camera Time'};
            data.textdata{3,15} = {'Camera Frame'}; 
            data.textdata{3,16} = {'Bead1X'};
            data.textdata{3,17} = {'Bead1Y'};            
            data.textdata{3,18} = {'Bead2X'};            
            data.textdata{3,19} = {'Bead2Y'};
            waitbar(1,h);
            close(h)
            
            idx = max(strfind(data.pathname,'\'));
            my_path = data.pathname(1:idx);
            my_file_name = data.pathname(idx+1:end);
            data.FileInfo = dir(data.pathname);     
            clear idx
            
            if exist('DualTrapCalValues.txt', 'file') == 2
                
                writeValAppend = dlmread('DualTrapCalValues.txt');
                timeMin = abs(writeValAppend(:,1)-data.FileInfo.datenum);
                [idx idx] = min(timeMin);
                keyboard;
                datVr.dualCst.Dex1 =  writeValAppend(idx,2);
                datVr.dualCst.cFreq1 = writeValAppend(idx,3); %Hz
                datVr.dualCst.Dex2 =  writeValAppend(idx,4);
                datVr.dualCst.cFreq2 = writeValAppend(idx,5); %Hz
                datVr.dualCst.Dex3 =  writeValAppend(idx,6);
                datVr.dualCst.cFreq3 = writeValAppend(idx,7); %Hz
                datVr.dualCst.Dex4 =  writeValAppend(idx,8);
                datVr.dualCst.cFreq4 = writeValAppend(idx,9); %Hz
            
            end
            
            return
        end
   
    else
        tData = cell(1,size(tFileName,2)); %#ok<*PREALL>
        h = waitbar(0,'Importing Data');

        disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName{ii})])
        tData = importdata(fullfile(tPathName, tFileName{ii}));
        [pathstr,name,ext]=fileparts(fullfile(tPathName, tFileName));
        FileNameShort=sprintf('%s%s',name,ext);
        FileNameNoExt=name;
        FileFolder=pathstr;
        tData.FileName=FileNameShort;
        tData.FileFolder=FileFolder;
        tData.pathname = fullfile(tPathName, tFileName{ii});

        waitbar(1,h);
        close(h)
        clear ii
    end
end
%% SELECT COMMENT FILE

[ cFileName,cPathName,cFilterIndex,cData ] = getCommentFile( tPathName );
if isempty(cFileName)  
    return
end

%% FORMAT ALL DATA
clear data
h = waitbar(0,'Formatting Data');
waitbar(0,h);
[file,time,comment] = textread(cData.pathname,'%s %d %s','delimiter','\t');	%reads comment file to the variables file, time, comment using format mask
filenumber = strtok(tFileName,'.');
dayfiles = strmatch(date,file);
trapcommentsindex = strmatch('**k',comment);	%gets comments' indices with values of trap constants
% Get the data from the datafiles and put them in arrays

% if not(iscell(tFileName))
%     tFileName = {tFileName};  %fileName becomes cell array
% end
% if not(iscell(filenumber))
%     filenumber = {filenumber};  %filenumber becomes cell array
% end



if isempty(strmatch(filenumber,file))
    commentindex = 1;
else
    commentindex = strmatch(filenumber,file);   %index of comments written during the saving of i-file !!gives a mistake if there are no commenta written during saving this file
end

data.commentindex = commentindex;   %creation of cell array with comments related to i-file
data.commenttime = time(commentindex);
data.commenttext = comment(commentindex);
data.filename = tFileName;

[pathstr,name,ext]=fileparts(fullfile(tPathName, tFileName));
FileNameShort=sprintf('%s%s',name,ext);
FileNameNoExt=name;
FileFolder=pathstr;
data.FileName=FileNameShort;
data.FileFolder=FileFolder;
data.pathname = tPathName;
data.date = date;
data.time = tData.data(:,1);
for jj = 2:size(tData.textdata,2)
    if exist('genvarname') == 2
        data.(genvarname(tData.textdata{jj})) = tData.data(:,jj);
    else
        data.(matlab.lang.makeValidName(tData.textdata{jj})) = tData.data(:,jj);
    end
end
% get trap constants from comment file and put them in struct.
kindex = trapcommentsindex(find(trapcommentsindex < commentindex(1),1,'last'));     % looks for last trap constant before this file
kindex = [kindex;trapcommentsindex(find(trapcommentsindex >= commentindex(1),3,'first'))];      % looks for three first trap constants during/after this file
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
data.kx = constants(1);
data.ky = constants(2);
data.dvdx = constants(3);
data.dvdy = constants(4);
clear kindex constantsstring nconstants constants

% calculate forces from trapping parameters and voltages
data.fx = data.vx*data.kx/data.dvdx;
data.fx0 = 0;
data.fy = data.vy*data.ky/data.dvdy;
data.fy0 = 0;
data.ext = -(data.xpz+data.vx/data.dvdx);
data.ext0 = 0;
data.shiftX = 0;
data.shiftY = 0;
waitbar(1,h);  
close(h) 

clear cData cFileName cFilterIndex cPathName comment commentindex dayfiles file filenumber h ii jj tData tFileName tFilterIndex tPathName time trapcommentsindex
end

function delete_graphic(hObject,eventdata,A,constantsstring)
    constants = sscanf(constantsstring{A}(4:end), '%f %f %f %f');
    setappdata(0, 'constants', constants);
    close('Choose Calibration Values');
end