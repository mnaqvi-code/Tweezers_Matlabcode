function [ data ] = getRaw ( varargin )
%% SELECT DATA

if nargin == 1;
    binary = varargin{1};    
elseif nargin == 2;
    binary = varargin{1};
    data = varargin{2};    
end

if exist('data')
    
    [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw Foldometer data files (*.dat)'},'Pick data files',strrep(data{1}.pathName, '\', filesep),'MultiSelect', 'off');
else
    [tFileName,tPathName,tFilterIndex] = uigetfile({'*.dat','Raw Foldometer data files (*.dat)'},'Pick data files',pwd,'MultiSelect', 'off');
end

if isequal(tFileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
    if exist('data')
        
    else
        
        data = 0;
    end
    return
else
    if ischar(tFileName)
        tData = cell(1,1);
        h = waitbar(0,'Importing Data');
        waitbar(0,h);
        disp([getenv('USERNAME'), ' selected ', fullfile(tPathName, tFileName)])
        if binary == 0;
            tData{1} = importdata(fullfile(tPathName, tFileName));
            tData{1}.pathName = fullfile(tPathName, tFileName);
        else
            clear f
            f = getBinary(fullfile(tPathName, tFileName));
            tData{1}.pathName = fullfile(tPathName, tFileName);
            if f.SampleFrequency > 1000;
                tData{1}.textdata{2,1} = {['Sample frequency: ' num2str(f.SampleFrequency/1000) ' kHz']};
            else
                tData{1}.textdata{2,1} = {['Sample frequency: ' num2str(f.SampleFrequency) 'Hz']};
            end
            %f.DmaData.Length/13
            tData{1}.textdata{3,1} = {'Timestamp'};
            tData{1}.data = f.Data;
            tData{1}.data(:,1) = tData{1}.data(:,1) - tData{1}.data(1,1);
%             for ii = 2:12
%                 tData{1}.data(:,ii) = tData{1}.data(:,ii+1);
%             end
            tData{1}.textdata{3,2} = {'FrameNumber'};
            tData{1}.textdata{3,3} = {'Psd1_VxDiff'};
            tData{1}.textdata{3,4} = {'Psd1_VxSum'};
            tData{1}.textdata{3,5} = {'Psd1_VyDiff'};
            tData{1}.textdata{3,6} = {'Psd1_VySum'};
            tData{1}.textdata{3,7} = {'Psd2_VxDiff'};
            tData{1}.textdata{3,8} = {'Psd2_VxSum'};
            tData{1}.textdata{3,9} = {'Psd2_VyDiff'};
            tData{1}.textdata{3,10} = {'Psd2_VySum'};
            tData{1}.textdata{3,11} = {'MirrorX'};
            tData{1}.textdata{3,12} = {'MirrorY'};            
            tData{1}.textdata{3,13} = {'Status'};
            tData{1}.textdata{3,14} = {'Bead1X'};
            tData{1}.textdata{3,15} = {'Bead1Y'};            
            tData{1}.textdata{3,16} = {'Bead2X'};            
            tData{1}.textdata{3,17} = {'Bead2Y'};
        end
        waitbar(1,h);
        close(h)
    else
        Data = cell(1,size(tFileName,2));
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
clear data

data = tData;

%% FORMAT ALL DATA
% clear data
% h = waitbar(0,'Formatting Data');
% waitbar(0,h);
% [file,time,comment] = textread(cData{1}.pathName,'%s %d %s','delimiter','\t');	%reads comment file to the variables file, time, comment using format mask
% filenumber = strtok(tFileName,'.');
% dayfiles = strmatch(date,file);
% trapcommentsindex = strmatch('**k',comment);	%gets comments' indices with values of trap constants
%
% % Get the data from the datafiles and put them in arrays
%
% if not(iscell(tFileName))
%     tFileName = {tFileName};  %fileName becomes cell array
% end
% if not(iscell(filenumber))
%     filenumber = {filenumber};  %filenumber becomes cell array
% end
%
% for ii = 1:size(tFileName,2)
%     waitbar((ii-1)/size(tFileName,2),h);
%     if isempty(strmatch(filenumber(ii),file))
%         commentindex{ii} = 1;
%     else
%         commentindex{ii} = strmatch(filenumber(ii),file);   %index of comments written during the saving of i-file !!gives a mistake if there are no commenta written during saving this file
%     end
%
%     data{ii}.commentindex = commentindex{ii};   %creation of cell array with comments related to i-file
%     data{ii}.commenttime = time(commentindex{ii});
%     data{ii}.commenttext = comment(commentindex{ii});
%     data{ii}.filename = tFileName{ii};
%     data{ii}.pathname = tPathName;
%     data{ii}.date = date;
%     data{ii}.time = tData{ii}.data(:,1);
%     for jj = 2:size(tData{ii}.textdata,2)
%
%         data{ii}.(genvarname(tData{ii}.textdata{jj})) = tData{ii}.data(:,jj);
%
%     end
%
%     % get trap constants from comment file and put them in struct.
%     kindex = trapcommentsindex(find(trapcommentsindex < commentindex{ii}(1),1,'last'));     % looks for last trap constant before this file
%     kindex = [kindex;trapcommentsindex(find(trapcommentsindex >= commentindex{ii}(1),3,'first'))];      % looks for three first trap constants during/after this file
%     % output things to get right trapping parameters
%     constantsstring = comment(kindex);
%     H1 = figure(200);
%     drawnow;
%     pause(0.1);
%     jFrame = get(handle(gcf),'JavaFrame');
%     jFrame.setMaximized(true);
%     pause(0.1);
%
%     colors      =   linspecer(8);
%     set(gcf,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name','Choose Calibration Values');
%     set(gcf,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
%         'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
%     constants = 0;
%     for jj = 1:1:size(constantsstring,1)
%
%         uicontrol('Style','Pushbutton', ...
%             'Position',[0.1 0.7-(jj-1)*(0.2) 0.8 0.1], ...
%             'String',constantsstring{jj},...
%             'FontSize',0.5,...
%             'BackgroundColor',colors(2,:),...
%             'ForegroundColor',colors(5,:),...
%             'Callback',{@delete_graphic,jj,constantsstring});
%     end
%     waitfor(H1);
%     figure(1);
%     constants = getappdata(0, 'constants');
%     data{ii}.kx = constants(1);
%     data{ii}.ky = constants(2);
%     data{ii}.dvdx = constants(3);
%     data{ii}.dvdy = constants(4);
%
%     clear kindex constantsstring nconstants constants
%
%     % calculate forces from trapping parameters and voltages
%     data{ii}.fx = data{ii}.vx*data{ii}.kx/data{ii}.dvdx;
%     data{ii}.fx0 = 0;
%     data{ii}.fy = data{ii}.vy*data{ii}.ky/data{ii}.dvdy;
%     data{ii}.fy0 = 0;
%     data{ii}.ext = -(data{ii}.xpz+data{ii}.vx/data{ii}.dvdx);
%     data{ii}.ext0 = 0;
%     data{ii}.shiftX = 0;
%     data{ii}.shiftY = 0;
%
% end
% waitbar(1,h);
% close(h)
%
% clear cData cFileName cFilterIndex cPathName comment commentindex dayfiles file filenumber h ii jj tData tFileName tFilterIndex tPathName time trapcommentsindex
end

function delete_graphic(hObject,eventdata,A,constantsstring)
constants = sscanf(constantsstring{A}(4:end), '%f %f %f %f');
setappdata(0, 'constants', constants);
close('Choose Calibration Values');
end