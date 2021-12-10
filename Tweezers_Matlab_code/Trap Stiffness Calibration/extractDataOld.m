% LOAD_TIME_SERIES reads the input data.

PATH_current = pwd;
global binaryD calData

clear Px Py P f chi2X File
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Ask for the input file
%[FileName,PathName,FilterIndex] = uigetfile({'*.dat','Tweezers calibration data files (*.dat)'},'Pick file','C:\Users\Roeland\Dropbox\AMOLF PhD\Foldometer\Data\');

if ~isempty(calData)
    [FileName,PathName,FilterIndex] = uigetfile({'*.dat;*.cal;* Power Spectrum.tdms','Calibration File (*.dat,*.cal,*.tdms)'},'Pick Calibration File',strrep(calData.pathName, '\', filesep),'MultiSelect', 'off');
else
    if exist('accessHistory.txt', 'file') == 2
        fid = fopen( 'accessHistory.txt', 'r' );
        hisRead = textscan(fid, '%s %s %s %s %s %s %s\n', 'delimiter', '|','collectoutput',true);
        fclose(fid);
        if exist(hisRead{1}{end,5},'dir') == 7
            [FileName,PathName,FilterIndex] = uigetfile({'*.dat;*.cal;* Power Spectrum.tdms','Calibration File (*.dat,*.cal,*.tdms)'},'Pick Calibration File',hisRead{1}{end,5},'MultiSelect', 'off');
        else
            [FileName,PathName,tFilterIndex] = uigetfile({'*.dat;*.cal;* Power Spectrum.tdms','Calibration File (*.dat,*.cal,*.tdms)'},'Pick Calibration File',pwd,'MultiSelect', 'off');
        end
    else
        [FileName,PathName,FilterIndex] = uigetfile({'*.dat;*.cal;* Power Spectrum.tdms','Calibration File (*.dat,*.cal,*.tdms)'},'Pick Calibration File',pwd,'MultiSelect', 'off');
    end
    
end

if isequal(FileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
else
    if exist('accessHistory.txt', 'file') == 2
        fid = fopen( 'accessHistory.txt', 'a+' );
        fprintf( fid, '%s| %s| %s| %s| %s| %s| %s\n', getenv('USERNAME'), 'selected', FileName, 'in', PathName, 'at', datestr(datetime('now')));
        fclose(fid);
    else
        fid = fopen( 'accessHistory.txt', 'a+' );
        fprintf( fid, '%s| %s| %s| %s| %s| %s| %s', getenv('USERNAME'), 'selected', FileName, 'in', PathName, 'at', datestr(datetime('now')));
        fclose(fid);
    end
    
    %     if ischar(FileName)
    calData = [];
    disp([getenv('USERNAME'), ' selected ', fullfile(PathName, FileName)])
    if prgVr.setupStt == 0
        calData.data = importdata(fullfile(PathName, FileName));
        calData.pathName = fullfile(PathName, FileName);
        %channel_x = 1;    
        %channel_y = 2;
        channel(1) = 1;
        channel(2) = 2;
        File = load(calData.pathName);
    elseif prgVr.setupStt == 1
        clear g
        g = getBinary(fullfile(PathName, FileName));
        calData.pathName = fullfile(PathName, FileName);
        if g.SampleFrequency > 1000
            calData.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency/1000) ' kHz']};
        else
            calData.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency) 'Hz']};
        end
        %f.DmaData.Length/13
        calData.textdata{3,1} = {'Timestamp'};
        calData.data = g.Data;
        calData.data(:,1) = calData.data(:,1) - calData.data(1,1);
        for ii = 2:12
            calData.data(:,ii) = calData.data(:,ii+1);
        end
        calData.data(:,13) = [];
        calData.textdata{3,2} = {'Psd1_VxDiff'};
        calData.textdata{3,3} = {'Psd1_VxSum'};
        calData.textdata{3,4} = {'Psd1_VyDiff'};
        calData.textdata{3,5} = {'Psd1_VySum'};
        calData.textdata{3,6} = {'Psd2_VxDiff'};
        calData.textdata{3,7} = {'Psd2_VxSum'};
        calData.textdata{3,8} = {'Psd2_VyDiff'};
        calData.textdata{3,9} = {'Psd2_VySum'};
        calData.textdata{3,10} = {'MirrorX'};
        calData.textdata{3,11} = {'MirrorY'};
        calData.textdata{3,12} = {'Status'};
        calData.textdata{3,13} = {'FrameNumber'};
        calData.textdata{3,16} = {'Bead1X'};
        calData.textdata{3,15} = {'Bead1Y'};
        calData.textdata{3,16} = {'Bead2X'};
        calData.textdata{3,17} = {'Bead2Y'};
        %channel_x = 3;    
        %channel_y = 5;
        channel(1) = 2;
        channel(2) = 4;
        channel(3) = 6;
        channel(4) = 8;
        File = calData.data;
    elseif prgVr.setupStt == 2
        calData = getCTrapDataDirect(false,fullfile(PathName, FileName));
        calData.data = calData.Data.MeasuredData;
        calData.pathName = fullfile(PathName, FileName);
        if calData.data(2).Property(6).Value > 1000
            calData.textdata{2,1} = {['Sample frequency: ' num2str(calData.data(2).Property(6).Value/1000) ' kHz']};
        else
            calData.textdata{2,1} = {['Sample frequency: ' num2str(calData.data(2).Property(6).Value) 'Hz']};
        end
        File = [calData.data(3:7).Data];
        channel(1) = 2;
        channel(2) = 3;
        channel(3) = 4;
        channel(4) = 5;
    end
    
    %     else
    %         calData{1} = cell(1,size(FileName,2));
    %         for ii = 1:size(FileName,2)
    %             disp([getenv('USERNAME'), ' selected ', fullfile(PathName, FileName{ii})])
    %             calData{1}{ii} = importdata(fullfile(PathName, FileName{ii}));
    %             calData{1}{ii}.pathName = fullfile(PathName, FileName{ii});
    %         end
    %         clear ii
    %     end
    
    
    %data.data(:,6) = data.data(:,6)./data.data(:,7); %Normalize Data
    %channel_x = 1;
    %data.data(:,8) = data.data(:,8)./data.data(:,9); %Normalize Data
    %channel_y = 2;
    %File = load(calData.pathName);
    
    idx = max(strfind(calData.pathName,'\'));
    my_path = calData.pathName(1:idx);
    my_file_name = calData.pathName(idx+1:end);
    FileInfo = dir(calData.pathName);
    [time.Y, time.M, time.D, time.H, time.MN, time.S] = datevec(FileInfo.datenum);    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % Ask for the other input parameters
    
    cd(PATH_current)
    
    input_para1;
    input_para2;
    input_para3;
    
end




