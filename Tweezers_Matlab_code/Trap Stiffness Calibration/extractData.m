% LOAD_TIME_SERIES reads the input data.

PATH_current = pwd;

clear Px Py P f chi2X File
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Ask for the input file
[FileName,PathName,FilterIndex] = uigetfile({'*.dat','Tweezers calibration data files (*.dat)'},'Pick file','C:\Users\Roeland\Dropbox\AMOLF PhD\Foldometer\Data\');

%[FileName,PathName,FilterIndex] = uigetfile({'*.dat','Tweezers data files (*.dat)'},'Pick file','\\alanine\Temp\');

if isequal(FileName,0)
    disp([getenv('USERNAME'), ' selected Cancel'])
else
    %     if ischar(FileName)
    calData = cell(1,1);
    disp([getenv('USERNAME'), ' selected ', fullfile(PathName, FileName)])
    if binary == 1;
        clear g
        g = Foldometer(fullfile(PathName, FileName));
        calData{1}.pathName = fullfile(PathName, FileName);
        if g.SampleFrequency > 1000;
            calData{1}.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency/1000) ' kHz']};
        else
            calData{1}.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency) 'Hz']};
        end
        %f.DmaData.Length/13
        calData{1}.textdata{3,1} = {'Timestamp'};
        calData{1}.data = double(g.DmaData);
        calData{1}.data(:,1) = calData{1}.data(:,1) - calData{1}.data(1,1);
        for ii = 2:12
            calData{1}.data(:,ii) = calData{1}.data(:,ii+1);
        end
        calData{1}.data(:,13) = [];
        calData{1}.textdata{3,2} = {'Psd1_VxDiff'};
        calData{1}.textdata{3,3} = {'Psd1_VxSum'};
        calData{1}.textdata{3,4} = {'Psd1_VyDiff'};
        calData{1}.textdata{3,5} = {'Psd1_VySum'};
        calData{1}.textdata{3,6} = {'Psd2_VxDiff'};
        calData{1}.textdata{3,7} = {'Psd2_VxSum'};
        calData{1}.textdata{3,8} = {'Psd2_VyDiff'};
        calData{1}.textdata{3,9} = {'Psd2_VySum'};
        calData{1}.textdata{3,10} = {'MirrorX'};
        calData{1}.textdata{3,11} = {'MirrorY'};
        calData{1}.textdata{3,12} = {'Status'};
        calData{1}.textdata{3,13} = {'FrameNumber'};
        calData{1}.textdata{3,14} = {'Bead1Y'};
        calData{1}.textdata{3,16} = {'Bead2x'};
        calData{1}.textdata{3,17} = {'Bead2Y'};
    else
        calData{1} = importdata(fullfile(PathName, FileName));
        calData{1}.pathName = fullfile(PathName, FileName);
    end
    %     else
    %         calData = cell(1,size(FileName,2));
    %         for ii = 1:size(FileName,2)
    %             disp([getenv('USERNAME'), ' selected ', fullfile(PathName, FileName{ii})])
    %             calData{ii} = importdata(fullfile(PathName, FileName{ii}));
    %             calData{ii}.pathName = fullfile(PathName, FileName{ii});
    %         end
    %         clear ii
    %     end
    
    
    %data{1}.data(:,6) = data{1}.data(:,6)./data{1}.data(:,7); %Normalize Data
    channel_x = 6;
    %data{1}.data(:,8) = data{1}.data(:,8)./data{1}.data(:,9); %Normalize Data
    channel_y = 8;
    
    File = calData{1}.data;
    
    idx = max(strfind(calData{1,1}.pathName,'\'));
    my_path = calData{1,1}.pathName(1:idx);
    my_file_name = calData{1,1}.pathName(idx+1:end);
    FileInfo = dir(calData{1,1}.pathName);
    [time.Y, time.M, time.D, time.H, time.MN, time.S] = datevec(FileInfo.datenum);
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % Ask for the other input parameters
    
    cd(PATH_current)
    
    input_para1;
    input_para2;
    input_para3;
    
end


