% LOAD_TIME_SERIES reads the input data.

PATH_current = pwd;
clear n tolx fNyq T nblock channel_x channel_y channel_z want_alias want_hydro R l nu rho rhol want_filter elec_filters want_hist parametersX parametersY scal_fitX scal_fitY parameters0 Plot_start Plot_end Ffit_start Ffit_end Px Py P f chi2X File;

global prgVr calData n tolx fNyq T nblock channel_x channel_y channel_z want_alias want_hydro R l nu rho rhol want_filter elec_filters want_hist parametersX parametersY scal_fitX scal_fitY parameters0 Plot_start Plot_end Ffit_start Ffit_end
color1 = [0.6 0.6 0.9];
color_x = [0.6 0.1 0.8];
color_y = [0.2 0.7 0.3];
cPos = 1;
nCal = 1;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Ask for the input file
%[FileName,PathName,FilterIndex] = uigetfile({'*.dat','Tweezers calibration data files (*.dat)'},'Pick file','C:\Users\Roeland\Dropbox\AMOLF PhD\Foldometer\Data\');

if ~isempty(calData)
    [FileName,PathName,FilterIndex] = uigetfile({'*.dat;*.cal','Calibration File (*.dat,*.cal)'},'Pick Calibration File',strrep(calData.pathName, '\', filesep),'MultiSelect', 'off');
else
    if exist('accessHistory.txt', 'file') == 2
        fid = fopen( 'accessHistory.txt', 'r' );
        hisRead = textscan(fid, '%s %s %s %s %s %s %s\n', 'delimiter', '|','collectoutput',true);
        fclose(fid);
        if exist(hisRead{1}{end,5},'dir') == 7
            [FileName,PathName,FilterIndex] = uigetfile({'*.dat;*.cal','Calibration File (*.dat,*.cal)'},'Pick Calibration File',hisRead{1}{end,5},'MultiSelect', 'off');
        else
            [FileName,PathName,tFilterIndex] = uigetfile({'*.dat;*.cal','Calibration File (*.dat,*.cal)'},'Pick Calibration File',pwd,'MultiSelect', 'off');
        end
    else
        [FileName,PathName,FilterIndex] = uigetfile({'*.dat;*.cal','Calibration File (*.dat,*.cal)'},'Pick Calibration File',pwd,'MultiSelect', 'off');
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
    if binaryD == 1;
        clear g
        g = getBinary(fullfile(PathName, FileName));
        calData.pathName = fullfile(PathName, FileName);
        if g.SampleFrequency > 1000;
            calData.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency/1000) ' kHz']};
        else
            calData.textdata{2,1} = {['Sample frequency: ' num2str(g.SampleFrequency) 'Hz']};
        end
        %f.DmaData.Length/13
        calData.textdata{3,1} = {'Timestamp'};
        calData.data = g.Data
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
        calData.textdata{3,14} = {'Bead1Y'};
        calData.textdata{3,16} = {'Bead2x'};
        calData.textdata{3,17} = {'Bead2Y'};
        channel_x = 3;    
        channel_y = 5;
        File = calData.data;
    else
        calData.data = importdata(fullfile(PathName, FileName));
        calData.pathName = fullfile(PathName, FileName);
        channel_x = 1;    
        channel_y = 2;
        File = load(calData.pathName);
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
    
    
    n           =   0;              %number of aliasing terms to include, will be set later to a nonzero value
    sampling_f  =   100;             % (kHz)  ---> change here if you want another default sampling frequency
    fNyq        =   1000*sampling_f/2;   % Nyquist frequency. In Hz.
    want_decorr =   0;
    channel_x = 2;    
    channel_y = 4;
    
    view_dataQ;
    
    nblock      =   150;   % number of data points in a block
    tolx        =   1e-7;   % tolerance in the fit
    Lfit_start  =   100;    % range of the Lorentzian fit (Hz)
    Lfit_end    =   20000;   %  -"-
    Ffit_start  =   100;    % range of the final fit (Hz)
    Ffit_end    =   20000;   %  -"-
    Plot_start  =   100;     % plotting range (Hz)
    Plot_end    =   20000;   %  -"-
    want_alias  =   1;      % want_alias = 1 for accounting for aliasing;
    want_diode  =   1;      % want_diode = 0 for no filtering by diode;
    want_alpha  =   1;      % want_alpha = 0 for filtering by diode with 1 parameter (f3dBeff); want_alpha = 1 for fitting filtering by diode with 2 parameters (alpha and f3dB); 

    
    want_hydro  =   0;              % want_hydro = 1 for including hydrodynamic corrections; otherwise 0 
    R           =   0.5e-6;         % radius of the bead (m)
    l           =   5e-6;           % distance to surface (m)
    nu          =   1e-12/1e-6;     % kinematic viscosity (m2/s)
    rho         =   1e-3/1e-6;      % density of the bead (kg/m3)
    
    fit_powerspectrum;
    
end


