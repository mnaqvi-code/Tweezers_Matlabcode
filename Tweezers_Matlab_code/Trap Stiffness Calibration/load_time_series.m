% LOAD_TIME_SERIES reads the input data.

PATH_current = pwd;

clear Px Py P f chi2X File
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Ask for the input file

[filename, pathname] = uigetfile('*.*','Choose file containing time series to be analyzed');  %Open file dialog box
cd(pathname);
File = load(filename);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Find channels X, Y and Z. Default values are:
if size(File,2) > 2,
% divide by values in z-channel, i.e. normalize!
    
    channel_x = 1;
    channel_y = 2;
    channel_z = 3;
    
elseif size(File,2) == 2,
    
    channel_x   =   1;
    channel_y   =   2;
    channel_z   =   0;
    
elseif size(File,2) == 1,
    
    channel_x   =   1;
    channel_y   =   0;
    channel_z   =   0;

end;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Ask for the other input parameters

cd(PATH_current)

input_para1;
input_para2;
input_para3;

