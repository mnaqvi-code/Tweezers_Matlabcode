%close all; clear all; clc;
clear n tolx fNyq T nblock channel_x channel_y channel_z want_alias want_hydro R l nu rho rhol want_filter elec_filters want_hist parametersX parametersY scal_fitX scal_fitY parameters0 Plot_start Plot_end Ffit_start Ffit_end;
cPos = 1;
nCal = 1;
global prgVr calData n tolx fNyq T nblock channel_x channel_y channel_z want_alias want_hydro R l nu rho rhol want_filter elec_filters want_hist parametersX parametersY scal_fitX scal_fitY parameters0 Plot_start Plot_end Ffit_start Ffit_end

%screensize = get(0,'ScreenSize');
figure(1); clf; hold on;
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
pause(0.1);
P2 = get(gcf,'pos');
lef = P2(1); bot = P2(2); wid = P2(3); hei = P2(4);
%set(gcf, 'Position', get(0,'Screensize'));
set(0,'DefaultFigurePosition', [lef bot wid hei],'DefaultAxesPosition', [0.1  0.15  0.8 0.75]);
set(gcf,'BackingStore','off','MenuBar','figure',...
    'NumberTitle','off','Name','Power spectrum analysis');
set(gcf,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
    'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
set(gca,'handlevisibility','off','visible','off');
% This creates the 'background' axes
ha = axes('units','normalized', ...
            'position',[0 0 1 1]);
% Move the background axes to the bottom
uistack(ha,'bottom');
% Load in a background image and display it using the correct colors
% The image used below, is in the Image Processing Toolbox.  If you do not have %access to this toolbox, you can use another image file instead.
I=imread('BackgroundImage.jpg');
hi = imagesc(I);
colormap gray;
% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
set(ha,'handlevisibility','off', ...
            'visible','off');


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Positions of text boxes in the main window

edge        =   0.03;
wid1        =   (1.2)*((1-9*edge)/4); 
wid2        =   (0.8)*((1-9*edge)/4); 
left1       =   2*edge;
left2       =   3*edge+wid1;
left3       =   6*edge+wid1+wid2;
left4       =   7*edge+2*wid1+wid2;
titlefs     =   0.5;
bot_title   =   0.8;
for i = 1 : 14, eval(['bot' num2str(i) '= bot_title - 0.01 - i * 0.06;']); end;
for i = 1 : 16, eval(['bbot' num2str(i) '= bot_title - 0.01 - i * 0.05;']); end;
height      =   0.045;
color       =   [0.7 0.7 0.7];
color1      =   [0.6 0.6 0.9];
color2      =   'w';
color_x     =   [0.6 0.1 0.8];
color_y     =   [0.2 0.7 0.3];

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Title and the LOAD_TIME_SERIES button

hpc  = uicontrol('Style','text','Position',[0.35 0.92 0.3 height+0.01],'String',' Power spectrum analysis',...
    'FontSize',titlefs+0.2,'ForegroundColor','w','BackgroundColor','k');
jh = findjobj(hpc);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
ht  = uicontrol('Style','text', ...
    'Position',[edge 0.02 wid1+wid2+3*edge bot_title+height], ...
    'String',' ',...
    'BackgroundColor',color);
% hload_time_series  = uicontrol('Style','Pushbutton', ...
%     'Position',[left1+0.5*wid1 bot_title 0.5*(wid1+wid2) height], ...
%     'String','Load time series',...
%     'BackgroundColor',[1 0.1 0.1],...
%     'Callback','load_time_series');

hload_time_series  = uicontrol('Style','Pushbutton', ...
    'Position',[left1+0.5*wid1 bot_title 0.5*(wid1+wid2) height], ...
    'String','Load time series',...
    'FontSize',titlefs+0.2,...
    'BackgroundColor',[1 0.1 0.1],...
    'Callback','extractDataQ');
hreturn_menu  = uicontrol('Style','Pushbutton', ...
    'Position',[left1+0.1*wid1 0.92 0.5*(wid1+wid2) height], ...
    'String','Return to Main Menu',...
    'FontSize',titlefs+0.2,...
    'BackgroundColor',[1 0.1 0.1],...
    'Callback','tweezerStart');
if prgVr.setupStt == 1
    cToggle  = uicontrol('Style','Pushbutton', ...
         'Position',[left1+0.1*wid1+0.6 0.92 0.5*(wid1+wid2) height], ...
         'String','Normal Calibration',...
         'FontSize',titlefs+0.2,...
         'BackgroundColor',[1 0.1 0.1],...
         'Callback','if(nCal == 1) nCal = 0; set(cToggle, ''String'', ''Sinusoidal Calibration''); else nCal = 1; set(cToggle, ''String'', ''Normal Calibration''); end;');
end
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
