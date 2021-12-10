% VIEW_DATA shows the position histogram by calling PLOT_HISTOGRAM 
% and the powerspectrum by calling CALC_POWERSP and PLOT_POWERSPECTRUM.
% DECORR_XY is called to decorrelate the channels if WANT_DECORR > 0.
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Performs various checks on input parameters
check1(nblock,color1); % check if nblock is larger than zero
check(sampling_f,0,1000,'sampling frequency',color1);
check2(channel_x,1,size(File,2),'Channel X',color1);
check2(channel_y,0,size(File,2),'Channel Y',color1);
check2(channel_z,0,size(File,2),'Channel Z',color1);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
disp(['channel x = ',num2str(channel_x)]);
disp(['channel y = ',num2str(channel_y)]);
disp(['channel z = ',num2str(channel_z)]);

screensize = get(0,'ScreenSize');
P2 = get(gcf,'pos');
lef = 0; bot = P2(2)+5; wid = screensize(3); hei = screensize(4) - P2(2)*2-20;

if channel_z ~= 0,
    X       = File(:,channel_x)./File(:,channel_z); 
else
    X       = File(:,channel_x);
end;  %(if channel_z ~= 0)
clear fS PxS PyS TS
% Make tMsr = N*tDrive Sinusoidal calibration at 32 Hz
if nCal == 0;
    %X = X(100001:200000);
    tRemove = mod(length(X)*1/(1000*sampling_f),1/32);
    pRemove = tRemove * 1000*sampling_f;
    X = X(1:end-pRemove);
    X = X - mean(X);
    for ii = 1:length(X)/(1000*sampling_f)
        [fS{ii},PxS{ii},TS{ii}]  = calc_powersp(X(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)),1000*sampling_f);
    end
    f = mean([fS{:}],2);
    Px = mean([PxS{:}],2);
    T = mean([TS{:}],2);  
else    
   X = X - mean(X);
   [f,Px,T]  = calc_powersp(X,1000*sampling_f);
end

%-------------------------

figure(6); clf; set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name','PSD with crosstalk');  plot_powerspectrum(f,Px,color_x,1);
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
hold on; ...
    title('Power spectra. Elimination of cross-talk not applied');
figure(8); clf; hold on; set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name','X PSD, few datapoints per block');  plot_powerspec_lin(f,Px,fNyq,color_x); 
    title('X PSD, few datapoints per block');
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);

if channel_y ~= 0,  
    if channel_z ~= 0,
        Y       = File(:,channel_y)./File(:,channel_z);
    else
        Y       = File(:,channel_y);
    end; %(if channel_z ~= 0)
%     if nCal == 0;
%         Y = Y(100001:200000);
%         Y = Y(1:end-pRemove);
%         %Y = Y(1:300000);
%     end
%     Y       = Y - mean(Y);
%     [f,Py]  = calc_powersp(Y,1000*sampling_f);

    if nCal == 0;
        Y = Y(1:end-pRemove);
        Y = Y - mean(Y);
        for ii = 1:length(Y)/(1000*sampling_f)
            [fS{ii},PyS{ii},TS{ii}]  = calc_powersp(Y(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)),1000*sampling_f);
        end
        f = mean([fS{:}],2);
        Py = mean([PyS{:}],2);
        T = mean([TS{:}],2);
    else
        Y = Y - mean(Y);
        [f,Py,T]  = calc_powersp(Y,1000*sampling_f);
    end
    
    figure(6); plot_powerspectrum(f,Py,color_y,2);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    plot_hx = findobj('color',color_x);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    plot_hy = findobj('color',color_y);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    legend([plot_hx(1) plot_hy(1)],{'X','Y'});
    figure(9); clf; hold on; set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name','Y PSD, few datapoints per block');  plot_powerspec_lin(f,Py,fNyq,color_y); 
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    title('Y PSD, few datapoints per block');
    IXY = [1 2];
    hold off
if want_decorr > 0, decorr_xy; end; %(if want_decorr > 0)
    figure(5); clf; set(gcf,'Numbertitle','off','Name','Distribution of Y-signal'); pos = Y; Tit = 'Y-'; plot_histogram;
    %set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name','Distribution of Y-signal'); pos = Y; Tit = 'Y-'; plot_histogram;
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
else IXY = 1;
end; %(if channel_y ~= 0)
figure(4); clf; set(gcf,'Numbertitle','off','Name','Distribution of X-signal'); pos = X; Tit = 'X-'; plot_histogram; 
%set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name','Distribution of X-signal'); pos = X; Tit = 'X-'; plot_histogram; 
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
