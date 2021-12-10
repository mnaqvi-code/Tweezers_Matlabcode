% VIEW_DATA shows the position histogram by calling PLOT_HISTOGRAM
% and the powerspectrum by calling CALC_POWERSP and PLOT_POWERSPECTRUM.
% DECORR_XY is called to decorrelate the channels if WANT_DECORR > 0.
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

global pltVr sampling_f

% Performs various checks on input parameters
check1(nblock,color1); % check if nblock is larger than zero
check(sampling_f,0,1000,'sampling frequency',color1);
check2(channel(1),1,size(File,2),'Channel X-1',color1);
check2(channel(2),0,size(File,2),'Channel Y-1',color1);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
disp(['channel X-1 = ',num2str(channel(1))]);
disp(['channel Y-1 = ',num2str(channel(2))]);
if prgVr.setupStt == 1 || prgVr.setupStt == 2
    if channel(3) ~= 0 && channel(4) ~= 0
        check2(channel(3),0,size(File,2),'Channel X-2',color1);
        check2(channel(4),0,size(File,2),'Channel Y-2',color1);
        disp(['channel X-2 = ',num2str(channel(3))]);
        disp(['channel Y-2 = ',num2str(channel(4))]);
        cIter = 2;
    else
        cIter = 1;
    end
    
else
    cIter = 1;
end


screensize = get(0,'ScreenSize');
P2 = get(gcf,'pos');
lef = 0; bot = P2(2)+5; wid = screensize(3); hei = screensize(4) - P2(2)*2-20;

% if channel_z ~= 0
%     X       = File(:,channel_x)./File(:,channel_z);
% else
%     X       = File(:,channel_x);
% end;
%(if channel_z ~= 0)
clear fS PxS PyS TS X Y

for jj = 1:cIter
    clear fS PxS PyS TS
    
    if isfield(pltVr.figs, 'PSDNorm') && jj == 1
        if ~isempty(findobj('name','PSD With No Crosstalk'))
            close(pltVr.figs.PSDNorm)
        end
        rmfield(pltVr.figs, 'PSDNorm');
    end
    
    if isfield(pltVr.figs, 'elCross') && jj == 1
        if ~isempty(findobj('name','Elimination of crosstalk'))
            close(pltVr.figs.elCross)
        end
        rmfield(pltVr.figs, 'elCross');
    end
    
    if isfield(pltVr.figs, 'PSDbin') && jj == 1
        if ~isempty(findobj('name','PSDs, few datapoints per block'))
            %close(findobj('name','PSDs, few datapoints per block'))
            close(pltVr.figs.PSDbin)
        end
        rmfield(pltVr.figs, 'PSDbin');
    end
    
    X{jj} = File(:,channel(jj*2-1));
    Y{jj} = File(:,channel(jj*2));
    
    % Make tMsr = N*tDrive Sinusoidal calibration at 32 Hz
    if nCal == 0
        %X{jj} = X{jj}(100001:200000);
        tRemove = mod(length(X{jj})*1/(1000*sampling_f),1/32);
        pRemove = tRemove * 1000*sampling_f;
        X{jj} = X{jj}(1:end-pRemove);
        X{jj} = X{jj} - mean(X{jj});
        for ii = 1:length(X{jj})/(1000*sampling_f)
            [fS{ii},PxS{ii},TS{ii}]  = calc_powersp(X{jj}(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)),1000*sampling_f);
        end
        f = mean([fS{:}],2);
        Px{jj} = mean([PxS{:}],2);
        T = mean([TS{:}],2);
    else
        X{jj} = X{jj} - mean(X{jj});
        [f,Px{jj},T]  = calc_powersp(X{jj},1000*sampling_f);
    end
    
    %-------------------------
    
    if isfield(pltVr.figs, 'PSDwCross') && jj == 1
        if ~isempty(findobj('name','PSD with crosstalk'))
            close(pltVr.figs.PSDwCross)
        end
        rmfield(pltVr.figs, 'PSDwCross');
    end
    
    if jj == 1
        pltVr.figs.PSDwCross = figure;
        set(gcf,'Numbertitle','off','Name','PSD with crosstalk');
        hf2 = suptitle('Power spectra. Elimination of cross-talk not applied');
        set(hf2,'FontSize',20,'FontWeight','Bold')
        %set(hf2,'Fontweight','Bold');
    elseif jj > 1
        set(0, 'currentfigure', pltVr.figs.PSDwCross);
    end
    subplot(1,cIter,jj);
    hold on;
    set(0, 'currentfigure', pltVr.figs.PSDwCross); shg; plot_powerspectrum(f,Px{jj},color_x,1);
    hold off;
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    hold on; ...
        
%title('Power spectra. Elimination of cross-talk not applied');
if jj == 1
    pltVr.figs.PSDbin = figure;
    set(gcf,'Numbertitle','off','Name','PSDs, few datapoints per block');
    hf2 = suptitle('PSDs, few datapoints per block');
    set(hf2,'FontSize',20,'FontWeight','Bold')
    %set(hf2,'Fontweight','Bold');
elseif jj > 1
    set(0, 'currentfigure', pltVr.figs.PSDbin);
end

hold on;
set(0, 'currentfigure', pltVr.figs.PSDbin); subplot(2,cIter,jj*2-1); shg; plot_powerspec_lin(f,Px{jj},fNyq,color_x);
hold off;
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);

if channel(jj*2) ~= 0
    %         if channel_z ~= 0
    %             Y{jj}       = File(:,channel_y)./File(:,channel_z);
    %         else
    %             Y{jj}       = File(:,channel_y);
    %         end;
    %(if channel_z ~= 0)
    %     if nCal == 0;
    %         Y{jj} = Y{jj}(100001:200000);
    %         Y{jj} = Y{jj}(1:end-pRemove);
    %         %Y{jj} = Y{jj}(1:300000);
    %     end
    %     Y{jj}       = Y{jj} - mean(Y{jj});
    %     [f,Py]  = calc_powersp(Y{jj},1000*sampling_f);
    
    if nCal == 0
        Y{jj} = Y{jj}(1:end-pRemove);
        Y{jj} = Y{jj} - mean(Y{jj});
        for ii = 1:length(Y{jj})/(1000*sampling_f)
            [fS{ii},PyS{ii},TS{ii}]  = calc_powersp(Y{jj}(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)),1000*sampling_f);
        end
        f = mean([fS{:}],2);
        Py{jj} = mean([PyS{:}],2);
        T = mean([TS{:}],2);
    else
        Y{jj} = Y{jj} - mean(Y{jj});
        [f,Py{jj},T]  = calc_powersp(Y{jj},1000*sampling_f);
    end
    hold on;
    set(0, 'currentfigure', pltVr.figs.PSDwCross); shg; subplot(1,cIter,jj); plot_powerspectrum(f,Py{jj},color_y,2);
    hold off;
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
    
    hold on;
    set(0, 'currentfigure', pltVr.figs.PSDbin); subplot(2,cIter,jj*2); plot_powerspec_lin(f,Py{jj},fNyq,color_y);
    hold off;
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    IXY = [1 2];
    hold off
    
    if want_decorr > 0
        decorr_xy;
    end %(if want_decorr > 0)
    if isfield(pltVr.figs, 'GaussbinX') && jj == 1
        if iscell(pltVr.figs.GaussbinX)
            for kk = 1:length(pltVr.figs.GaussbinX)
                if ~isempty(findobj('name',['Distribution of X-', num2str(kk), ' signal']))
                    close(pltVr.figs.GaussbinX{kk})
                end
            end
        else
            if ~isempty(findobj('name',['Distribution of X-', num2str(jj), ' signal']))
                close(pltVr.figs.GaussbinX{jj})
            end
        end
        rmfield(pltVr.figs, 'GaussbinX');
    end
    
    if isfield(pltVr.figs, 'GaussbinY')  && jj == 1
        if iscell(pltVr.figs.GaussbinY)
            for kk = 1:length(pltVr.figs.GaussbinY)
                if ~isempty(findobj('name',['Distribution of Y-', num2str(kk), ' signal']))
                    close(pltVr.figs.GaussbinY{kk})
                end
            end
        else
            if ~isempty(findobj('name',['Distribution of Y-', num2str(jj), ' signal']))
                close(pltVr.figs.GaussbinY{jj})
            end
        end
        rmfield(pltVr.figs, 'GaussbinY');
    end
    
    pltVr.figs.GaussbinY{jj} = figure;
    set(gcf,'Numbertitle','off','Name',['Distribution of Y-', num2str(jj), ' signal']);
    hold on;
    set(0, 'currentfigure', pltVr.figs.GaussbinY{jj}); shg; pos = Y{jj}; Tit = 'Y-'; clf; plot_histogram;
    hold off;
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
else
    IXY = 1;
end
pltVr.figs.GaussbinX{jj} = figure;
set(gcf,'Numbertitle','off','Name',['Distribution of X-', num2str(jj), ' signal']);
hold on;
set(0, 'currentfigure', pltVr.figs.GaussbinX{jj}); shg; pos = X{jj}; Tit = 'X-'; clf; plot_histogram;
hold off;

drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
end