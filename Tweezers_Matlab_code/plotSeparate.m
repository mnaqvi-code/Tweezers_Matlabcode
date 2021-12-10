
function [] = plotSeparate( )

global data wndVr pltVr datVr; %#ok<REDEF>

fh = findall(0,'type','figure');
formatTimePointOut = 'yymmdd HHMMSS';
for ii = 1:length(fh)
    if fh(ii) ~= 1
        close(fh(ii))
    end
end
if pltVr.splitPulls == 1
    C = distinguishable_colors(size(data.scombined,2)*2);
    %C = linspecer(size(data.scombined,2)*2,'qualitative');
else
    C = distinguishable_colors(size(data.scombined,2));
    %C = linspecer(size(data.scombined,2),'qualitative');
end

kk = 1;
for jj = 1:size(data.scombined,2)
    if data.x(kk) < 0
        selectedInterval = '';
        selectedInterval = strcat('0-',num2str(round(data.x(kk+1))))
    else
        selectedInterval = '';
        selectedInterval = strcat(num2str(round(data.x(kk))),'-',num2str(round(data.x(kk+1))))
    end
    kk = kk + 2;
    f = figure(jj+1);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    pause(0.1);
    P2 = get(gcf,'pos');
    set(0,'DefaultFigurePosition', P2,'DefaultAxesPosition', [0.1,0.1,0.8,0.8]);
    
    hold on
        if pltVr.splitPulls == 1
            for gg = 1:2
                plotSepf(jj*2-1+gg-1) = plot(data.scombinedSplit{jj*2-1+gg-1}+data.shiftX, data.forceSplit{jj*2-1+gg-1}+data.shiftY,'color',C(jj*2-1+gg-1,:));
            end
        else    
            plotSepf(jj) = plot(data.scombinedSplit{jj}+data.shiftX, data.force{jj}+data.shiftY,'color',C(jj,:));
        end
    
    set(gca,'Fontsize',16)
    title({strrep(data.filename, '_', ' ');strcat('Time interval: ',32,selectedInterval)})
    set(f,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name',data.filename);
    set(f,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
        'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
    xlabel('Extension (µm)')
    ylabel('Force (pN)');
    
    if pltVr.legOn == 1
        if pltVr.splitPulls == 1            
            legend(plotSepf(jj*2-1:jj*2),pltVr.legArry(jj*2-1:jj*2))            
        else
            legend(plotSepf(jj),pltVr.legArry(jj))
        end        
        legend('show')
    else
        legend('hide')
    end
    
    plotRuler();
    
    xlim([str2num(get(wndVr.sampling_xlim,'String'))]);
    ylim([str2num(get(wndVr.sampling_ylim,'String'))]);
    hold off
    [data(:).ExpInfo,data(:).ExpNumber] = strtok(data.filename,'#');
    data(:).ExpNumber = strtok(data.ExpNumber,'.');
    [data(:).SaveFolder] = strcat(data.FileFolder,'\Analysis\',data.ExpNumber,'\');
    mkdir(data.SaveFolder);
    saveas(gcf, fullfile(data.SaveFolder,strcat(datestr(now,formatTimePointOut),32,selectedInterval)),'fig');
    saveas(gcf, fullfile(data.SaveFolder,strcat(datestr(now,formatTimePointOut),32,selectedInterval)),'epsc');
end
saveConstantData = datVr.dualCst;
save(fullfile(data.SaveFolder,strcat(datestr(now,formatTimePointOut),32,selectedInterval,' ConstantData.mat')),'saveConstantData');
end