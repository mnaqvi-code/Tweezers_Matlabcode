function [] = plotAll( )

global data wndVr pltVr; %#ok<REDEF>

fh = findall(0,'type','figure');
for ii = 1:length(fh)
    if fh(ii) ~= 1
        close(fh(ii))
    end
end
gg = figure(2);
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
pause(0.1);
P2 = get(gcf,'pos');
set(0,'DefaultFigurePosition', P2,'DefaultAxesPosition', [0.1,0.1,0.8,0.8]);
if pltVr.splitPulls == 1
    C = distinguishable_colors(size(data.scombined,2)*2);
    %C = linspecer(size(data.scombined,2)*2,'qualitative');
else
    C = distinguishable_colors(size(data.scombined,2)*2);
    %C = linspecer(size(data.scombined,2),'qualitative');
end
if pltVr.grphAct == 1
    for jj =  1:length(data.x)/2 
        hold on
        plotAllft(jj) = plot(data.figs.ft.axh(jj).XData,data.figs.ft.axh(jj).YData,'color',data.figs.ft.axh(jj).Color);
    end
    hold off
    title(strrep(data.filename, '_', ' '));
    xlabel('Time (s)')
    ylabel('Force (pN)');
    set(gca,'Tag',num2str(1))
    
    if pltVr.legOn == 1        
        legend(plotAllft(1:size(data.scombined,2)),pltVr.legArry)
        legend('show')
    else
        legend('hide')
    end
end

if pltVr.grphAct == 2
    for jj =  1:length(data.x)/2 
        hold on
        plotAllte(jj) = plot(data.figs.te.axh(jj).XData,data.figs.te.axh(jj).YData,'color',data.figs.te.axh(jj).Color);
    end
    hold off
    title(strrep(data.filename, '_', ' '));
    xlabel('Time (s)')
    ylabel('Extension (\mum)');
    set(gca,'Tag',num2str(1))
    
    
    data.figs.teUnf.axh(ii)
    
    if pltVr.legOn == 1        
        legend(plotAllte(1:size(data.scombined,2)),pltVr.legArry)
        legend('show')
    else
        legend('hide')
    end
    
end

if pltVr.grphAct == 3
    for jj = 1:size(data.scombined,2)
        if pltVr.splitPulls == 1
            for gg = 1:2
                plotAllf(jj*2-1+gg-1) = plot(data.scombinedSplit{jj*2-1+gg-1}+data.shiftX, data.forceSplit{jj*2-1+gg-1}+data.shiftY,'color',C(jj*2-1+gg-1,:));
                hold on
            end
        else
            plotAllf(jj) = plot(data.scombinedSplit{jj}+data.shiftX, data.force{jj}+data.shiftY,'color',C(jj,:));
            hold on
        end
    end
    set(gca,'Fontsize',16)
    title(strrep(data.filename, '_', ' '));
    set(gg,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name',data.filename);
    set(gg,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
        'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
    xlabel('Extension (µm)')
    ylabel('Force (pN)');
    
    
    if pltVr.legOn == 1
        if pltVr.splitPulls == 1            
            legend(plotAllf(1:(size(data.scombined,2)*2)),pltVr.legArry)            
        else
            legend(plotAllf(jj),pltVr.legArry(jj))
        end        
        legend('show')
    else
        legend('hide')
    end
    
    plotRuler();
    
    hold off;
    
end

if pltVr.grphAct == 4
    for jj =  1:length(data.x)/2 
        hold on
        plotAllfs(jj) = plot(data.figs.fs.axh(jj).XData,data.figs.fs.axh(jj).YData,'color',data.figs.fs.axh(jj).Color);
        plotAllfsfit(jj) = plot(data.figs.fsfit.axh(jj).XData,data.figs.fsfit.axh(jj).YData,'color','k');
    end
    hold off
    title(data.filename);
    xlabel('Xspt2 (\mum)')
    ylabel('Voltage (V)');
    set(gca,'Tag',num2str(1))
    
    if pltVr.legOn == 1        
        legend(plotAllfs(1:size(data.scombined,2)),pltVr.legArry)
        legend('show')
    else
        legend('hide')
    end
end

set(gca,'box','off')
xlim([str2num(get(wndVr.sampling_xlim,'String'))]);
ylim([str2num(get(wndVr.sampling_ylim,'String'))]);
hold off

end