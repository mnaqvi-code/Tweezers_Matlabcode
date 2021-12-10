function [] = plotSeparate( )

global data wndVr pltVr; %#ok<REDEF>

fh = findall(0,'type','figure');
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


for jj = 1:size(data.scombined,2)
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
    title(strrep(data.filename, '_', ' '));
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
end
end