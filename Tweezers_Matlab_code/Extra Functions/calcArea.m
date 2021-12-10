clc;

fh = findall(0,'type','figure');
for ii = 1:length(fh)
    if fh(ii) ~= 1
        close(fh(ii))
    end
end
f = figure(2);
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
pause(0.1);
P2 = get(gcf,'pos');
set(0,'DefaultFigurePosition', P2,'DefaultAxesPosition', [0.1,0.1,0.8,0.8]);

%extL = 0
aFit = 70;
maxlim1 = 1;

C = linspecer(size(data.scombined,2)*2);
for jj = 1:size(data.scombined,2)
    f = figure(2);
    [i,j] = max(data.force{jj}+data.shiftY);
    xShift = mean(data.scombined{jj}(1:j)+data.shiftX);
    tStep = (data.scombined{jj}(j)-data.scombined{jj}(1))./length(data.scombined{jj}(1:j));
    %xP = [[data.scombined{jj}(1)-extL*tStep:tStep:data{1}.scombined(1)-tStep]'; data.scombined{jj}(1:j); [data.scombined{jj}(j)+tStep:tStep:data.scombined{jj}(j)+extL*tStep]']+data.shiftX;
    %yP = [data.force{jj}(1).*ones(extL,1); data.force{jj}(1:j); data.force{jj}(j).*ones(extL,1)]+data.shiftY;
    xPf = data.scombined{jj}(1:j)+data.shiftX;
    [xPf,I] = sort(xPf);
    yP = data.force{jj}(1:j)+data.shiftY;
    yP = yP(I);
    
    
    
    xPf = xPf-xShift;
    plot(data.scombined{jj}+data.shiftX, data.force{jj}+data.shiftY,'color',C(jj*2-1,:));    
    hold on
    %Pfit1 = polyfit(xPf,yP,aFit);
    %Pdev{jj} =  polyval(Pfit1,xPf);
    f1 = fit(xPf+xShift, yP,  'smoothingspline', 'SmoothingParam', 0.9999999998982048);
    plot(xPf+xShift, f1(xPf+xShift),'color',C(jj*2,:),'linewidth',2);
    
    xPb = data.scombined{jj}(j:end)+data.shiftX;
    [xPb,I] = sort(xPb);
    yP = data.force{jj}(j:end)+data.shiftY;
    yP = yP(I);    
    
    xPb = xPb-xShift;
    %plot(xPb+xShift, yP,'color',C(jj*2-1,:));    
    hold on
    %Pfit2 = polyfit(xPb,yP,aFit);
    %Pdev{jj} =  polyval(Pfit2,xPb);
    f2 = fit(xPb+xShift, yP,  'smoothingspline', 'SmoothingParam', 0.9999999998982048);
    plot(xPb+xShift, f2(xPb+xShift),'color',C(jj*2,:),'linewidth',2);
    %hold off;
    
    f = figure(3);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    pause(0.1);
    P2 = get(gcf,'pos');
    set(0,'DefaultFigurePosition', P2,'DefaultAxesPosition', [0.1,0.1,0.8,0.8]);
    hold on;
    maxlim1 = max(maxlim1,max(f1(xPf+xShift)-f2(xPf+xShift)));
    plot(xPf+xShift,f1(xPf+xShift)-f2(xPf+xShift),'color',C(jj*2,:),'linewidth',2);
    ylim([-1 maxlim1*1.1]) 
    set(gca,'Fontsize',16)
    title(data.filename);
    set(f,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name',data.filename);
    set(f,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
        'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
    xlabel('Extension (µm)')
    ylabel('Force (pN)');
    %set(gca,'ylim',[-1 YL(2)+50]) % Adjust only the lower limit.
    
    %hold off;
    %plot(data.scombined{jj}(j:end)+data.shiftX, data.force{jj}(j:end)+data.shiftY,'color',C(jj*2,:));
    f = figure(4);
    drawnow;
    pause(0.1);
    jFrame = get(handle(gcf),'JavaFrame');
    jFrame.setMaximized(true);
    pause(0.1);
    P2 = get(gcf,'pos');
    set(0,'DefaultFigurePosition', P2,'DefaultAxesPosition', [0.1,0.1,0.8,0.8]);
    hold on;
    intX = fliplr((xPf+xShift)')';
    intY = -cumtrapz(fliplr((xPf+xShift)')',fliplr((f1(xPf+xShift)-f2(xPf+xShift))')');
    
    %plot(((xPf+xShift)'),-fliplr((cumtrapz(fliplr((xPf+xShift)')',fliplr((f1(xPf+xShift)-f2(xPf+xShift)))')')'),'color',C(jj*2,:),'linewidth',2);
    plot(intX,intY,'color',C(jj*2,:),'linewidth',2);
    set(gca,'Fontsize',16)
    title(data.filename);
    set(f,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name',data.filename);
    set(f,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
        'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
    xlabel('Extension (µm)')
    ylabel('Energy (pN*µm)');
    [intI, intJ] = min(abs(intX-0.9));
    disp(intY(intJ))
    
end

f = figure(2);
hold on;

set(gca,'Fontsize',16)
title(data.filename);
set(f,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name',data.filename);
set(f,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
    'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
xlabel('Extension (µm)')
ylabel('Force (pN)');

if rulVal == 1
    
elseif rulVal == 2
    
    load Luciferasewlc.mat %Luciferase ruler to get the unfolding length
    plot(ext, f0, 'k');
    
elseif rulVal == 3
    
    load sMBPrulerMinimal.mat %MBP ruler to get the unfolding length
    plot( rulerext,  rulerf, 'k')
    
    load wormlikechains % wormlike chains
    plot( ext,  singleMBPForce,'r', ext, DNAforce, 'k')
    
elseif rulVal == 4
    
    load Luciferasewlc.mat %Luciferase ruler to get the unfolding length
    plot( intermediates(:,1), intermediates(:,2), 'k', ext, OneLuci,'r', ext, f0, 'k');
    
end
xlim([str2num(get(sampling_xlim,'String'))]);
ylim([str2num(get(sampling_ylim,'String'))]);
hold off

