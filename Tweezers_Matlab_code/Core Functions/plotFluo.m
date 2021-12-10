% Import global values used troughout the program.
global data fluodata fluof wndVr
% Useful function for searching in cell structures.
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
%% Draw Window
if ishandle(fluof)
    close(fluof)
end
if ishandle(wndVr.fluoPanel)
    close(wndVr.fluoPanel)
end
if ishandle(wndVr.axesPosFluo)
    close(wndVr.axesPosFluo)
end
if ishandle(wndVr.b)
    close(wndVr.b)
end

fluof = figure;

if size(wndVr.monPos,1) > 1
    set(fluof,'Position',wndVr.monPos(size(wndVr.monPos,1)-1,:));
end

drawnow;
pause(0.1);
jFrame = get(handle(fluof),'JavaFrame');
jFrame.setMaximized(true);
pause(0.1);
P2 = get(fluof,'pos');
%set(fluof,'DefaultFigurePosition', P2,'DefaultAxesPosition', [0.4,0.1,0.5,0.7]);
set(fluof,'pos', P2,'DefaultAxesPosition', [0.4,0.1,0.5,0.7]);
colors      =   linspecer(8);
set(fluof,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name','Fluorescence Data','Color','g');
set(fluof,'DefaultUicontrolUnits','normalized','DefaultUicontrolFontUnits','normalized',...
    'DefaultUicontrolFontSize',0.45,'DefaultUicontrolFontWeight','bold')
ax = gca;
set(ax,'handlevisibility','off','visible','off');
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
exportplot = -1;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Positions of text boxes in the main window

edge        =   0.01;
wid1        =   (0.6)*((1-9*edge)/4);
wid2        =   (0.6)*((1-9*edge)/4);
left1       =   2*edge;
left2       =   3*edge+wid1;
left3       =   6*edge+wid1+wid2;
left4       =   7*edge+2*wid1+wid2;
left5       =   10*edge+3*wid1+wid2;
left6       =   11*edge+4*wid1+wid2;
titlefs     =   0.5;
bot_title   =   0.835;

for i = 1 : 14, eval(['bot' num2str(i) '= bot_title - (i-1) * 0.06;']); end;
for i = 1 : 16, eval(['bbot' num2str(i) '= bot_title + 0.02 - (i-1) * 0.065;']); end;
height      =   0.045;
color       =   [0.7 0.7 0.7];
color1      =   colors(2,:);%[0.6 0.8 0.9];
color2      =   'w';
color3      =   colors(5,:);%[0.6550 0.75 0.4722];
color_x     =   [0.6 0.1 0.8];
color_y     =   [0.2 0.7 0.3];
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

set(0,'Currentfigure',fluof)
wndVr.fluoPanel = uipanel('BorderType','etchedin',...
		'BackgroundColor',colors(8,:),...
		'Position',[edge*3,0.65-edge*3,1-edge*6,0.35+edge*3],...
		'Parent',fluof);
    figBord = [0.05 0.2 0.9 0.6];
    wndVr.b = axes('parent',wndVr.fluoPanel,'Position', figBord, 'FontSize', 16);
    
wndVr.axesPosFluo = get(wndVr.fluoPanel,'pos');
%% Process Data

% Get inaccurate time offset from the Scanary file and Tweezers file by
% subtracting their metadata time.
fluodata.time = fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Time (ms)'), {fluodata.Data.MeasuredData.Name}), 1)).Data*1E-3;
fluodata.timeStamp = fluodata.Data.Root.Property(find(cellfun(cellfind('Metadata.Date/time'), {fluodata.Data.Root.Property.Name}))).Value{1};
data.timeStamp = data.Data.Root.Property(find(cellfun(cellfind('Date/time'), {data.Data.Root.Property.Name}))).Value{1};
fluodata.offset = (str2num(fluodata.timeStamp(end-5:end-4))-str2num(data.timeStamp(end-5:end-4)))*3600+(str2num(fluodata.timeStamp(end-3:end-2))-str2num(data.timeStamp(end-3:end-2)))*60+(str2num(fluodata.timeStamp(end-1:end))-str2num(data.timeStamp(end-1:end)));

% Load Time and Position Data from the Tweezers data file.
time = unique(fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Time (ms)'), {fluodata.Data.MeasuredData.Name}),1)).Data).*1E-3;
posx = unique(fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Set position X (um)'), {fluodata.Data.MeasuredData.Name}),1)).Data);

% Use WMIC to get accurate timestamp of times both the Scanary and Tweezers
% file were last written. WMIC format is as follows:
% yyyymmddHHMMSS.mmmmmmsUUU
% where:
% 
% yyyy Four-digit year (0000 through 9999).
% mm Two-digit month (01 through 12).
% dd Two-digit day of the month (01 through 31).
% HH Two-digit hour of the day using the 24-hour clock (00 through 23).
% MM Two-digit minute in the hour (00 through 59).
% SS Two-digit number of seconds in the minute (00 through 59).
% xxxxxx Six-digit number of microseconds in the second (000000 through 999999)
% s Plus sign (+) or minus sign (-) to indicate a positive or negative offset from Coordinated Universal Times (UTC).
% UUU Three-digit offset indicating the number of minutes that the originating time zone deviates from UTC.
data.pathnameWmic = strrep(data.pathname,'\\','+++');
data.pathnameWmic = strrep(data.pathnameWmic,'\','\\');
data.pathnameWmic = strrep(data.pathnameWmic,'+++','\\\');
fluodata.pathnameWmic = strrep(fluodata.pathname,'\\','+++');
fluodata.pathnameWmic = strrep(fluodata.pathnameWmic,'\','\\');
fluodata.pathnameWmic = strrep(fluodata.pathnameWmic,'+++','\\\');
[~,str] = dos(['wmic datafile where name="' data.pathnameWmic '" get lastmodified | findstr /brc:[0-9]']); %get CreationDate | findstr /brc:[0-9]
data.timeStampAcc = str;
data.timeStampAccstart = str2num(data.timeStampAcc(end-20:end-19))*3600+str2num(data.timeStampAcc(end-18:end-17))*60+str2num(data.timeStampAcc(end-16:end-15))+str2num(data.timeStampAcc(end-13:end-8))*1E-6;
[~,str] = dos(['wmic datafile where name="' fluodata.pathnameWmic '" get lastmodified | findstr /brc:[0-9]']);
fluodata.timeStampAcc = str;
fluodata.timeStampAccstart = str2num(fluodata.timeStampAcc(end-20:end-19))*3600+str2num(fluodata.timeStampAcc(end-18:end-17))*60+str2num(fluodata.timeStampAcc(end-16:end-15))+str2num(fluodata.timeStampAcc(end-13:end-8))*1E-6;
fluodata.offsetAcc = (fluodata.timeStampAccstart - data.timeStampAccstart) + ((data.time(end)-data.time(1))-(fluodata.time(end)-fluodata.time(1)));

% If selected, create shifted values for the R G B data
beadShift = 1;
undBndShift = 65;
%if beadShift == 1
clear extShiftU clear extShiftD
extShiftU = data.sbead2 - data.sbead1 - min(data.sbead2 - data.sbead1);%- data.sbead2 - min(-data.sbead2);%- data.sbead2 - min(-data.sbead2); %ones(1,length(- data.sbead2 - min(-data.sbead2)));%
extShiftD = data.sbead2 - data.sbead1 +data.straps - min(data.sbead2 - data.sbead1 +data.straps);%data.sbead2+data.straps - min(data.sbead2+data.straps);
posxDelta = mean(diff(posx));
clear idx
% idx = zeros(1,length(time));
% for ii = 1:length(time)
%     [idx(ii), idx(ii)] = min(abs((data.time-fluodata.offsetAcc)-(time(ii))));
% end

idx = nearestpoint(time,abs(data.time-fluodata.offsetAcc));
clear pixShiftU pixShiftUa pixShiftD pixShiftDa
pixShiftU = floor(extShiftU(idx)./posxDelta);
pixShiftUa = (extShiftU(idx)./posxDelta);
pixShiftD = floor(extShiftD(idx)./posxDelta);
pixShiftDa = (extShiftD(idx)./posxDelta);
% pixShiftU = floor(extShiftU(idx)./posxDelta);
% pixShiftUa = (extShiftU(idx)./posxDelta);
% pixShiftD = floor(extShiftD(idx)./posxDelta);
% pixShiftDa = (extShiftD(idx)./posxDelta);
%     for ii = 1:length(pixShift)
% %         imR.CData(1+pixShift(ii):end,ii) = imR.CData(1:end-pixShift(ii),ii);
% %         imR.CData(1:pixShift(ii),ii) = zeros(pixShift(ii),1);
% %         imG.CData(1+pixShift(ii):end,ii) = imG.CData(1:end-pixShift(ii),ii);
% %         imG.CData(1:pixShift(ii),ii) = zeros(pixShift(ii),1);
%         imR.CData(1:end-pixShift(ii),ii) = imR.CData(1+pixShift(ii):end,ii);
%         imR.CData(end-pixShift(ii)+1:end,ii) = zeros(pixShift(ii),1);
%         imG.CData(1:end-pixShift(ii),ii) = imG.CData(1+pixShift(ii):end,ii);
%         imG.CData(end-pixShift(ii)+1:end,ii) = zeros(pixShift(ii),1);
%     end
%end
% Calculate the RGB values from the Scanary data, these can be scaled using
% a polynomial scale bScal to improve contrast. Also accounts for two file
% layouts, one where only channel 1 is used (the scenario where only one
% color is installed in the CTRAP) and when all three channels are active
% (The current scenario where both the green and red color are installed on
% the CTRAP.
bScal = 1;
if isempty(find(cellfun(cellfind('Data/Pixel ch 2'), {fluodata.Data.MeasuredData.Name}), 1))
    cvalueg = reshape(fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Pixel ch 1'), {fluodata.Data.MeasuredData.Name}), 1)).Data,length(posx),length(time));
    cvaluegs = cvalueg;
    for ii = 1:length(pixShiftD)
        cvaluegs(undBndShift:end-pixShiftD(ii),ii) = cvaluegs(undBndShift+pixShiftD(ii):end,ii);
        cvaluegs(end-pixShiftD(ii)+1:end,ii) = zeros(pixShiftD(ii),1);
        cvaluegs(1:undBndShift-pixShiftU(ii),ii) = cvaluegs(pixShiftU(ii)+1:undBndShift,ii);
        %cvalueg(end-pixShiftD(ii)+1:end,ii) = zeros(pixShiftU(ii),1);
    end
    maxg = max(max(cvalueg));
    for ii = 1:max(max(cvalueg))+1
        greenMap(ii,1) = 0.0;
        greenMap(ii,2) = (ii-1).^bScal/maxg.^bScal;
        greenMap(ii,3) = 0.0;
    end
    colormap(greenMap)
    %contourf(time,posx,cvalue,15);
    colorbar
    cvaluegCor = ((cvalueg+1).^bScal/(maxg+1).^bScal);
    cvaluegsCor = ((cvaluegs+1).^bScal/(maxg+1).^bScal);
    set(0,'Currentfigure',fluof)
    imG = imagesc(time,posx,cvaluegCor);
    %imG = imagesc(time,posx,cvaluegCor);
    hold on;
    cvaluerCor = zeros(size(cvaluegCor,1),size(cvaluegCor,2));
    imR = imagesc(time,posx,cvaluerCor);
    cvaluebCor = zeros(size(cvaluegCor,1),size(cvaluegCor,2));
    imB = imagesc(time,posx,cvaluebCor);      
    hold off;
else
    cvaluer = reshape(fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Pixel ch 1'), {fluodata.Data.MeasuredData.Name}), 1)).Data,length(posx),length(time));
    cvaluers = cvaluer;
%     for ii = 1:length(pixShiftD)
%         cvaluers(undBndShift:end-pixShiftD(ii),ii) = cvaluers(undBndShift+pixShiftD(ii):end,ii);
%         cvaluers(end-pixShiftD(ii)+1:end,ii) = zeros(pixShiftD(ii),1);
%         cvaluers(1:undBndShift-pixShiftU(ii),ii) = cvaluers(pixShiftU(ii)+1:undBndShift,ii);
%     end
    for ii = 1:length(pixShiftD)
%         cvaluers(undBndShift:end-pixShiftD(ii),ii) = cvaluers(undBndShift+pixShiftD(ii):end,ii);
%         cvaluers(undBndShift:end-pixShiftD(ii),ii) = conv(cvaluers(undBndShift:end-pixShiftD(ii),ii),[pixShiftDa(ii)-pixShiftD(ii) 1-(pixShiftDa(ii)-pixShiftD(ii))], 'same');
%         cvaluers(end-pixShiftD(ii)+1:end,ii) = zeros(pixShiftD(ii),1);
%         cvaluers(1:undBndShift-pixShiftU(ii),ii) = cvaluers(pixShiftU(ii)+1:undBndShift,ii);
%         cvaluers(1:undBndShift-pixShiftU(ii),ii) = conv(cvaluers(1:undBndShift-pixShiftU(ii),ii),[pixShiftUa(ii)-pixShiftU(ii) 1-(pixShiftUa(ii)-pixShiftU(ii))], 'same');
        cvaluers(1:undBndShift,ii) = conv(cvaluers(1:undBndShift,ii),[1-(pixShiftUa(ii)-pixShiftU(ii)) pixShiftUa(ii)-pixShiftU(ii)], 'same');        
    cvaluers(1:undBndShift,ii) = circshift(cvaluers(1:undBndShift,ii),pixShiftU(ii));
    end
    cvalueg = reshape(fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Pixel ch 2'), {fluodata.Data.MeasuredData.Name}), 1)).Data,length(posx),length(time));
    cvaluegs = cvalueg;
    for ii = 1:length(pixShiftD)
        %cvaluegs(undBndShift:end-pixShiftD(ii),ii) = cvaluegs(undBndShift+pixShiftD(ii):end,ii);  
        %cvaluegs(:,ii) = circshift(cvaluegs(:,ii),pixShiftD(ii));
        %cvaluegs(undBndShift:end-pixShiftD(ii),ii) = conv(cvaluegs(undBndShift:end-pixShiftD(ii),ii),[pixShiftDa(ii)-pixShiftD(ii) 1-(pixShiftDa(ii)-pixShiftD(ii))], 'same');
        %cvaluegs(end-pixShiftD(ii)+1:end,ii) = zeros(pixShiftD(ii),1);
        %cvaluegs(1:undBndShift-pixShiftU(ii),ii) = cvaluegs(pixShiftU(ii)+1:undBndShift,ii);
        cvaluegs(1:undBndShift,ii) = conv(cvaluegs(1:undBndShift,ii),[1-(pixShiftUa(ii)-pixShiftU(ii)) pixShiftUa(ii)-pixShiftU(ii)], 'same');        
        cvaluegs(1:undBndShift,ii) = circshift(cvaluegs(1:undBndShift,ii),pixShiftU(ii));
        cvaluegs(undBndShift:end,ii) = conv(cvaluegs(undBndShift:end,ii),[pixShiftDa(ii)-pixShiftD(ii) 1-(pixShiftDa(ii)-pixShiftD(ii))], 'same');        
        cvaluegs(undBndShift:end,ii) = circshift(cvaluegs(undBndShift:end,ii),-pixShiftD(ii));
        
        
        %cvalueg(1:undBndShift-pixShiftU(ii),ii) = (pixShiftUa(ii)).*cvalueg(pixShiftU(ii)+1:undBndShift,ii) + (1-pixShiftUa(ii)).*cvalueg(pixShiftU(ii)+2:undBndShift+1,ii);
    end

    cvalueb = reshape(fluodata.Data.MeasuredData(find(cellfun(cellfind('Data/Pixel ch 3'), {fluodata.Data.MeasuredData.Name}), 1)).Data,length(posx),length(time));
    cvaluebs = cvalueb;
    for ii = 1:length(pixShiftD)
%         cvaluebs(undBndShift:end-pixShiftD(ii),ii) = cvaluebs(undBndShift+pixShiftD(ii):end,ii);
%         cvaluebs(undBndShift:end-pixShiftD(ii),ii) = conv(cvaluebs(undBndShift:end-pixShiftD(ii),ii),[pixShiftDa(ii)-pixShiftD(ii) 1-(pixShiftDa(ii)-pixShiftD(ii))], 'same');
%         cvaluebs(end-pixShiftD(ii)+1:end,ii) = zeros(pixShiftD(ii),1);
%         cvaluebs(1:undBndShift-pixShiftU(ii),ii) = cvaluebs(pixShiftU(ii)+1:undBndShift,ii);
%         cvaluebs(1:undBndShift-pixShiftU(ii),ii) = conv(cvaluebs(1:undBndShift-pixShiftU(ii),ii),[pixShiftUa(ii)-pixShiftU(ii) 1-(pixShiftUa(ii)-pixShiftU(ii))], 'same');
    cvaluebs(1:undBndShift,ii) = conv(cvaluebs(1:undBndShift,ii),[1-(pixShiftUa(ii)-pixShiftU(ii)) pixShiftUa(ii)-pixShiftU(ii)], 'same');        
    cvaluebs(1:undBndShift,ii) = circshift(cvaluebs(1:undBndShift,ii),pixShiftU(ii));
        
    end

    maxr = max(max(cvaluer));
    maxg = max(max(cvalueg));
    maxb = max(max(cvalueb));
    
    for ii = 1:max(max(cvaluer))+1
        redMap(ii,1) = (ii-1).^bScal/max(max(maxg,maxr),maxb).^bScal;
        redMap(ii,2) = 0.0;
        redMap(ii,3) = 0.0;
    end
    hold on;
    colormap(redMap)
    %colorbar
    set(0,'Currentfigure',fluof)
    cvaluerCor = ((cvaluer+1).^bScal/(max(max(maxg,maxr),maxb)+1).^bScal);
    %cvaluerCor = cvaluerCor./sum(cvaluerCor,1).*mean(sum(cvaluerCor,1));
    cvaluersCor = ((cvaluers+1).^bScal/(max(max(maxg,maxr),maxb)+1).^bScal);
    cvaluersCor = cvaluersCor./sum(cvaluersCor,1).*mean(sum(cvaluersCor,1));
    imR = imagesc(time,posx,cvaluerCor);
    
    for ii = 1:max(max(cvalueg))+1
        greenMap(ii,1) = 0.0;
        greenMap(ii,2) = (ii-1).^bScal/max(max(maxg,maxr),maxb).^bScal;
        greenMap(ii,3) = 0.0;
    end
    colormap(greenMap)
    %contourf(time,posx,cvalue,15);
    %colorbar
    set(0,'Currentfigure',fluof)
    cvaluegCor = ((cvalueg+1).^bScal/(max(max(maxg,maxr),maxb)+1).^bScal);
    %cvaluegCor = cvaluegCor./sum(cvaluegCor,1).*mean(sum(cvaluegCor,1));
    cvaluegsCor = ((cvaluegs+1).^bScal/(max(max(maxg,maxr),maxb)+1).^bScal);
    cvaluegsCor = cvaluegsCor./sum(cvaluegsCor,1).*mean(sum(cvaluegsCor,1));
    imG = imagesc(time,posx,cvaluegCor);
    
    for ii = 1:max(max(cvalueg))+1
        blueMap(ii,1) = 0.0;
        blueMap(ii,2) = (ii-1).^bScal/max(max(maxg,maxr),maxb).^bScal;
        blueMap(ii,3) = 0.0;
    end
    colormap(blueMap)
    %contourf(time,posx,cvalue,15);
    %colorbar
    set(0,'Currentfigure',fluof)
    cvaluebCor = ((cvalueb+1).^bScal/(max(max(maxg,maxr),maxb)+1).^bScal);
    %cvaluebCor = cvaluebCor./sum(cvaluebCor,1).*mean(sum(cvaluebCor,1));
    cvaluebsCor = ((cvaluebs+1).^bScal/(max(max(maxg,maxr),maxb)+1).^bScal);
    cvaluebsCor = cvaluebsCor./sum(cvaluebsCor,1).*mean(sum(cvaluebsCor,1));
    imB = imagesc(time,posx,cvaluebCor);
end
%Find Bead edges
sumcvalueTotal = [];

sumcvalueTotal = sum(cvaluersCor+cvaluegsCor+cvaluebsCor,2);
%% Find Upper and Lower bounds of beads

%figure;
clear pks locs peakwidth pkLoc pkIdx
%[pks,locs,peakwidth] = findpeaks(sumcvalueTotal(1:end-min(pixShiftD))-min(sumcvalueTotal(1:end-min(pixShiftD))),posx(1:end-min(pixShiftD)),'Annotate','extents','SortStr','descend','MinPeakWidth',0.1,'WidthReference','halfheight');%halfprom
peakSum = sumcvalueTotal(1:end-min(pixShiftD))-min(sumcvalueTotal(1:end-min(pixShiftD)));
%peakSum = peakSum-min(peakSum(30:end-30));
peakSum = peakSum-max(min(peakSum(1:20)),min(peakSum(end-20:end)));
peakLoc = posx(1:end-min(pixShiftD));
[locs, pks] = peakseek(peakSum,5);
XXXX = posx(1:end-min(pixShiftD));
locs = XXXX(locs);
[pks, sortI] = sort(pks,'descend');
locs = locs(sortI);
clear XXXX sortI;
if length(locs)>=4
    [pkLoc, pkIdx] = sort(locs(1:4));
else
    [pkLoc, pkIdx] = sort(locs);
end
pks = pks(pkIdx);
clear peakIdx;
pkhfVal = 2;
for ii = 1:length(pks)
    peakMin = peakLoc(find(peakSum<=max(pks(ii)/pkhfVal)));
    peakIdx(ii) = nearestpoint(pkLoc(ii),peakLoc(find(peakSum<=max(pks(ii)/pkhfVal))));
    peakIdx(ii) = peakMin(peakIdx(ii));
    clear peakMin;
end
clear peakwidth
peakwidth = abs(pkLoc-peakIdx');

% figure; plot(peakLoc,peakSum);
% hold on;
% plot(pkLoc,pks,'+');
% plot(peakIdx,pks/pkhfVal, '+');
% hold off;

%peakwidth = peakwidth(pkIdx);
%findpeaks(sumcvalueTotal(1:end-min(pixShiftD))-min(sumcvalueTotal(1:end-min(pixShiftD))),posx(1:end-min(pixShiftD)),'Annotate','extents','SortStr','descend','MinPeakWidth',0.1,'WidthReference','halfheight');%,'MinPeakWidth',0.4

% Calculate the time window in which there is Tweezers data in the
% fluorescence data. This assumes that the Tweezers file covers a larger
% time window than the Scanary file.
timeMin = find(data.time-fluodata.offsetAcc < min(time),1,'last');
if isempty(timeMin)
    timeMin = 1;
end
timeMax = find(data.time-fluodata.offsetAcc > max(time),1,'first');
if isempty(timeMax)
    timeMax = length(data.time);
end
% Determine the upper and lower bounds of the tether, incorporating the time shift between the Tweezers and Scanary data files.
fluodata.sOffset = 41.8;
%upBound = fluodata.sOffset-data.sbead2(timeMin:timeMax)+data.shiftX;
%lowBound = fluodata.sOffset+data.sbead2(timeMin:timeMax)+data.straps(timeMin:timeMax)+data.shiftX-1.7;

% extShiftU = data.sbead2 - data.sbead1 - min(data.sbead2 - data.sbead1);%- data.sbead2 - min(-data.sbead2);%- data.sbead2 - min(-data.sbead2); %ones(1,length(- data.sbead2 - min(-data.sbead2)));%
% extShiftD = data.sbead2 - data.sbead1 +data.straps - min(data.sbead2 - data.sbead1 +data.straps);%data.sbead2+data.straps - min(data.sbead2+data.straps);
%
upBound = pkLoc(2) + 1*peakwidth(2) - (data.sbead2(timeMin:timeMax) - data.sbead1(timeMin:timeMax) - min(data.sbead2(timeMin:timeMax) - data.sbead1(timeMin:timeMax)));%+0.5*peakwidth(2)
lowBound = pkLoc(3) - 1*peakwidth(3) + (data.sbead2(timeMin:timeMax) - data.sbead1(timeMin:timeMax) +data.straps(timeMin:timeMax) - min(data.sbead2(timeMin:timeMax) - data.sbead1(timeMin:timeMax) +data.straps(timeMin:timeMax)));
%
% upBound = pkLoc(2)+0.5*peakwidth(2)-data.sbead2(timeMin:timeMax)-min(-data.sbead2(timeMin:timeMax));%+0.5*peakwidth(2)
% lowBound =  pkLoc(3)-0.5*peakwidth(3)+data.sbead2(timeMin:timeMax)+data.straps(timeMin:timeMax)-min(data.sbead2(timeMin:timeMax)+data.straps(timeMin:timeMax));%-0.5*peakwidth(3)
%% Draw Figure
tweezTime = data.time(timeMin:timeMax)-fluodata.offsetAcc;

%ts1 = resample(ts, time) 
% imR.CData = imR.CData./sum(imG.CData,1).*size(cvaluers,1)./3;
% imG.CData = imG.CData./sum(imG.CData,1).*size(cvaluegs,1)./3;
% imB.CData = imB.CData./sum(imG.CData,1).*size(cvaluebs,1)./3;
% Create a true-color image from the compiled R G B data (colors that are
% not used are initialized to zero.
imRGB = cat(3, imR.CData, imG.CData, imB.CData);
% Redraw the axes to reset the image command, generating the previous RGB
% data causes them to be set to the wrong values
wndVr.fluoPanel = uipanel('BorderType','etchedin',...
		'BackgroundColor',colors(8,:),...
		'Position',[edge*3,0.65-edge*3,1-edge*6,0.35+edge*3],...
		'Parent',fluof);
    figBord = [0.05 0.2 0.9 0.6];
    wndVr.b = axes('parent',wndVr.fluoPanel,'Position', figBord, 'FontSize', 16);    
wndVr.axesPosFluo = get(wndVr.fluoPanel,'pos');
set(0,'Currentfigure',fluof)
%Plot the true-color image.
procimRGB = image(time,posx,imRGB);
%draw the boundaries.
hold on;
%plot(time, (fluodata.sOffset)*ones(1,length(time))+(data.scombined{:})','LineWidth',2,'Color','r');
plot(tweezTime, upBound,'LineWidth',1,'Color','b');
plot(tweezTime, lowBound,'LineWidth',1,'Color','r');
hold off;


sumcvalueCor = [];
fluotimeUnique = unique(fluodata.time);
sIndex = nearestpoint(fluotimeUnique,tweezTime);
cupBound = upBound(sIndex);
supBound = nearestpoint(cupBound,posx);
clowBound = lowBound(sIndex);
slowBound = nearestpoint(clowBound,posx);

for ii = 1:length(fluotimeUnique)
%[ss(ii) sIndex(ii)] = min(abs(tweezTime-fluotimeUnique(ii)));
%cupBound(ii) = upBound(sIndex(ii));
%[pp(ii) supBound(ii)] = min(abs(posx-cupBound(ii)));
%clowBound(ii) = lowBound(sIndex(ii));
%[qq(ii) slowBound(ii)] = min(abs(posx-clowBound(ii)));
sumcvalueCor(ii) = sum(cvaluerCor(supBound(ii):slowBound(ii),ii)+cvaluegCor(supBound(ii):slowBound(ii),ii)+cvaluebCor(supBound(ii):slowBound(ii),ii))./length(supBound(ii):slowBound(ii));
end

range = 57:108;
%figure; plot(time+fluodata.offsetAcc,(sum(cvaluegCor(range,:)./mean(cvaluegCor,1),1)-min(sum(cvaluegCor(range,:)./mean(cvaluegCor,1),1)))./max(sum(cvaluegCor(range,:)./mean(cvaluegCor,1),1)-min(sum(cvaluegCor(range,:)./mean(cvaluegCor,1),1))),data.time(pltVr.lIdx(1):pltVr.rIdx(1)), abs(((data.scombined{1}-min(data.scombined{1}))./max(data.scombined{1}-min(data.scombined{1})))))
%figure; plot(posx,cvaluegCor(1,:));

%figure; plot(unique(fluodata.time),sumcvalueCor);
sumcvalueStretch = [];
%276:538; 1250:1437; 937:1312; 

% pixStretch = max(slowBound(1:length(fluotimeUnique))-supBound(1:length(fluotimeUnique)));
% for ii = 1:length(fluotimeUnique)
% sumcvalueStretch(1:pixStretch,ii) = (imresize(cvaluer(supBound(ii):slowBound(ii),ii), [pixStretch 1], 'nearest')+imresize(cvalueg(supBound(ii):slowBound(ii),ii), [pixStretch 1], 'nearest')+imresize(cvalueb(supBound(ii):slowBound(ii),ii), [pixStretch 1], 'nearest'));%./length(supBound(ii):slowBound(ii));
% end
% baseLine = sum(sumcvalueStretch(:,1:300),2)/300;
% %sumcvalueStretch = sumcvalueStretch-baseLine;
% figure;
% 
% baseCor = sum(sumcvalueStretch,2)/max(sum(sumcvalueStretch,2));
% datCor = sumcvalueStretch./max(max(sumcvalueStretch,2));
% CCC = bsxfun(@minus,datCor,baseCor);
% surf(fluotimeUnique,posxDelta:posxDelta:pixStretch*posxDelta,CCC,'EdgeColor','none');
% %surf(fluotimeUnique,posxDelta:posxDelta:pixStretch*posxDelta,sumcvalueStretch./max(max(sumcvalueStretch)),'EdgeColor','none');
% colormap(greenMap)