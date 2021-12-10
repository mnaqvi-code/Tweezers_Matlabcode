function [ data ] = plotFECombined( varargin ) % draws force-extension graphs for all chosen files
%plotFE - plots Force vs. Extension data for selected intervals.
%   plotFE( data, p_CCD, p_QPD_piezo, rulVal, datVr.sigVal, forVal, prgVr.forStt, fitVal, fitPar, index )

global p_CCD p_QPD_piezo f data pltVr datVr prgVr wndVr; %#ok<REDEF>

% if nargin == 1 %number of arguments
%     index = 1:size(data,2);
% end

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

pltVr.grphAct = 3;

edge = 0.01;
arrayfun(@cla,f)%findall(0,'type','axes')

%cla(findall(0,'type','axes'))

for ii = 1:size(data,2)
    if isfield(data,'figs')
        if isfield(data.figs,'fext')
            data.figs = rmfield(data.figs,'fext');
        end
    end
end
if ~isempty(f)
    set(0,'Currentfigure',f)
end
if isfield(data,'scombinedSplit') && isfield(data,'forceSplit')
    rmfield(data,{'scombinedSplit', 'forceSplit'});
end

if prgVr.setupStt == 0
    datVr.forVal = data.kx;
    
elseif prgVr.setupStt == 1
    data.time = data.data(:,1);
    
elseif prgVr.setupStt == 2
    if strcmp(data.Data.Root.Property(2).Value{1},'5')
        data.time = data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Time (ms)'), {data.Data.MeasuredData.Name}), 1)).Data*1E-3;
    elseif strcmp(data.Data.Root.Property(2).Value{1},'6')
        data.time = data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Time (ms)'), {data.Data.MeasuredData.Name}), 1)).Data*1E-3;
    end
end

for ii = 1:length(data.time)-1
    if data.time(ii+1) < data.time(ii)
        data.time(ii+1:end) = data.time(ii+1:end) + data.time(ii) + mean(diff(data.time));
    end
end

data.time = decimate(data.time,datVr.decFact);

for ii = 1:size(data,2)
    
    
    if ~isfield(data, 'x')
        data.x(1) = min(data.time);
        data.x(2) = max(data.time);
    elseif isfield(data, 'x')
        if isempty(data.x)
            data.x(1) = min(data.time);
            data.x(2) = max(data.time);
        end
    end
    
    if pltVr.legOn == 1 && pltVr.splitPulls == 1
        pltVr.legArry = cell(length(data.x),1);
        C = linspecer(size(1:length(data.x),2),'qualitative');
    elseif pltVr.legOn == 1 && pltVr.splitPulls ~= 1
        pltVr.legArry = cell(length(data.x)/2,1);
        C = linspecer(size(1:length(data.x)/2,2),'qualitative');
    end
    for jj =  1:length(data.x)/2
        [pltVr.lIdx(jj),pltVr.lIdx(jj)] = min(abs(data.x(2*jj-1)-data.time));
        [pltVr.rIdx(jj),pltVr.rIdx(jj)] = min(abs(data.x(2*jj)-data.time));
        if prgVr.setupStt == 0
            
            
            if prgVr.forStt == 1 %Fsmooth
                
                [X,I] = sort(data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                Y = data.xspt2(pltVr.lIdx(jj):pltVr.rIdx(jj))*0.0615;
                Y = Y(I);
                if length(pltVr.fitVal) == 2
                    J = intersect(find(Y >= (min(Y) + pltVr.fitVal(1)*(max(Y) - min(Y)))), find(Y <= (min(Y) + pltVr.fitVal(2)*(max(Y) - min(Y)))));
                    pp = splinefit(X(J),Y(J),5,4);
                else
                    pp = splinefit(X,Y,5,4);
                end
                
                [X,I] = sort(data.vx);
                Pdev{jj} = ppval(pp,X);
                unsorted = 1:size(data.vx,1);
                newInd(I) = unsorted;
                Pdev{jj} = Pdev{jj}(newInd);
                
                extension = -(data.xpz(pltVr.lIdx(jj):pltVr.rIdx(jj))-Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                
                if datVr.forVal == 0
                    data.force{jj} = data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj))*430.5./(-Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj))+6.52);
                else
                    data.force{jj} = -Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)).*datVr.forVal;
                end
                if length(pltVr.fitVal) == 1
                    Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)) = -pltVr.fitVal;
                    extension = -(data.xpz(pltVr.lIdx(jj):pltVr.rIdx(jj))+data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj))./-Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                    data.force{jj} = data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj))*datVr.forVal./(-Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                end
                
            elseif prgVr.forStt == 2 % Flinear
                [X,I] = sort(data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                Y = data.xspt2(pltVr.lIdx(jj):pltVr.rIdx(jj))*0.0615;
                Y = Y(I);
                if length(pltVr.fitVal) == 2
                    J = intersect(find(Y >= (min(Y) + pltVr.fitVal(1)*(max(Y) - min(Y)))), find(Y <= (min(Y) + pltVr.fitVal(2)*(max(Y) - min(Y)))));
                    Pfit = polyfit(Y(J),X(J),1);
                else
                    Pfit = polyfit(Y,X,1);
                end
                
                [X,I] = sort(data.vx);
                Y = data.xspt2*0.0615;
                Y = Y(I);
                Pdev{jj} =  polyval(polyder(Pfit),Y);
                if length(pltVr.fitVal) == 1
                    Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)) = -pltVr.fitVal;
                end
                extension = -(data.xpz(pltVr.lIdx(jj):pltVr.rIdx(jj))+data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj))./-Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                data.force{jj} = data.vx(pltVr.lIdx(jj):pltVr.rIdx(jj))*datVr.forVal./(-Pdev{jj}(pltVr.lIdx(jj):pltVr.rIdx(jj)));
                
                
                
            end
            
            % Now make W function
            ls = length(data.force{jj});
            lp = length(p_CCD);
            fss = 50;
            fsp = 50;
            p1 = p_CCD;
            p2 = p_QPD_piezo;
            a = 5;
            W = zeros(1, ceil(ls/2+1));
            for i = 1:lp
                floor(i*ceil(ls/2+1)*fsp/fss/lp);
                if floor(i*ceil(ls/2+1)*fsp/fss/lp) <= ceil(ls/2)-1
                    W(floor(i*ceil(ls/2+1)*fsp/fss/lp)+1) = 1./(exp(a)-1).*((exp(a)+1)./(1+exp(a*(p2(i)-p1(i))/(p1(i)+p2(i))))-1); %1/2*(1+1/exp(1)*sign(p1(i)-p2(i))*exp(abs(p1(i)-p2(i))/(p1(i)+p2(i))));
                end
            end
            
            %remove zeros by linear interpolation
            j = 0;
            if W(1) == 0
                while(W(1+j)==0)
                    j = j+1;
                end
                for m = 1:j
                    W(m) = W(j+1);
                end
            end
            
            j=0;
            
            if W(ceil(ls/2+1)) == 0
                while(W(ceil(ls/2+1)-j)==0)
                    j = j+1;
                end
                for m = 1:j
                    W(ceil(ls/2+1)+1-m) = W(ceil(ls/2+1)-j);
                end
            end
            
            for i = 1:ceil(ls/2+1)
                j = 0;
                if W(i) == 0
                    while(W(i+j)==0)
                        j = j+1;
                    end
                    for m = 0:j-1
                        W(i+m) = (W(i+j)-W(i-1))*(m)/(j)+W(i-1);
                    end
                end
            end
            
            % Now combine signals
            
            extensionvideo = (data.xspt2(pltVr.lIdx(jj):pltVr.rIdx(jj)) - data.xspt1(pltVr.lIdx(jj):pltVr.rIdx(jj))).*0.0615;
            s1 = extensionvideo;
            s2 = extension;
            
            ms1 = s1 - mean(s1);
            ms2 = s2 - mean(s2);
            
            fs1 = fft(ms1);
            fs2 = fft(ms2);
            
            if datVr.sigVal == 1
                W = ones(1,length(W));
            elseif datVr.sigVal == 2
                W = zeros(1,length(W));
            end
            
            % two cases for odd and even length
            
            if mod(ls,2) == 0
                for i = 2:ls/2
                    fs1(i) = (1-W(i))*fs1(i);
                    fs1(ls-i+2) = (1-W(i))*fs1(ls-i+2);
                    fs2(i) = W(i)*fs2(i);
                    fs2(ls-i+2) = W(i)*fs2(ls-i+2);
                end
                fs1(1) = (1-W(1))*fs1(1);
                fs1(ls/2+1) = (1-W(ls/2+1))*fs1(ls/2+1);
                fs2(1) = W(1)*fs2(1);
                fs2(ls/2+1) = W(ls/2+1)*fs2(ls/2+1);
            else
                for i = 2:(ls+1)/2
                    fs1(i) = (1-W(i))*fs1(i);
                    fs1(ls-i+2) = (1-W(i))*fs1(ls-i+2);
                    fs2(i) = W(i)*fs2(i);
                    fs2(ls-i+2) = W(i)*fs2(ls-i+2);
                end
                fs1(1) = (1-W(1))*fs1(1);
                fs2(1) = W(1)*fs2(1);
            end
            
            % finally reverse fourier transform
            fscombined = fs1 + fs2;
            data.scombined{jj} = ifft(fscombined) + mean(s1);
            
        elseif prgVr.setupStt == 1
            
            %             datVr.dualCst.gam1 = 6*pi*datVr.dualCst.eta*datVr.dualCst.R1; %N*s/m
            %             datVr.dualCst.k1 = 2*pi*datVr.dualCst.gam1*datVr.dualCst.cFreq1*10^6; %pN/um
            %             datVr.dualCst.Dth1 = datVr.dualCst.kB*datVr.dualCst.T/datVr.dualCst.gam1; %m^2/s
            %             datVr.dualCst.alpha1 = sqrt(datVr.dualCst.Dth1/datVr.dualCst.Dex1)*10^6; %um/a.u.
            %
            %             datVr.dualCst.gam2 = 6*pi*datVr.dualCst.eta*datVr.dualCst.R2; %N*s/m
            %             datVr.dualCst.k2 = 2*pi*datVr.dualCst.gam2*datVr.dualCst.cFreq2*10^6; %pN/um
            %             datVr.dualCst.Dth2 = datVr.dualCst.kB*datVr.dualCst.T/datVr.dualCst.gam2; %m^2/s
            %             datVr.dualCst.alpha2 = sqrt(datVr.dualCst.Dth2/datVr.dualCst.Dex2)*10^6; %um/a.u.
            
            %Mirror xCoordinates with the correct correction factor (Check this factor experimentally one more time later?).
            xMir = data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),11).*0.00016575398;%0.00016575398;%0.00016575398;%$0.00016575398;%0.000145;%0.00016575398;%0.00016575398 0.000123697
            %[Y,I] = sort(xMir);
            %xValU = unique(Y);
            %xPSD values.
            
            if ~exist('xPSD1Xmin')
                xPSD1Xmin = min(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),3));
            end
            if ~exist('xPSD1Ymin')
                xPSD1Ymin = max(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),5));
            end
            if ~exist('xPSD2Xmin')
                xPSD2Xmin = max(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),7));
            end
            if ~exist('xPSD2Ymin')
                xPSD2Ymin = max(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),9));
            end
            
            xPSD1X = (data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),3)-xPSD1Xmin).*datVr.dualCst.alpha1;
            xPSD1Y = -(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),5)-xPSD1Ymin).*datVr.dualCst.alpha3;
            xPSD1 = sqrt(xPSD1X.^2);
            %xPSD1 = sqrt(xPSD1X.^2+xPSD1Y.^2);
            xPSD1 = xPSD1 - min(xPSD1);
            xPSD2X = -(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),7)-xPSD2Xmin).*datVr.dualCst.alpha2;
            xPSD2Y = -(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),9)-xPSD2Ymin).*datVr.dualCst.alpha4;
            xPSD2 = sqrt(xPSD2X.^2);
            %xPSD2 = sqrt(xPSD2X.^2+xPSD2Y.^2);
            xPSD2 = xPSD2 - min(xPSD2);
           
            
            if ~exist('yDat1Xmin')
                yDat1Xmin = min(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),3));
            end
            if ~exist('yDat1Ymin')
                yDat1Ymin = max(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),5));
            end
            if ~exist('yDat2Xmin')
                yDat2Xmin = max(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),7));
            end
            if ~exist('yDat2Ymin')
                yDat2Ymin = max(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),9));
            end
            
            %yDat1 and yDat2 are swapped because approach is now from the right.
            yDat1X = (data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),3)-yDat1Xmin).*datVr.dualCst.k1.*datVr.dualCst.alpha1;
            yDat1Y = -(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),5)-yDat1Ymin).*datVr.dualCst.k3.*datVr.dualCst.alpha3;
            yDat1 = sqrt(yDat1X.^2);
            %yDat1 = sqrt(yDat1X.^2+yDat1Y.^2);
            yDat2X = -(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),7)-yDat2Xmin).*datVr.dualCst.k2.*datVr.dualCst.alpha2;
            yDat2Y = -(data.data(pltVr.lIdx(jj):pltVr.rIdx(jj),9)-yDat2Ymin).*datVr.dualCst.k4.*datVr.dualCst.alpha4;
            yDat2 = sqrt(yDat2X.^2);
            %yDat2 = sqrt(yDat2X.^2+yDat2Y.^2);
            
            if ~exist('yDat')
                minyDat = min((yDat1+yDat2)/2);
            end
            
            %keyboard;
            yDat = (yDat1+yDat2)/2;
            %minyDat
            data.force{jj} = smooth(yDat - minyDat,1);%,'sgolay',5);
            %xPSD1 = -(xPSD1-max(xPSD1));
            %xPSD2 = -(xPSD2-min(xPSD2));
            
            xDat = xMir - (xPSD2 + xPSD1);
            data.scombined{jj} = smooth(xDat +datVr.dualCst.R1*1E6 + datVr.dualCst.R2*1E6,1);%,'sgolay',5);
            
            if isfield(data,'beadData')
                %%interpolate PSD data to cameraData
                data.dataInt{jj} = [];
                data.range{jj} = [];
                data.range{jj} = find(data.beadData(:,1)>data.time(pltVr.lIdx(jj)),1):find(data.beadData(:,1)<data.time(pltVr.rIdx(jj)),1,'last');
                data.dataInt{jj}(:,1) = interp1( data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), data.beadData(data.range{jj},1 ));
                data.dataInt{jj}(:,2) = interp1( data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), data.force{jj}, data.dataInt{jj}(:,1) ); %% PSD Force
                data.dataInt{jj}(:,3) = interp1( data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), data.scombined{jj}, data.dataInt{jj}(:,1) ); %%PSD Extension
                data.dataInt{jj}(:,4) = interp1( data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), xPSD1, data.dataInt{jj}(:,1) ); %%PSD1
                data.dataInt{jj}(:,5) = interp1( data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), xPSD2, data.dataInt{jj}(:,1) ); %%PSD2
                
                data.dataInt{jj}(:,6) = data.beadData(data.range{jj},5); %%Camera1 x-coordinate
                data.dataInt{jj}(:,7) = data.beadData(data.range{jj},3); %%Camera2 x-coordinate
            end
            
        elseif prgVr.setupStt == 2
            pltVr.lIdx(jj) = pltVr.lIdx(jj)*datVr.decFact-datVr.decFact+1;
            pltVr.rIdx(jj) = pltVr.rIdx(jj)*datVr.decFact-datVr.decFact+1;
            %data.force{jj} = data.Data.MeasuredData(3).Data(pltVr.lIdx(jj):pltVr.rIdx(jj));
            if strcmp(data.Data.Root.Property(2).Value{1},'5')
                data.force{jj} = (data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj)))./2;
                data.scombined{jj} = data.Data.MeasuredData(8).Data(pltVr.lIdx(jj):pltVr.rIdx(jj));
            elseif strcmp(data.Data.Root.Property(2).Value{1},'6')
                data.force{jj} = sqrt(((data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value.*(datVr.dualCst.k1.*datVr.dualCst.alpha1)-data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value.*(datVr.dualCst.k2.*datVr.dualCst.alpha2))./2).^2 + ((data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 3 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 3 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value.*(datVr.dualCst.k4.*datVr.dualCst.alpha4)-data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 1 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 1 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value.*(datVr.dualCst.k3.*datVr.dualCst.alpha3))./2).^2);
                if ~isempty(find(cellfun(cellfind('FD Data/Trap X Position (um)'), {data.Data.MeasuredData.Name}), 1))
                    %data.scombined{jj} = data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Trap X Position (um)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./datVr.dualCst.k1+data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./datVr.dualCst.k2;
                    data.sbead1 = ((data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value+data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(2).Value).*(datVr.dualCst.alpha1));
                    data.sbead2 = ((data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value+data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(2).Value).*(datVr.dualCst.alpha2));
                    data.straps = data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Distance 1 (um)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj));

                    data.scombined{jj} = data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Distance 1 (um)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-((data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value+data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 0 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(2).Value).*(datVr.dualCst.alpha1))+((data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(1).Value+data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Force Channel 2 (pN)'), {data.Data.MeasuredData.Name}), 1)).Property(2).Value).*(datVr.dualCst.alpha2));
                else
                    data.scombined{jj} = data.Data.MeasuredData(find(cellfun(cellfind('FD Data/Distance 1 (um)'), {data.Data.MeasuredData.Name}), 1)).Data(pltVr.lIdx(jj):pltVr.rIdx(jj));
                end
            end
            pltVr.lIdx(jj) = (pltVr.lIdx(jj)+datVr.decFact-1)/datVr.decFact;
            pltVr.rIdx(jj) = (pltVr.rIdx(jj)+datVr.decFact-1)/datVr.decFact;
            
            data.force{jj} = decimate(data.force{jj},datVr.decFact);
            data.scombined{jj} = decimate(data.scombined{jj},datVr.decFact);
            
            
            
        end
        %         if pltVr.splitPulls == 1
        %             [valtemp,Itemp] = max(data.scombined{jj});
        %             if Itemp == length(data.scombined{jj})
        %                 Itemp = Itemp - 1;
        %             end
        %             data.scombinedSplit{1,jj} = data.scombined{jj}(1:Itemp);
        %             data.scombinedSplit{2,jj} = data.scombined{jj}(Itemp+1:end);
        %             data.forceSplit{1,jj} = data.force{jj}(1:Itemp);
        %             data.forceSplit{2,jj} = data.force{jj}(Itemp+1:end);
        %         end
        
        %Calculate Bias Correction
        
        if ~isempty(datVr.xBias) && ~isempty(datVr.yBias)
            biasP = polyfit(datVr.xBias,datVr.yBias,1);
            biasFit = biasP(1).*data.scombined{jj};
            data.force{jj} = data.force{jj}-biasFit;
        end
        
    end
    
    
    
    if prgVr.setupStt == 0
        minTemp = min(cellfun(@(x)min(x(:)), data.force));
        for jj =  1:length(data.x)/2
            data.force{jj} = data.force{jj} - minTemp;
        end
    end
    
    
    if pltVr.splitPulls == 1
        data.scombinedSplit = [];
        data.forceSplit = [];
        for jj =  1:length(data.x)/2
            [valtemp,Itemp] = max(data.scombined{jj});
            if Itemp == length(data.scombined{jj})
                Itemp = Itemp - 1;
            end
            data.scombinedSplit{1,jj} = data.scombined{jj}(1:Itemp);
            data.scombinedSplit{2,jj} = data.scombined{jj}(Itemp+1:end);
            data.forceSplit{1,jj} = data.force{jj}(1:Itemp);
            data.forceSplit{2,jj} = data.force{jj}(Itemp+1:end);
            
        end
    else
        data.scombinedSplit = [];
        data.forceSplit = [];
        data.forceSplit = data.force;
        data.scombinedSplit = data.scombined;
    end
    
    if(wndVr.hTrig == 0)
        wndVr.hTrig = 1;
        maxTemp = max(cellfun(@(x)max(x(:)), data.scombined));
        minTemp = min(cellfun(@(x)min(x(:)), data.scombined));
        set(wndVr.hSlider, 'Min', min(1*(0-maxTemp), 1*(0-minTemp)), 'Max', max(1*(0-minTemp)+0.9, 1*(0-maxTemp)+0.9), 'Value', mean([min(1*(0-maxTemp)+0.9, 1*(0-minTemp)+0.9) max(1*(0-minTemp)+0.9, 1*(0-maxTemp)+0.9)]), 'SliderStep', abs([1 1].*((1*(0-maxTemp))-(1*(0-minTemp)))./1000));
        data.shiftX = get(wndVr.hSlider,'Value');
        data.shiftY = get(wndVr.vSlider,'Value');
        set(wndVr.sampling_shift,'String', num2str([data.shiftX data.shiftY]));
    end
    
    if(wndVr.vTrig == 0)
        wndVr.vTrig = 1;
        maxTemp = max(cellfun(@(x)max(x(:)), data.force));
        minTemp = min(cellfun(@(x)min(x(:)), data.force));
        set(wndVr.vSlider, 'Min', min(1*(0-maxTemp)-0.9, 1*(0-minTemp)-0.9), 'Max', max(1*(0-minTemp)+0.9, 1*(0-maxTemp)+0.9), 'Value', max([min(1*(0-maxTemp), 1*(0-minTemp)) max(1*(0-minTemp), 1*(0-maxTemp))]), 'SliderStep', abs([1 1].*((1*(0-maxTemp))-(1*(0-minTemp)))./1000));
        data.shiftX = get(wndVr.hSlider,'Value');
        data.shiftY = get(wndVr.vSlider,'Value');
        set(wndVr.sampling_shift,'String', num2str([data.shiftX data.shiftY]));
    end
    
    datVr.fePoints = [];
    data.figs.feUnf.axh = []; %#ok<*STRNU>
    
    for jj =  1:length(data.x)/2
        hold on
        if pltVr.splitPulls ~= 1
            data.figs.fext.axh(jj) = plot(data.scombined{jj}+data.shiftX, data.force{jj}+data.shiftY,'color',C(jj,:), 'XDataSource',strcat('data{',num2str(ii),'}.scombined'),'YDataSource', strcat('data{',num2str(ii),'}.force'));
            %data.figs.fext.axh(jj) = plot(((data.dataInt{jj}(:,6)-data.dataInt{jj}(:,7))-min(data.dataInt{jj}(:,6)-data.dataInt{jj}(:,7))).*0.1+0.6+data.shiftX,data.dataInt{jj}(:,2)+data.shiftY,'color',C(jj,:), 'XDataSource',strcat('data{',num2str(ii),'}.scombined'),'YDataSource', strcat('data{',num2str(ii),'}.force'));
            title(strrep(data.filename, '_', ' '));
            xlabel('Extension (µm)')
            ylabel('Force (pN)');
            set(gca,'Tag',num2str(ii))
            hold off
            graphArray = cell(length(data.x)/2,1);
            pltVr.legArry{jj,1} = ['Cycle ' num2str(jj)];
        else
            for gg = 1:2
                data.figs.fext.axh(jj*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg)) = plot(data.scombinedSplit{gg,jj}+data.shiftX, data.forceSplit{gg,jj}+data.shiftY,'color',C(jj*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:), 'XDataSource',strcat('data{',num2str(ii),'}.scombinedSplit'),'YDataSource', strcat('data{',num2str(ii),'}.forceSplit'));
                %keyboard;
                %data.figs.fext.axh(jj*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg)) = plot(data.scombinedSplit{gg,jj}+data.shiftX, data.forceFit{gg,jj}(data.scombinedSplit{gg,jj}+data.shiftX-data.forceShift{gg,jj}),'color',C(jj*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:), 'XDataSource',strcat('data{',num2str(ii),'}.scombinedSplit'),'YDataSource', strcat('data{',num2str(ii),'}.forceSplit'));

                
                if prgVr.setupStt == 0
                    title(strrep(data.filename, '_', ' '));
                elseif prgVr.setupStt == 2
                    title(data.Data.Root.Property(1).Value{1});
                end
                xlabel('Extension (µm)')
                ylabel('Force (pN)');
                set(gca,'Tag',num2str(ii))
                
                graphArray = cell(length(data.x)/2,1);
                if gg == 1
                    pltVr.legArry{jj*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),1} = ['Pull ' num2str(jj)];
                else
                    pltVr.legArry{jj*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),1} = ['Relaxation ' num2str(jj)];
                end
            end
            hold off
        end
        
        if pltVr.unfDet == 1
            
            force{jj} = data.force{jj}+data.shiftY;
            ext{jj} = data.scombined{jj}+data.shiftX;
            l = length(force{jj});
            e1 = 0.78; %0.78 % threshold for force jumps defining start and end points [pN]
            e2 = 0.5; %0.5 % threshold for step selection, minimal force difference (between min and max in interval [start end]) [pN]
            n = 25; %25 % number of frames that correspond to half a second, upper limit for duration of unfolding event.
            a= [];
            % First find steps where the force decreases and then encircle the
            % unfolding event [startpoint endpoint]
            
            kk = 2;
            while kk < l-1
                if force{jj}(kk)-force{jj}(kk+1)> e1
                    endpoint = kk;
                    while force{jj}(endpoint)-force{jj}(endpoint+1)>0
                        endpoint = endpoint + 1;
                        if endpoint-kk >n || endpoint>l-1
                            display('unfolding event too long, or end of trace reached')
                            break
                        end
                    end
                    startpoint = kk;
                    while force{jj}(startpoint-1)-force{jj}(startpoint)>0
                        startpoint = startpoint - 1;
                        if kk - startpoint>n || startpoint <2
                            display('unfolding event too long, or end of trace reached')
                            break
                        end
                    end
                    a = [a; startpoint endpoint];
                    kk = endpoint+1;
                    
                else
                    kk = kk+1;
                end
            end
            %check if points were found
            if isempty(a)
                datVr.fePoints{jj} = [];
                disp('no steps found')
            end
            
            la = size(a);
            
            % now try to figure out which steps are actual steps
            % unfolding should consist of at least two frames, and the force step must
            % be big enough, and force before unfolding above 6 pN
            datVr.fePoints{jj} = [];
            
            for kk =1:la
                if abs(a(kk,2) - a(kk,1)) > 1 && abs(force{jj}(a(kk,1)) - force{jj}(a(kk,2))) > e2 && force{jj}(a(kk,2))>6
                    datVr.fePoints{jj} = [datVr.fePoints{jj}; a(kk,:)];
                end
            end
            
            if ~isempty(datVr.fePoints{jj})
                for kk = 1:size(datVr.fePoints{jj},1)
                    hold on;
                    data.figs.feUnf.axh{kk,jj} = plot([ext{jj}(datVr.fePoints{jj}(kk,1)) ext{jj}(datVr.fePoints{jj}(kk,2))], [force{jj}(datVr.fePoints{jj}(kk,1)) force{jj}(datVr.fePoints{jj}(kk,2))], 'linewidth', 2, 'color', 'k', 'linestyle', '-');
                    %(ext{jj}(points(ii,2))-ext{jj}(points(ii,1)))*1000
                end
            end
        end
        
        if pltVr.WLCfit == 1 && pltVr.unfDet == 1
            fMinOptions = optimset('MaxFunEvals', 600, 'TolFun', 1E-6, 'TolX', 1E-11); % fminsearch options.
            kB = datVr.dualCst.kB; %J/K
            T = datVr.dualCst.T; %K
            Ld = datVr.basepairs*datVr.baseL*1E-9;%920E-9; %nm
            Lp = datVr.residues*datVr.resL*1E-9;%120E-9;%370*0.%120E-9; %nm
            pP = datVr.pPROT; %nm
            pD = datVr.pDNA; %nm
            K = datVr.KDNA;
            P = pD;
            L = Ld;
            
            if size(datVr.fePoints{jj},1)>=2
                datVr.LpCalc = [];
                h = waitbar(0,'Initializing waitbar...');
                set(h,'Name','Finding Unfolding Events');
                tTot = 0;
                tCount = 0;
                for kk = 2:size(datVr.fePoints{jj},1)
                    
                    wlcForTemp = data.force{jj}(datVr.fePoints{jj}(kk-1,2):datVr.fePoints{jj}(kk,1))+data.shiftY;
                    wlcExtTemp = data.scombined{jj}(datVr.fePoints{jj}(kk-1,2):datVr.fePoints{jj}(kk,1))+data.shiftX;
                    xDNA = [];
                    for mm = 1:length(wlcForTemp)
                        xDNA(mm) = fminsearch(@(x) fitEWLC( wlcForTemp(mm)*1E-12, x(1), kB, T, P, L, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))
                    end
                    
                    %Calculate WLC ruler for getting protein contour length vs. time/force.
                    
                    for mm = 1:length(data.force{jj}(:))
                        tic
                        datVr.WLCRuler(mm,kk) = fminsearch(@(x) fitEWLC( (data.force{jj}(mm)+data.shiftY)*1E-12, x(1), kB, T, P, L, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))
                        datVr.contourL(mm,kk) = max(0,(data.scombined{jj}(mm)+data.shiftX) - datVr.WLCRuler(mm,kk)*1E6);
                        t = toc;
                        tCount = tCount + 1;
                        tTot = t + tTot;
                        hours = floor((tTot/tCount)*(length(data.force{jj}(:))-tCount)/3600);
                        mins = floor(((tTot/tCount)*(length(data.force{jj}(:))-tCount)-hours*3600)/60);
                        secs = (tTot/tCount)*(length(data.force{jj}(:))-tCount)-hours*3600-mins*60;
                        waitbar((mm)/(length(data.force{jj}(:))),h,[num2str(round(hours)) ':' num2str(round(mins)) ':' num2str(round(secs)) ' remaining.'])
                    end
                    
                    
                    xprot = (wlcExtTemp*1E-6-xDNA');
                    dnaCalcforce = calcEWLC( xDNA, kB, T, P, L, K  );
                    %plot((xDNA+xprot')*1E6,dnaCalcforce.*1E12,'linewidth',2);
                    K = datVr.KPROT;%1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
                    P = datVr.pPROT;   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
                    L = Lp;    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa
                    datVr.LpCalc(jj,kk-1) = fminsearch(@(x) calcProtLp( wlcForTemp*1E-12, xprot, kB, T, P, x(1)), 300E-9, fMinOptions);
                    xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(dnaCalcforce)-1):1.5*Lp*1E6].*(1E-6);
                    for mm = 1:length(dnaCalcforce)
                        xprot(mm) = fminsearch(@(x) fitEWLC( dnaCalcforce(mm), x(1), kB, T, P, datVr.LpCalc(jj,kk-1), K ), xprotGuess(mm), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
                    end
                    hold on;
                    plot((xDNA+xprot')*1E6,dnaCalcforce.*1E12,'linewidth',3)
                    hold off;
                    
                end
                close(h);
                
            end
        end
        
        
    end
    
    hold on
    
    plotRuler();
    
    hold off
    
    if pltVr.legOn == 1
        if pltVr.splitPulls ~= 1
            legend(data.figs.fext.axh(1:length(data.x)/2),pltVr.legArry)
        else
            legend(data.figs.fext.axh(1:length(data.x)),pltVr.legArry)
        end
        legend('show')
    else
        legend('hide')
    end
  
    if pltVr.rstAx == 1||isempty(pltVr.feAx)
        pltVr.rstAx = 0;
        pltVr.feAx = [0.5 1.5 -10  85];
        axis(pltVr.feAx)
        set(wndVr.sampling_xlim,'String','0.5 1.5')
        set(wndVr.sampling_ylim,'String','-10 85')
        set(wndVr.sampling_shift,'String',num2str([data.shiftX data.shiftY]));
    else
        axis(pltVr.feAx)
        set(wndVr.sampling_xlim,'String',num2str(get(gca,'xlim')));
        set(wndVr.sampling_ylim,'String',num2str(get(gca,'ylim')));
        set(wndVr.sampling_shift,'String',num2str([data.shiftX data.shiftY]));
    end
    
    
end
if pltVr.LpPlot == 1
    if isfield(pltVr.figs, 'contF')
        if ~isempty(findobj('name','Protein Contour Lengths'))
            close(pltVr.figs.contF)
        end
        rmfield(pltVr.figs, 'contF');
    end
    
    if isfield(pltVr.figs, 'contF')
        if iscell(pltVr.figs.contF)
            for kk = 1:length(pltVr.figs.contF)
                if ~isempty(findobj('name','Protein Contour Lengths'))
                    close(pltVr.figs.contF{kk})
                end
            end
        else
            if ~isempty(findobj('name','Protein Contour Lengths'))
                close(pltVr.figs.contF)
            end
        end
        rmfield(pltVr.figs, 'contF');
    end    
    pltVr.figs.contF = figure;
    if isfield(pltVr.figs, 'compF')
        if iscell(pltVr.figs.compF)
            for kk = 1:length(pltVr.figs.compF)
                if ~isempty(findobj('name','Protein Compaction vs. Force'))
                    close(pltVr.figs.compF{kk})
                end
            end
        else
            if ~isempty(findobj('name','Protein Compaction vs. Force'))
                close(pltVr.figs.compF)
            end
        end
        rmfield(pltVr.figs, 'compF');
    end
    pltVr.figs.compF = figure;
     if isfield(pltVr.figs, 'compFbar')
        if iscell(pltVr.figs.compFbar)
            for kk = 1:length(pltVr.figs.compFbar)
                if ~isempty(findobj('name','Protein Compaction vs. Force Mean'))
                    close(pltVr.figs.compFbar{kk})
                end
            end
        else
            if ~isempty(findobj('name','Protein Compaction vs. Force Mean'))
                close(pltVr.figs.compFbar)
            end
        end
        rmfield(pltVr.figs, 'compFbar');
    end
    pltVr.figs.compFbar = figure;
    h = waitbar(0,'Initializing waitbar...');
    set(h,'Name','Calculating Contour Lengths');
    tTot = 0;
    tCount = 0;
    datVr.contourL = [];
    datVr.forceSort = [];
    datVr.contourSort = [];
    data.intEnergyProt = [];
    data.intEnergyWLC = [];
    data.intEnergyTot = [];   
    fMinOptions = optimset('MaxFunEvals', 22600, 'TolFun', 1E-16, 'TolX', 1E-16); % fminsearch options.
    fMinOptions2 = optimset('MaxFunEvals', 22600, 'TolFun', 1E-16, 'TolX', 1E-16); % fminsearch options.
    kB = datVr.dualCst.kB; %J/K
    T = datVr.dualCst.T; %K
    Ld = datVr.basepairs*datVr.baseL*1E-9;%920E-9; %nm
    Lp = datVr.residues*datVr.resL*1E-9;%120E-9;%370*0.%120E-9; %nm
    pP = datVr.pPROT;%0.65E-9;%0.8E-9;%1.5E-9; %nm  % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
    pD = datVr.pDNA; %nm
    K = datVr.KDNA;
    KP = datVr.KPROT;%1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
    %Calculate WLC ruler for getting protein contour length vs. time/force.
    xDNA = [];
    set(0, 'currentfigure', f);
    FT = 0.01:0.01:100;
    FT = FT'.*1E-12;
    FTC = (FT.*pD./(kB.*T));
    xDNApl = Ld.*(4./3-4./(3.*sqrt(FTC+1))-10.*exp((900.*(1./FTC)).^(1./4))./(sqrt(FTC).*(exp((900.*(1./FTC)).^(1./4))-1).^2)+(FTC).^0.782843227494225./(8.728882712132439+1.948228355775715.*(FTC).^1.636511536819267) + (FTC./pD.*(kB.*T))./K );
    FTC = (FT.*pP./(kB.*T));
    xPROTpl = Lp.*(4./3-4./(3.*sqrt(FTC+1))-10.*exp((900.*(1./FTC)).^(1./4))./(sqrt(FTC).*(exp((900.*(1./FTC)).^(1./4))-1).^2)+(FTC).^0.782843227494225./(8.728882712132439+1.948228355775715.*(FTC).^1.636511536819267) + (FTC./pP.*(kB.*T))./KP );
    
    %keyboard;
    data.extFit = fit((xDNApl+xPROTpl).*1E6, FT.*1E12,  'smoothingspline', 'SmoothingParam', 0.9999999998982048);
    hold on;
    plot(xDNApl.*1E6,FT.*1E12,(xDNApl+xPROTpl).*1E6,FT.*1E12);
    hold off;
    
    for gg = 1:size(data.forceSplit,1)
        for oo = 1:size(data.forceSplit,2)
            if datVr.WLCacc == 1
                for mm = 1:length(data.forceSplit{gg,oo}(:))
                    tic
                    if data.forceSplit{gg,oo}(mm)+data.shiftY == 0
                        data.forceSplit{gg,oo}(mm) = data.forceSplit{gg,oo}(mm) + 0.001;
                    end
                    
                    xDNA{gg,oo}(mm,1) = fminsearch(@(x) fitEWLC( (data.forceSplit{gg,oo}(mm)+data.shiftY)*1E-12, x(1), kB, T, pD, Ld, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))
                    datVr.contourL{gg,oo}(mm,1) = fminbnd(@(x) fitEWLC( (data.forceSplit{gg,oo}(mm)+data.shiftY)*1E-12, (data.scombinedSplit{gg,oo}(mm)+data.shiftX)*1E-6-xDNA{gg,oo}(mm), kB, T, pP, x(1), KP ), 0E-9, Lp*1.3, fMinOptions2); %linspace(0, 0.9*Lp, length(extension))
                    
                    %datVr.contourL(mm,1) = max(0,(data.scombinedSplit{jj}(mm)+data.shiftX) - datVr.WLCRuler(mm,1)*1E6);
                    t = toc;
                    tCount = tCount + 1;
                    tTot = t + tTot;
                    hours = floor((tTot/tCount)*(length(data.forceSplit{gg,oo}(:))-tCount)/3600);
                    mins = floor(((tTot/tCount)*(length(data.forceSplit{gg,oo}(:))-tCount)-hours*3600)/60);
                    secs = (tTot/tCount)*(length(data.forceSplit{gg,oo}(:))-tCount)-hours*3600-mins*60;
                    waitbar((mm)/(length(data.forceSplit{gg,oo}(:))),h,[num2str(round(hours)) ':' num2str(round(mins)) ':' num2str(round(secs)) ' remaining.'])
                end
            else
                FT = (data.forceSplit{gg,oo}+data.shiftY).*1E-12;
                FTC = (FT.*pD./(kB.*T));
                xDNA{gg,oo} = Ld.*(4./3-4./(3.*sqrt(FTC+1))-10.*exp((900.*(1./FTC)).^(1./4))./(sqrt(FTC).*(exp((900.*(1./FTC)).^(1./4))-1).^2)+(FTC).^0.782843227494225./(8.728882712132439+1.948228355775715.*(FTC).^1.636511536819267) + (FTC./pD.*(kB.*T))./K );
                FT = (data.forceSplit{gg,oo}+data.shiftY).*1E-12;
                FTC = (FT.*pP./(kB.*T));
                datVr.contourL{gg,oo} = real(((data.scombinedSplit{gg,oo}+data.shiftX).*1E-6-xDNA{gg,oo})./(4./3-4./(3.*sqrt(FTC+1))-10.*exp((900.*(1./FTC)).^(1./4))./(sqrt(FTC).*(exp((900.*(1./FTC)).^(1./4))-1).^2)+(FTC).^0.782843227494225./(8.728882712132439+1.948228355775715.*(FTC).^1.636511536819267) + (FTC./pD.*(kB.*T))./KP )); %#ok<*NODEF>
            end
            for mm = 1:length(data.forceSplit{gg,oo}(:))
                if ((data.scombinedSplit{gg,oo}(mm)+data.shiftX)*1E-6-xDNA{gg,oo}(mm)) <= 0
                    datVr.contourL{gg,oo}(mm) = 0;
                end
            end
            
            [datVr.forceSort{gg,oo},I] = sort(data.forceSplit{gg,oo}+data.shiftY,1,'descend');
            datVr.contourSort{gg,oo} = datVr.contourL{gg,oo}(I);
            IIII = find(datVr.forceSort{gg,oo} > 2);
            datVr.forceSort{gg,oo} = datVr.forceSort{gg,oo}(IIII);
            datVr.contourSort{gg,oo} = datVr.contourSort{gg,oo}(IIII);
            if gg == 2
                set(0, 'currentfigure', pltVr.figs.compF);
                hold on;
                plot(decimate(datVr.contourL{gg,oo}(find(data.forceSplit{gg,oo}+data.shiftY > 2))*1E9./(datVr.residues*datVr.resL),1),decimate(data.forceSplit{gg,oo}(find(data.forceSplit{gg,oo}+data.shiftY > 2))+data.shiftY,1),'color',C(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:))
                hold off;
                set(0, 'currentfigure', pltVr.figs.compFbar);
                hold on;
                for tt = 3:2:floor(max(datVr.forceSort{gg,oo}))
                    errorbar(tt,mean(datVr.contourSort{gg,oo}(find(datVr.forceSort{gg,oo}<tt+1 & datVr.forceSort{gg,oo}>tt-1))).*1E9/(datVr.residues*datVr.resL),sqrt(var(datVr.contourSort{gg,oo}(find(datVr.forceSort{gg,oo}<tt+1 & datVr.forceSort{gg,oo}>tt-1)).*1E9/(datVr.residues*datVr.resL))),'-o','MarkerSize',5,'color',C(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:))
                    if tt == 3
                        datVr.FBar{1,oo}(1) = tt;
                        datVr.contBar{1,oo}(1) = mean(datVr.contourSort{gg,oo}(find(datVr.forceSort{gg,oo}<tt+1 & datVr.forceSort{gg,oo}>tt-1))).*1E9/(datVr.residues*datVr.resL);
                        datVr.errBar{1,oo}(1) = sqrt(var(datVr.contourSort{gg,oo}(find(datVr.forceSort{gg,oo}<tt+1 & datVr.forceSort{gg,oo}>tt-1)).*1E9/(datVr.residues*datVr.resL)));
                        
                    else
                        datVr.FBar{1,oo}(end+1) = tt;
                        datVr.contBar{1,oo}(end+1) = mean(datVr.contourSort{gg,oo}(find(datVr.forceSort{gg,oo}<tt+1 & datVr.forceSort{gg,oo}>tt-1))).*1E9/(datVr.residues*datVr.resL);
                        datVr.errBar{1,oo}(end+1) = sqrt(var(datVr.contourSort{gg,oo}(find(datVr.forceSort{gg,oo}<tt+1 & datVr.forceSort{gg,oo}>tt-1)).*1E9/(datVr.residues*datVr.resL)));
                        
                    end
                end
                
                hold off;
            end
            
            
             %%
            xPf = data.scombinedSplit{gg,oo}+data.shiftX;
            xShift = mean(xPf);
            [xPf,I] = sort(xPf);
            yPf = data.forceSplit{gg,oo}+data.shiftY;
            yPf = yPf(I);
            xPf = xPf-xShift;
            data.forceShift{gg,oo} = xShift;
            data.forceFit{gg,oo} = fit(xPf, yPf,  'smoothingspline', 'SmoothingParam', 0.9999999998982048);
            set(0, 'currentfigure', f);
            hold on;
            data.figs.fext.axh(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg)) = plot(sort(data.scombinedSplit{gg,oo})+data.shiftX, data.forceFit{gg,oo}(sort(data.scombinedSplit{gg,oo})+data.shiftX-data.forceShift{gg,oo}),'color',C(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:), 'XDataSource',strcat('data{',num2str(ii),'}.scombinedSplit'),'YDataSource', strcat('data{',num2str(ii),'}.forceSplit'));
      
            set(0, 'currentfigure', f);
            hold on;
            XX = sort(data.scombinedSplit{gg,oo})+data.shiftX;
            YY1 = data.extFit(sort(data.scombinedSplit{gg,oo})+data.shiftX);
            YY2 = data.forceFit{gg,oo}(sort(data.scombinedSplit{gg,oo})+data.shiftX-data.forceShift{gg,oo}); 
            
            %%% Draw integrated area and perform integrations
            fill([XX;flipud(XX)]',[YY1;flipud(max(YY1,YY2))]',C(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:))            
            data.intEnergyProt{gg,oo} = max(diff(integrate(data.forceFit{gg,oo},sort(data.scombinedSplit{gg,oo})+data.shiftX-data.forceShift{gg,oo},min(sort(data.scombinedSplit{gg,oo})+data.shiftX-data.forceShift{gg,oo}))),0);
            data.intEnergyWLC{gg,oo} = max(diff(integrate(data.extFit,sort(data.scombinedSplit{gg,oo})+data.shiftX,min(sort(data.scombinedSplit{gg,oo})+data.shiftX))),0);
            data.intEnergyTot{gg,oo} = sum(max(data.intEnergyWLC{gg,oo},data.intEnergyProt{gg,oo})-data.intEnergyWLC{gg,oo}).*1E-6.*1E-12./(kB.*T);          
            %%%
            if datVr.WLCacc == 1
                %%%% Fitting approximate WLC to accurate WLC
                FT = (data.forceSplit{gg,oo}+data.shiftY).*1E-12;
                FTC = (FT.*pD./(kB.*T));
                FT2 = (data.forceSplit{gg,oo}+data.shiftY).*1E-12;
                FTC2 = (FT.*pP./(kB.*T));
                XVar{gg,oo} = fminsearch(@(x) sum((xDNA{gg,oo} - (Ld.*(4./3-4./(3.*sqrt(FTC+1))-10.*exp((900.*(1./FTC)).^(1/4))./(sqrt(FTC).*(exp((900.*(1./FTC)).^(1./4))-1).^2)+(FTC).^x(1)./(x(2)+x(3).*(FTC).^x(4)) + FT./K ))).^2)+sum((datVr.contourL{gg,oo} - ((data.scombinedSplit{gg,oo}+data.shiftX).*1E-6-xDNA{gg,oo})./((4./3-4./(3.*sqrt(FTC2+1))-10.*exp((900.*(1./FTC2)).^(1/4))./(sqrt(FTC2).*(exp((900.*(1./FTC2)).^(1./4))-1).^2)+(FTC2).^x(1)./(x(2)+x(3).*(FTC2).^x(4)) + FT2./KP ))).^2), [1.62 3.55 3.8 2.2], fMinOptions);
                %%%%
            end
            set(0, 'currentfigure', pltVr.figs.contF);
            hold on;
            IIII = find(data.forceSplit{gg,oo}+data.shiftY > 2);
            if ~isempty(IIII)
                data.figs.contF.axh(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg)) = plot(datVr.contourL{gg,oo}(IIII).*1E9,data.forceSplit{gg,oo}(IIII)+data.shiftY,'color',C(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:));
            else
                data.figs.contF.axh(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg)) = plot(datVr.contourL{gg,oo}(:).*1E9,data.forceSplit{gg,oo}+data.shiftY,'color',C(oo*size(data.forceSplit,1)-size(data.forceSplit,1)+(gg),:));
            end
            clear IIII
            
            hold off;
        end
    end
    title(strrep(data.filename, '_', ' '));
    xlabel('Contour length (nm)')
    ylabel('Force (pN)');
    xlim([0 Lp.*1.3E9]);
    set(pltVr.figs.contF,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name','Protein Contour Lengths');
    if pltVr.legOn == 1
        if pltVr.splitPulls ~= 1
            legend(data.figs.contF.axh(1:length(data.x)/2),pltVr.legArry)
        else
            legend(data.figs.contF.axh(1:length(data.x)),pltVr.legArry)
        end
        legend('show')
    else
        legend('hide')
    end
    close(h)
    set(0, 'currentfigure', pltVr.figs.compF);
    set(pltVr.figs.compF,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name','Protein Compaction vs. Force');
    title(strrep(data.filename, '_', ' '));
    xlabel('Relative compaction (fraction uncompacted)')
    ylabel('Force (pN)');
    set(0, 'currentfigure', pltVr.figs.compFbar);
    hold on;
    plot(0:max([datVr.FBar{:,:}]),ones(1,length(0:max([datVr.FBar{:,:}]))),'k--');
    xlim([0 max([datVr.FBar{:,:}])+1]);
    hold off;
    set(pltVr.figs.compFbar,'BackingStore','off','MenuBar','figure','NumberTitle','off','Name','Protein Compaction vs. Force Mean');
    title(strrep(data.filename, '_', ' '));
    xlabel('Force (pN)');
    ylabel('Relative compaction (fraction uncompacted)')
    
end
pause(0.1);
end