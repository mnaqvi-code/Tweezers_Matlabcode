function [ data ] = plotFECTrap( varargin ) % draws force-extension graphs for all chosen files
%plotFE - plots Force vs. Extension data for selected intervals.
%   plotFE( data, p_CCD, p_QPD_piezo, rulVal, datVr.sigVal, forVal, prgVr.forStt, fitVal, fitPar, index )

global p_CCD p_QPD_piezo f data pltVr datVr prgVr wndVr; %#ok<REDEF>
%fitVal = fitVal;

% if nargin == 1 %number of arguments
%     index = 1:size(data,2);
% end
pltVr.grphAct = 3;

edge = 0.01;
arrayfun(@cla,f);
%cla(findall(0,'type','axes'))
for ii = 1:size(data,2)
    if isfield(data,'figs')
        if isfield(data.figs,'fext')
            data.figs = rmfield(data.figs,'fext');
        end
    end
end
if ~isempty(f)
    set(0,'CurrentFigure',1)
end
                

if strcmp(data.Data.Root.Property(2).Value{1},'5')
    data.time = data.Data.MeasuredData(7).Data*1E-3;
elseif strcmp(data.Data.Root.Property(2).Value{1},'6')
    data.time = data.Data.MeasuredData(13).Data*1E-3;
end
data.time = decimate(data.time,datVr.decFact);

for ii = 1:size(data,2)    
    if ~isfield(data, 'x')
        data.x(1) = min(data.time);
        data.x(2) = max(data.time);
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
%             if ~exist('ppOffset')
%                 ppOffset = max(Pdev{jj});                
%             else
%                 ppOffset = max(ppOffset,max(Pdev{jj}));                
%             end
%             Pdev{jj} = Pdev{jj} -  ppOffset;
            
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
            pltVr.lIdx(jj) = pltVr.lIdx(jj)*datVr.decFact-datVr.decFact+1;
            pltVr.rIdx(jj) = pltVr.rIdx(jj)*datVr.decFact-datVr.decFact+1; 
            %data.force{jj} = data.Data.MeasuredData(3).Data(pltVr.lIdx(jj):pltVr.rIdx(jj));
            if strcmp(data.Data.Root.Property(2).Value{1},'5')
                data.force{jj} = (data.Data.MeasuredData(3).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-data.Data.MeasuredData(5).Data(pltVr.lIdx(jj):pltVr.rIdx(jj)))./2;
                data.scombined{jj} = data.Data.MeasuredData(8).Data(pltVr.lIdx(jj):pltVr.rIdx(jj));
            elseif strcmp(data.Data.Root.Property(2).Value{1},'6')
                data.force{jj} = sqrt(((data.Data.MeasuredData(5).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-data.Data.MeasuredData(7).Data(pltVr.lIdx(jj):pltVr.rIdx(jj)))./2).^2 + ((data.Data.MeasuredData(8).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-data.Data.MeasuredData(6).Data(pltVr.lIdx(jj):pltVr.rIdx(jj)))./2).^2);
                data.scombined{jj} = data.Data.MeasuredData(3).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))-data.Data.MeasuredData(5).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(5).Property(4).Value.*1E6+data.Data.MeasuredData(7).Data(pltVr.lIdx(jj):pltVr.rIdx(jj))./data.Data.MeasuredData(7).Property(4).Value.*1E6;
            end
            pltVr.lIdx(jj) = (pltVr.lIdx(jj)+datVr.decFact-1)/datVr.decFact;
            pltVr.rIdx(jj) = (pltVr.rIdx(jj)+datVr.decFact-1)/datVr.decFact;
            
            data.force{jj} = decimate(data.force{jj},datVr.decFact);
            data.scombined{jj} = decimate(data.scombined{jj},datVr.decFact);
            
            
            %Calculate Bias Correction
            
            if ~isempty(datVr.xBias) && ~isempty(datVr.yBias)
                biasP = polyfit(datVr.xBias,datVr.yBias,1);
                biasFit = biasP(1).*data.scombined{jj};
                data.force{jj} = data.force{jj}-biasFit;
            end
            
            if pltVr.splitPulls == 1
                [valtemp,Itemp] = max(data.scombined{jj});
                data.scombinedSplit{1,jj} = data.scombined{jj}(1:Itemp);
                data.scombinedSplit{2,jj} = data.scombined{jj}(Itemp+1:end);
                data.forceSplit{1,jj} = data.force{jj}(1:Itemp);
                data.forceSplit{2,jj} = data.force{jj}(Itemp+1:end);
            end   
            

        end
 
    end
    if(wndVr.hTrig == 0)
        wndVr.hTrig = 1;
        maxTemp = max(cellfun(@(x)max(x(:)), data.scombined));
        minTemp = min(cellfun(@(x)min(x(:)), data.scombined));
        set(wndVr.hSlider, 'Min', min(1*(0-maxTemp), 1*(0-minTemp)), 'Max', max(1*(0-minTemp)+0.9, 1*(0-maxTemp)+0.9), 'Value', mean([min(1*(0-maxTemp)+0.9, 1*(0-minTemp)+0.9) max(1*(0-minTemp)+0.9, 1*(0-maxTemp)+0.9)]), 'SliderStep', abs([1 1]./((1*(0-maxTemp))-(1*(0-minTemp))*100)));
        data.shiftX = get(wndVr.hSlider,'Value');
        data.shiftY = get(wndVr.vSlider,'Value');
        set(wndVr.sampling_shift,'String', num2str([data.shiftX data.shiftY]));
    end;

    if(wndVr.vTrig == 0)
        wndVr.vTrig = 1;
        maxTemp = max(cellfun(@(x)max(x(:)), data.force));
        minTemp = min(cellfun(@(x)min(x(:)), data.force));
        set(wndVr.vSlider, 'Min', min(1*(0-maxTemp)-0.9, 1*(0-minTemp)-0.9), 'Max', max(1*(0-minTemp)+0.9, 1*(0-maxTemp)+0.9), 'Value', max([min(1*(0-maxTemp), 1*(0-minTemp)) max(1*(0-minTemp), 1*(0-maxTemp))]), 'SliderStep', abs([1 1]./((1*(0-maxTemp))-(1*(0-minTemp))*100)));     
        data.shiftX = get(wndVr.hSlider,'Value');
        data.shiftY = get(wndVr.vSlider,'Value');
        set(wndVr.sampling_shift,'String', num2str([data.shiftX data.shiftY]));
    end;
    
    datVr.fePoints = [];
    data.figs.feUnf.axh = [];
    
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
                data.figs.fext.axh(jj*2-1+(gg-1)) = plot(data.scombinedSplit{gg,jj}+data.shiftX, data.forceSplit{gg,jj}+data.shiftY,'color',C(jj*2-1+(gg-1),:), 'XDataSource',strcat('data{',num2str(ii),'}.scombinedSplit'),'YDataSource', strcat('data{',num2str(ii),'}.forceSplit'));
                %data.figs.fext.axh(jj) = plot(((data.dataInt{jj}(:,6)-data.dataInt{jj}(:,7))-min(data.dataInt{jj}(:,6)-data.dataInt{jj}(:,7))).*0.1+0.6+data.shiftX,data.dataInt{jj}(:,2)+data.shiftY,'color',C(jj,:), 'XDataSource',strcat('data{',num2str(ii),'}.scombined'),'YDataSource', strcat('data{',num2str(ii),'}.force'));
                title(data.Data.Root.Property(1).Value{1});
                xlabel('Extension (µm)')
                ylabel('Force (pN)');
                set(gca,'Tag',num2str(ii))        
                
                graphArray = cell(length(data.x)/2,1);
                if gg == 1
                    pltVr.legArry{jj*2-1+(gg-1),1} = ['Pull ' num2str(jj)];
                else
                    pltVr.legArry{jj*2-1+(gg-1),1} = ['Relaxation ' num2str(jj)];
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
            fMinOptions = optimset('MaxFunEvals', 2000, 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
            kB = 1.38064852*10^(-23); %J/K
            T = 293.15; %K
            Ld = datVr.basepairs*datVr.baseL*1E-9;%920E-9; %nm
            Lp = datVr.residues*datVr.resL*1E-9;%120E-9;%370*0.%120E-9; %nm
            pP = 1.5E-9; %nm
            pD = 45E-9; %nm
            K = 1.2E-9;
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
                    xprot = xprot;
                    dnaCalcforce = calcEWLC( xDNA, kB, T, P, L, K  );
                    %plot((xDNA+xprot')*1E6,dnaCalcforce.*1E12,'linewidth',2);
                    K = 1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
                    P = 1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
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
    hold on
    
    plotRuler();
    
    hold off
    if pltVr.rstAx == 1||isempty(pltVr.feAx)   
        pltVr.rstAx = 0;
        pltVr.feAx = [0.5 1.5 -10  85];
        axis(pltVr.feAx) 
        set(wndVr.sampling_xlim,'String','0.5 1.5')
        set(wndVr.sampling_ylim,'String','-10 85')
    else
        axis(pltVr.feAx)        
    end   
    
   
end
    if pltVr.LpPlot == 1
        h = waitbar(0,'Initializing waitbar...');
        set(h,'Name','Calculating Contour Lengths');
        tTot = 0;
        tCount = 0;
        datVr.contourL = [];
        fMinOptions = optimset('MaxFunEvals', 2000, 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
        fMinOptions2 = optimset('MaxFunEvals', 2000, 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
        kB = 1.38064852*10^(-23); %J/K
        T = 293.15; %K
        Ld = datVr.basepairs*datVr.baseL*1E-9;%920E-9; %nm
        Lp = datVr.residues*datVr.resL*1E-9;%120E-9;%370*0.%120E-9; %nm
        pP = 0.65E-9;%0.8E-9;%1.5E-9; %nm  % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
        pD = 45E-9; %nm
        K = 1.2E-9;
        KP = 2.0E-9;%1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
        L = Lp;           
        %Calculate WLC ruler for getting protein contour length vs. time/force.
        
        xDNA = [];    
        for mm = 1:length(data.force{1}(:))
            tic
            xDNA(mm,1) = fminsearch(@(x) fitEWLC( (data.force{1}(mm)+data.shiftY)*1E-12, x(1), kB, T, pD, Ld, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))

            datVr.contourL(mm,1) = fminbnd(@(x) fitEWLC( (data.force{1}(mm)+data.shiftY)*1E-12, (data.scombined{1}(mm)+data.shiftX)*1E-6-xDNA(mm,1), kB, T, pP, x(1), KP ), 0E-9, Lp.*1.3, fMinOptions2); %linspace(0, 0.9*Lp, length(extension))
            
            if ((data.scombined{1}(mm)+data.shiftX)*1E-6-xDNA(mm,1)) <= 0
                datVr.contourL(mm,1) = 0;
            end
            
            %datVr.contourL(mm,1) = max(0,(data.scombined{jj}(mm)+data.shiftX) - datVr.WLCRuler(mm,1)*1E6);
            t = toc;
            tCount = tCount + 1;
            tTot = t + tTot;
            hours = floor((tTot/tCount)*(length(data.force{1}(:))-tCount)/3600);
            mins = floor(((tTot/tCount)*(length(data.force{1}(:))-tCount)-hours*3600)/60);
            secs = (tTot/tCount)*(length(data.force{1}(:))-tCount)-hours*3600-mins*60;
            waitbar((mm)/(length(data.force{1}(:))),h,[num2str(round(hours)) ':' num2str(round(mins)) ':' num2str(round(secs)) ' remaining.'])
        end
        figure;
        hold on;
        IIII = find(data.force{1}+data.shiftY > 6);
        if ~isempty(IIII)
            plot(datVr.contourL(IIII,1).*1E9,data.force{1}(IIII,1)+data.shiftY)
        else
           plot(datVr.contourL(:,1).*1E9,data.force{1}+data.shiftY) 
        end
        clear IIII
        title(strrep(data.filename, '_', ' '));
        xlabel('Contour length (nm)')
        ylabel('Force (pN)');
        xlim([0 130]);
        hold off;
        close(h)
    end

pause(0.1);
end