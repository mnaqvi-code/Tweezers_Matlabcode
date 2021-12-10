function [ data ] = plotSDual( varargin )

global data pltVr datVr prgVr; %#ok<REDEF>

pltVr.grphAct = 4;
% if nargin == 4 %number of arguments
%     index = 1:size(data,2);
% end
arrayfun(@cla,f);
%cla(findall(0,'type','axes'))
for ii = 1:size(data,2);
    if isfield(data,'figs')
        if isfield(data.figs,'fs')
            data.figs = rmfield(data.figs,'fs');
        end
    end
end

for ii = 1:size(data,2);
    if ~isfield(data, 'x')
        data.x(1) = min(data.time);
        data.x(2) = max(data.time);
    end
    %if isfield(data, 'x')
    
    C = linspecer(size(1:2:length(data.x),2));
    for jj =  1:2:length(data.x)
        
        
        [pltVr.lIdx,pltVr.lIdx] = min(abs(data.x(jj)-data.time));
        [pltVr.rIdx,pltVr.rIdx] = min(abs(data.x(jj+1)-data.time));
        
        if prgVr.forStt == 1 %Fsmooth
            
            %xPos = smooth((data.xspt1*0.0615) - mean(data.xspt2*0.0615), 50);
            [X,I] = sort(data.vx(pltVr.lIdx:pltVr.rIdx));
            Y = data.xspt2(pltVr.lIdx:pltVr.rIdx)*0.0615;
            Y = Y(I);
            
            if length(pltVr.fitVal) == 2
                J = intersect(find(Y >= (min(Y) + pltVr.fitVal(1)*(max(Y) - min(Y)))), find(Y <= (min(Y) + pltVr.fitVal(2)*(max(Y) - min(Y)))));
                pp = splinefit(X(J),Y(J),5,4); %pp = csaps(X(J), Y(J), datVr.fitPar); %f1 = fit(X(J), Y(J),  'smoothingspline', 'SmoothingParam', datVr.fitPar);
            else
                pp = splinefit(X,Y,5,4); %pp = csaps(X, Y, datVr.fitPar); %f1 = fit(X, Y,  'smoothingspline', 'SmoothingParam', datVr.fitPar);
            end        
        elseif prgVr.forStt == 2 %Flinear
            
            [X,I] = sort(data.vx(pltVr.lIdx:pltVr.rIdx)); %sort(data.xspt2(pltVr.lIdx:pltVr.rIdx)*0.0615);
            Y = data.xspt2(pltVr.lIdx:pltVr.rIdx)*0.0615;
            Y = Y(I);            
            if length(pltVr.fitVal) == 2
                J = intersect(find(Y >= (min(Y) + pltVr.fitVal(1)*(max(Y) - min(Y)))), find(Y <= (min(Y) + pltVr.fitVal(2)*(max(Y) - min(Y)))));                
                Pfit = polyfit(Y(J),X(J),1);
            else
                Pfit = polyfit(Y,X,1);
            end
            [X,I] = sort(data.vx);
          %  X = data.vx(I);
            Y = data.xspt2*0.0615;
            Y = Y(I);
            Pdev{jj} =  polyval(polyder(Pfit),Y);
        end
        
        
        
        %         Pfit = polyfit(X(I),Y(I),fitOrder);
        %         [X,I] = sort(data.xspt2*0.0615);
        %         Y = data.vx(I);
        %         Pdev{ii} = polyval(polyder(Pfit),X);%3*Pfit(1)*X.^2+2*Pfit(2).*X +Pfit(3);
        %         [d1,d2] = differentiate(f1,X);
        %         Pdev{ii} = d1;
        %
        %         unsorted = 1:size(data.xspt2,1);
        %         newInd(I) = unsorted;
        %         Pdev{ii} = Pdev{ii}(unsorted);
        % Voltage sensitivities. [default 2.74?]
        %dvdx(ii) = data.dvdx./0.4;
        %-data.vx*100./Pdev{ii}
        %data.figs.fext.axh = plot(-data.vx./dvdx(ii)-xPos,data.vx*160/dvdx(ii),'XDataSource',strcat('data{',num2str(ii),'}.ext'),'YDataSource', strcat('data{',num2str(ii),'}.fx'));
        %data.figs.fext.axh =
        %plot(-data.vx./dvdx(ii)-xPos,data.fx); %'XDataSource',strcat('data{',num2str(ii),'}.ext'),'YDataSource', strcat('data{',num2str(ii),'}.fx')
        data.figs.fs.axh((jj+1)/2) = plot(data.xspt2(pltVr.lIdx:pltVr.rIdx).*0.0615,data.vx(pltVr.lIdx:pltVr.rIdx),'XDataSource',strcat('data{',num2str(ii),'}.xspt2'),'YDataSource', strcat('data{',num2str(ii),'}.vx'),'color',C((jj+1)/2,:));

        hold on;
        if prgVr.forStt == 1
            data.figs.fsfit.axh((jj+1)/2) = plot(ppval(pp,X),X,'color','k','linewidth',1);
        elseif prgVr.forStt == 2   
            data.figs.fsfit.axh((jj+1)/2) = plot(sort(data.xspt2(pltVr.lIdx:pltVr.rIdx).*0.0615),polyval(Pfit,sort(data.xspt2(pltVr.lIdx:pltVr.rIdx).*0.0615)),'color', 'k');
        end
        
        %-(data.xpz+data.vx/data.dvdx);
        %plot(data.xspt2(pltVr.lIdx:pltVr.rIdx).*0.0615,polyval(Pfit,data.xspt2(pltVr.lIdx:pltVr.rIdx).*0.0615),'color', 'k')
        title(strrep(data.filename, '_', ' '))
        xlabel('Xspt2 (µm)')
        ylabel('Voltage (V)');
        %datacursormode on  
        set(gca,'Tag',num2str(ii))
        hold on
        %figure
        %plot(data.yspt2(pltVr.lIdx:pltVr.rIdx).*0.0615,data.vy(pltVr.lIdx:pltVr.rIdx))
        %plot(sqrt((data.xspt2(pltVr.lIdx:pltVr.rIdx)-data.xspt2(1)).^2+(data.yspt2(pltVr.lIdx:pltVr.rIdx)-data.yspt2(1)).^2).*0.0615,sqrt((data.vx(pltVr.lIdx:pltVr.rIdx)-data.vx(1)).^2+(data.vy(pltVr.lIdx:pltVr.rIdx)-data.vy(1)).^2))
    end
    
    hold off
    %     else
    %         %data{i}.figs.fext.figh = figure('Name',strcat(data{i}.filename,' F-Ext'),...
    %         %    'NumberTitle','off','Renderer','zbuffer'); % draws figure frame with the name "filenumber.dat F-Ext"
    %         data.figs.fext.axh = plot(data.xspt2.*0.0615,data.vx,'XDataSource',strcat('data{',num2str(ii),'}.xspt2'),'YDataSource', strcat('data{',num2str(ii),'}.vx'));
    %         title(data.filename);
    %         [X,I] = sort(data.xspt2*0.0615);
    %         Y = data.vx(I);
    %         Pfit = polyfit(X,Y,2);
    %         disp(Pfit)
    %         hold on;
    %         plot(X,polyval(Pfit,X),'color', 'r')
    %         xlabel('Xspt2 (µm)')
    %         ylabel('Voltage (V)');
    %         %datacursormode on
    %         set(gca,'Tag',num2str(ii))
    %         hold off
    %     end
end
    if pltVr.rstAx == 1||isempty(pltVr.ssAx)
        pltVr.rstAx = 0;
        pltVr.ssAx = [get(gca,'xlim') get(gca,'ylim')];
        axis(pltVr.ssAx) 
    else
        axis(pltVr.ssAx)        
    end   

end

