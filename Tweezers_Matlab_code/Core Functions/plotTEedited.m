function [ data ] = plotTE( varargin ) % draws force-extension graphs for all chosen files

% if nargin == 1 %number of arguments
%     index = 1:size(data,2);
% end

global data f pltVr datVr; %#ok<REDEF>

pltVr.grphAct = 2;
arrayfun(@cla,f);
%cla(findall(0,'type','axes'))

if isfield(data,'figs')
    if isfield(data.figs,'te')
        data.figs = rmfield(data.figs,'te');
    end
end

C = linspecer(size(1:2:length(data.x),2));
if pltVr.legOn == 1
    pltVr.legArry = cell(length(data.x)/2,1);
end
%data{i}.figs.fext.figh = figure('Name',strcat(data{i}.filename,' F-Ext'),...
%    'NumberTitle','off','Renderer','zbuffer'); % draws figure frame with the name "filenumber.dat F-Ext"
%data.figs.fta.axh = plot(data.time,data.fx,'XDataSource',strcat('data{',num2str(1),'}.time'),...
%    'YDataSource', strcat('data{',num2str(1),'}.fx'));
 

%datacursormode on
datVr.tePoints = [];
data.figs.teUnf.axh = [];


for jj =  1:length(data.x)/2       
    data.figs.te.axh(jj) = plot(data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), data.scombined{jj}+data.shiftX,'color',C(jj,:), 'XDataSource',strcat('data{',num2str(1),'}.time'),'YDataSource', strcat('data{',num2str(1),'}.scombined'));
    hold on;
    pltVr.legArry{jj,1} = ['Cycle ' num2str(jj)];
    
    if pltVr.unfDet == 1;
    
        time{jj} = data.time(pltVr.lIdx(jj):pltVr.rIdx(jj));
        ext{jj} = data.scombined{jj}+data.shiftX;
        l = length(ext{jj});
        e1 = 0.01; % threshold for extension jumps defining start and end points [um]
        e2 = 0.00010; % threshold for step selection, minimal extension difference (between min and max in interval [start end]) [um]
        n = 10; % number of frames that correspond to half a second, upper limit for duration of unfolding event.
        a= [];
        % First find steps where the force decreases and then encircle the
        % unfolding event [startpoint endpoint]

        kk = 2;
        while kk <l-1
            if abs(ext{jj}(kk+1)-ext{jj}(kk))> e1
               startpoint = kk;
               endpoint = kk;
               while abs(ext{jj}(endpoint+1)-ext{jj}(endpoint))>0 && abs(ext{jj}(endpoint+1)-ext{jj}(startpoint)) > abs(ext{jj}(endpoint)-ext{jj}(startpoint))
                   endpoint = endpoint + 1;
                   if endpoint-kk >n || endpoint>l-1
                       endpoint = endpoint - 1;
                       display('unfolding event too long, or end of trace reached')
                       break
                   end
               end
               
%                while ext{jj}(startpoint-1)-ext{jj}(startpoint)>0
%                    startpoint = startpoint - 1;
%                    if kk - startpoint>n || startpoint <2
%                        startpoint = startpoint + 1;
%                         display('unfolding event too long, or end of trace reached')
%                        break
%                    end
%                end
               a = [a; startpoint endpoint];
               kk = endpoint+1;

            else
               kk = kk+1;
            end
        end
        %check if points were found
        if isempty(a)
            datVr.tePoints{jj} = [];
            disp('no steps found')        
        end   

        la = size(a);

        % now try to figure out which steps are actual steps
        % unfolding should consist of at least two frames, and the force step must
        % be big enough, and force before unfolding above 6 pN
        datVr.tePoints{jj} = [];    

        for kk =1:la
            if abs(a(kk,2) - a(kk,1)) >= 1 && abs(ext{jj}(a(kk,1)) - ext{jj}(a(kk,2))) > e2 %&& ext{jj}(a(kk,2))>0.5
                datVr.tePoints{jj} = [datVr.tePoints{jj}; a(kk,:)];
            end
        end

        if ~isempty(datVr.tePoints{jj})
            for ii = 1:size(datVr.tePoints{jj},1)
                hold on;
                data.figs.teUnf.axh{ii,jj} = plot([time{jj}(datVr.tePoints{jj}(ii,1)) time{jj}(datVr.tePoints{jj}(ii,2))], [ext{jj}(datVr.tePoints{jj}(ii,1)) ext{jj}(datVr.tePoints{jj}(ii,2))], 'linewidth', 2, 'color', 'k');
                %(ext{jj}(points(ii,2))-ext{jj}(points(ii,1)))*1000
            end
        end
    end

    
    
end
hold off;
title(data.filename);
xlabel('Time (s)')
ylabel('Extension (\mum)');
set(gca,'Tag',num2str(1))

if pltVr.legOn == 1        
    legend(data.figs.te.axh(1:length(data.x)/2),pltVr.legArry)
    legend('show')
else
    legend('hide')
end


hold off;
if pltVr.rstAx == 1||isempty(pltVr.teAx)
    pltVr.rstAx = 0;
    pltVr.teAx = [get(gca,'xlim') get(gca,'ylim')];
    axis(pltVr.teAx) 
else
    axis(pltVr.teAx)        
end   
end

