function [ data ] = plotFT( varargin ) % draws force-extension graphs for all chosen files

% if nargin == 1 %number of arguments
%     index = 1:size(data,2);
% end

global data f pltVr prgVr; %#ok<REDEF>

pltVr.grphAct = 1;

arrayfun(@cla,f)%arrayfun(@cla,findall(0,'type','axes'))
%cla(findall(0,'type','axes'))

if isfield(data,'figs')
    if isfield(data.figs,'ft')
        data.figs = rmfield(data.figs,'ft');
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
set(0, 'currentfigure', f);
for jj =  1:length(data.x)/2  
    data.figs.ft.axh(jj) = plot(data.time(pltVr.lIdx(jj):pltVr.rIdx(jj)), data.force{jj}+data.shiftY,'color',C(jj,:), 'XDataSource',strcat('data{',num2str(1),'}.time'),'YDataSource', strcat('data{',num2str(1),'}.force'));
    hold on;
    pltVr.legArry{jj,1} = ['Cycle ' num2str(jj)];
end
if prgVr.setupStt == 2
    title(strrep(data.filename, '_', ' '));
else
    title(strrep(data.filename, '_', ' '));
    %title(data.Data.Root.Property(1).Value{1})
end

xlabel('Time (s)')
ylabel('Force (pN)');
set(gca,'Tag',num2str(1))

if pltVr.legOn == 1        
    legend(data.figs.ft.axh(1:length(data.x)/2),pltVr.legArry)
    legend('show')
else
    legend('hide')
end


hold off;
if pltVr.rstAx == 1||isempty(pltVr.ftAx)
    pltVr.rstAx = 0;
    pltVr.ftAx = [get(gca,'xlim') get(gca,'ylim')];
    axis(pltVr.ftAx) 
else
    axis(pltVr.ftAx)        
end   
end

