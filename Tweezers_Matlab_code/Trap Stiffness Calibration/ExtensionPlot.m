function [ exportplot ] = ExtensionPlot( data, stepSize, k, alpha, startp, endp, exportplot )

if sum(ishandle(exportplot)) > 0
    delete(exportplot(:));
    clear exportplot;
end

hold on;
extColor = linspecer(length(data));

for ii = 1:length(data)

%xDat = -data{ii}.data(endp:stepSize:startp,10).*0.000012369726100557; %0.00001815556572892595; 
xDat = -data{ii}.data(end:stepSize:1,10).*0.000012369726100557; %0.00001815556572892595; 
[Y,I] = sort(xDat);
xValU = unique(Y);
%yDat = abs(-data{ii}.data(endp:stepSize:startp,8)-min(-data{ii}.data(endp:stepSize:startp,8))).*k.*alpha;    
    
yDat = abs(-data{ii}.data(end:stepSize:1,8)-min(-data{ii}.data(end:stepSize:1,8))).*k.*alpha; 

exportplot(ii) = plot(xDat,yDat,'Color',extColor(ii,:)); title(data{ii}.textdata(3,8));

end
hold off;