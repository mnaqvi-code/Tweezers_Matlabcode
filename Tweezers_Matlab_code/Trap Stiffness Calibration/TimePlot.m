function [ exportplot ] = TimePlot( data, stepSize, k, alpha, exportplot )

if sum(ishandle(exportplot)) > 0
    delete(exportplot(:));
    clear exportplot;
end

hold on;
extColor = linspecer(length(data));

for ii = 1:length(data)

exportplot(ii) = plot(data{ii}.data(1:stepSize:end,1),abs(-data{ii}.data(1:stepSize:end,8)-min(-data{ii}.data(1:stepSize:end,8))).*k.*alpha,'Color',extColor(ii,:)); title(data{ii}.textdata(3,8));

end
hold off;