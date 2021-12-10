function plot_corrxy(f,Px,Py,Pxy,color,s)

% Px,Py, and Pxy are already binned
% s indicate which symbol

global fNyq

if nargin < 5
    color = [0 0 0];
end;

fun = Pxy ./ sqrt(Px.*Py);

nblock  =   round(length(f)/10); % to get about 10 final data points
nbin    =   floor(length(f)/nblock);
funb = []; sb = [];
for i = 1 : nbin
    fb(i)   = mean_nan(f((i-1)*nblock+1 : i*nblock));
    funb(i) = mean_nan(fun((i-1)*nblock+1 : i*nblock));
    sb(i)   = std_nan(fun((i-1)*nblock+1 : i*nblock));
end;
fz(1) = 0.;
fz(2) = fNyq/1e3;
yz(1) = 0.;
yz(2) = 0.;

%h = plot(fb/1e3,funb,'k-'); set(h,'Color',color);
if s == 1
    h = errorbar(fb/1e3,funb,sb,'ko'); %set(h,'Color',color,'MarkerFaceColor','w','MarkerEdgeColor',color,'MarkerSize',6);
    h = plot(fz,yz,'k-');
else 
    % plot these points slightly displaced;
    h = errorbar((fb+0.01*fNyq)/1e3,funb,sb,'ks'); set(h,'Color',color,'MarkerFaceColor','w','MarkerEdgeColor',color,'MarkerSize',6);
end; % if s == 1
h = xlabel('Frequency (kHz)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
h = ylabel('P_{xy} / sqrt(P_xP_y)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
set(gca,'XLim',[0 fNyq/1e3],'FontWeight','bold');

