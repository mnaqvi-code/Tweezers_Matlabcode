function position_histogram_e(x,Nbin,Title)

global tolx     %Tolerance in fit

if nargin < 3, Title = ' '; end;
if nargin < 2, Nbin  = 50; end;

edges   =   [min(x) : (max(x)-min(x))/Nbin : max(x)];
N       =   histc(x,edges); 
Nred = N(1:end-1);                      %Remove last bin which is not really a bin (see help on histc)
delX    =   edges(2)-edges(1);          %Bin width
Xh      =   edges(1:end-1)' + delX/2; 
x_fit   =   min(edges)-2*delX : delX/100 : max(edges)+2*delX;
beta1=[mean(x) std(x) sum(Nred)];
beta1(3)=sum(Nred)*delX;
y_fit   =   free_gauss(beta1,x_fit); 
nfree = length(edges) - length(beta1);



axes;hold on;
h = bar(edges,N,'histc');
set(h,'facecolor','y')
h = plot(Xh,Nred,'k.'); set(h,'MarkerSize',6);
h = plot(x_fit,y_fit,'k-'); set(h,'LineWidth',2);
set(gca,'XLim',[min(x_fit) 1.5*max(x_fit)],'YLim',[0 max(N+abs(max(y_fit)-max(N)))],'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
hx = xlabel([Title 'position']); hy = ylabel('Count');
set(hx,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
set(hy,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');

curFigTemp = gcf;
title(['Distribution of ', Title, curFigTemp.Name(19), ' signal']);
clear curFigTemp;

axes('position',[.6 .55 .30 .35])
h = bar(edges,N,'histc');
hold on
set(h,'facecolor','y');
h = plot(Xh,Nred,'k.'); set(h,'MarkerSize',6);
h = plot(x_fit,y_fit,'k-'); set(h,'LineWidth',1);
set(gca,'XLim',[1.2*min(x_fit) 1.2*max(x_fit)],'YLim',[.5 max(N+abs(max(y_fit)-max(N)))],'FontUnits','normalized','FontSize',0.06,'FontWeight','bold');
hx = xlabel([Title 'position']); hy = ylabel('Count');
set(hx,'FontUnits','normalized','FontSize',0.06,'FontWeight','bold');
set(hy,'FontUnits','normalized','FontSize',0.06,'FontWeight','bold');
set(gca,'Yscale','log');

STR(1) = {[char(10) char(10) char(10) 'Mean = ' num2str(beta1(1),'%5.3f')  ]};
STR(2) = {['Std = ' num2str(beta1(2),'%5.3f')   ]};
caption(STR,0.1,0.15);






