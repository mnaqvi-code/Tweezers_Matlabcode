function plot_Phist(f,P,scal_fit,parameters,color);
global fNyq Ffit_start Ffit_end

if nargin < 4, color = 'k'; end;

Nbinp = 100;
   
    ind     = find(f > Ffit_start & f <= Ffit_end);
    f_plot    = f(ind);
    P_plot    = P(ind);
    Pfit = zeros(size(f_plot));
    di  =  10;
    nex =  floor(length(f_plot)/di);
    dfex = (Ffit_end - Ffit_start)/(nex-1);
    f_ex    = Ffit_start : dfex : Ffit_end + 2*dfex;
    P_ex    = P_theor(scal_fit,parameters,f_ex,0);
    f_plot    = f(ind);
    P_plot    = P(ind);
    k = 1;
    for i = 1:di:length(f_plot)-di
        for j = 0:di-1
            Pfit(i + j)  =  P_ex(k) + (P_ex(k+1)-P_ex(k))*(f_plot(i + j)-f_ex(k))/dfex;
        end %j
        k = k + 1;
    end %i
    
ind0 = find(Pfit ~= 0);
Pp = P_plot(ind0);
fp = f_plot(ind0);
Pf = Pfit(ind0);

r = Pp ./ Pf;
edges   =   [min(r) : (max(r)-min(r))/Nbinp : max(r)]; 
N       =   histc(r,edges);

delX    =   edges(2)-edges(1);
Ntot    =   sum(N);
N       =   N./(Ntot*delX);
Xh      =   edges(1:length(edges))' + delX/2; 
x_fit   =   min(edges)-2*delX : delX/100 : max(edges)+2*delX;
y_fit   =   exp(-x_fit);
y_err   =   exp(-x_fit/2)./sqrt(Ntot*delX);

hh = plot(Xh,N,'k.'); set(hh,'MarkerSize',6); hold on;
h1 = plot(x_fit,y_fit,'k-'); set(h1,'LineWidth',1);
h2 = plot(x_fit,y_fit + y_err, 'k:'); set(h2,'LineWidth',0.5,'Color',color);
index = find((y_fit - y_err) > 0);
h3 = plot(x_fit(index),y_fit(index) - y_err(index), 'k:'); set(h3,'LineWidth',0.5,'Color',color);
hx = xlabel('P^{(exp)}/P^{(fit)}'); set(hx,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
hy = ylabel(['Distribution']); set(hy,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
set(gca,'YScale','log','XLim',[0 max(x_fit)],'YLim',[0.1/(Ntot*delX) ceil(max(y_fit+sqrt(y_fit)))]);

drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);