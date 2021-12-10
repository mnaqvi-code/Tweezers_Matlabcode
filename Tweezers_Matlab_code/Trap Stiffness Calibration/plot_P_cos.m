function plot_P_cos(f,P,parameters,parameters0,color);

global fNyq nblock Plot_start Plot_end Ffit_start Ffit_end

if nargin < 4, color = 'k'; end;

%nbin = 70;
nbin = floor(length(f)/nblock);
[count, centers] = hist(f,nbin);
delta   = centers(2) - centers(1);

clear fb_plot Pb_plot nblock_plot

ind = quickfind(f,nbin);
for i = 1 : nbin
    %ind         = find(f >= centers(i)-delta/2 & f < centers(i)+delta/2);
    fb_plot(i)  = mean_nan(f(ind(i,:)));
    Pb_plot(i)  = mean_nan(P(ind(i,:)));
    nblock_plot(i) = sum(isfinite(P(ind(i,:))));
end;

fb_plot = fb_plot(:);
Pb_plot = Pb_plot(:);
nblock_plot = nblock_plot(:);
clear ind
ind     = find(fb_plot > Ffit_start & fb_plot <= Ffit_end);
fbp    = fb_plot(ind);
Pbp    = Pb_plot(ind);
%nbp    = nblock_plot(ind);

Pfit = P_theor(parameters,parameters0,fb_plot,0);
hf = plot(cos(pi*fb_plot/fNyq), 1./Pfit,'r-'); set(hf,'LineWidth',2,'Color',color); hold on;       

he1 = plot(cos(pi*fb_plot/fNyq), 1./Pfit+(1./Pfit)./sqrt(nblock_plot), 'k:'); set(he1,'LineWidth',0.5,'Color',color);
he2 = plot(cos(pi*fb_plot/fNyq), 1./Pfit-(1./Pfit)./sqrt(nblock_plot), 'k:'); set(he2,'LineWidth',0.5,'Color',color);

hd = plot(cos(pi*fb_plot/fNyq), 1./Pb_plot, 'ko'); set(hd,'MarkerFaceColor','w','MarkerEdgeColor',color,'MarkerSize',2); 
hd2 = plot(cos(pi*fbp/fNyq), 1./Pbp, 'ko'); set(hd2,'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',2); 
set(gca,'FontUnits','Normalized','FontSize',0.08); 

h = xlabel('cos(pi*f/f_{Nyq})'); set(h,'FontUnits','normalized','FontSize',0.08,'FontWeight','bold');
h = ylabel(['1/P(f) (arbitrary units)']); set(h,'FontUnits','normalized','FontSize',0.08,'FontWeight','bold');
