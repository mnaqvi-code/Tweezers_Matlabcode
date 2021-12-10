function hf = plot_fit(f,P,parameters,parameters0,color,turn);

global nblock Plot_start Plot_end Ffit_start Ffit_end
persistent indp

if nargin < 4, color = 'k'; end;
logf=log(f);
%nbin = 70;
nbin = floor(length(f)/nblock);
[count, centers] = hist(logf,nbin);
delta   = centers(2) - centers(1);

if (turn == 1)
    for i = 1 : nbin,
        indp(i).indx         = find(logf >= centers(i)-delta/2 & logf < centers(i)+delta/2);
    end
end

clear fb_plot Pb_plot nblock_plot
for i = 1 : nbin
    fb_plot(i)  = mean_nan(f(indp(i).indx));
    Pb_plot(i)  = mean_nan(P(indp(i).indx));
    nblock_plot(i) = sum(isfinite(P(indp(i).indx)));
end;

index = find(nblock_plot~=0);
fb_plot = fb_plot(:); fb_plot = fb_plot(index);
Pb_plot = Pb_plot(:); Pb_plot = Pb_plot(index);
nblock_plot = nblock_plot(:); nblock_plot = nblock_plot(index);


ind    = find(fb_plot > Ffit_start & fb_plot <= Ffit_end);
fbp    = fb_plot(ind);
Pbp    = Pb_plot(ind);
%nbp    = nblock_plot(ind);


Pfit = P_theor(parameters,parameters0,fb_plot,0);

hf = plot(fb_plot, Pfit, 'r-'); set(hf,'LineWidth',2,'Color',color); hold on;       

he1 = plot(fb_plot, Pfit + Pfit./sqrt(nblock_plot), 'k:'); set(he1,'LineWidth',0.5,'Color',color);
he2 = plot(fb_plot, Pfit - Pfit./sqrt(nblock_plot), 'k:'); set(he2,'LineWidth',0.5,'Color',color);

hd = plot(fb_plot, Pb_plot, 'ko'); set(hd,'MarkerFaceColor','w','MarkerEdgeColor',color,'MarkerSize',3); 
hd2 = plot(fbp, Pbp, 'ko'); set(hd2,'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',3);

h = xlabel('Frequency (Hz)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
h = ylabel(['P(f) (arbitrary units) ']); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
set(gca,'XScale','log','YScale','log','FontUnits','normalized','FontSize',0.04,'FontWeight','bold',...
    'XLim',[Plot_start Plot_end]);