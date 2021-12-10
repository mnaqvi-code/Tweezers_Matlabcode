function plot_data_div_fit(f,P,parameters,parameters0,color);

global fNyq nblock Ffit_start Ffit_end

if nargin < 4, color = 'k'; end;

%nbin = 70;
nbin = floor(length(f)/nblock);
[count, centers] = hist(f,nbin);
delta   = centers(2) - centers(1);

clear fb_plot Pb_plot nblock_plot
ind = quickfind(f,nbin);
    for i = 1 : nbin
        fb_plot(i)  = mean_nan(f(ind(i,:)));
        Pb_plot(i)  = mean_nan(P(ind(i,:)));
        nblock_plot(i) = sum(isfinite(P(ind(i,:))));
    end;

    fb_plot = fb_plot(:);
    Pb_plot = Pb_plot(:);
    nblock_plot = mean(nblock_plot);
    clear ind
    ind     = find(fb_plot > Ffit_start & fb_plot <= Ffit_end);
    fbp    = fb_plot(ind);
    Pbp    = Pb_plot(ind);
    %nbp    = nblock_plot(ind);


Pfit    =   P_theor(parameters,parameters0,fb_plot,0);
Pfit2   =   P_theor(parameters,parameters0,fbp,0);
E       =   Pb_plot ./ Pfit;
E2      =   Pbp ./ Pfit2;

hf = plot([min(fb_plot) max(fb_plot)], [1 1], 'k-'); set(hf,'LineWidth',2); hold on;

he1 = plot([min(fb_plot) max(fb_plot)], [1+1/sqrt(nblock_plot) 1+1/sqrt(nblock_plot)], 'k:'); set(he1,'LineWidth',1);
he2 = plot([min(fb_plot) max(fb_plot)], [1-1/sqrt(nblock_plot) 1-1/sqrt(nblock_plot)], 'k:'); set(he2,'LineWidth',1);

hd = plot(fb_plot, E, 'ko'); set(hd,'MarkerFaceColor','w','MarkerEdgeColor',color,'MarkerSize',3); 
hd2 = plot(fbp, E2, 'ko'); set(hd2,'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',3);
set(gca,'FontUnits','Normalized','FontSize',0.08); 

hx1 = xlabel('Frequency (Hz)','FontWeight','bold'); set(hx1,'FontUnits','Normalized','FontSize',0.08,'FontWeight','bold');
hx2 = ylabel('Data / Fit','FontWeight','bold'); set(hx2,'FontUnits','Normalized','FontSize',0.08,'FontWeight','bold');
