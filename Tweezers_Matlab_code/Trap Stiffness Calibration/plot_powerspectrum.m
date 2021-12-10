function hd = plot_powerspectrum(f,P,color,turn)

if nargin < 3, color = 'k'; end;
persistent indp         %Reuse indp next time this function is called

%f = f(f>=0.5 & f<=1132.5);
%P = P(f>=0.5 & f<=1132.5);
logf=log(f);
nbin = 70; %ceil(length(f)/5000);
%keyboard
%max(P(f>=31.5 & f<=32.5))
[count, centers] = hist(logf,nbin);
delta   = centers(2) - centers(1);
clear fb_plot Pb_plot s_plot
if  (turn == 1)
    for i = 1 : nbin
        indp(i).indx  = find(logf >= centers(i)-delta/2 & logf < centers(i)+delta/2);
    end
end

for i = 1:nbin    
    fb_plot(i)  = mean_nan(f(indp(i).indx));
    Pb_plot(i)  = mean_nan(P(indp(i).indx));
    if sum(isfinite(P(indp(i).indx))) > 0, s_plot(i) = Pb_plot(i)/sqrt(sum(isfinite(P(indp(i).indx)))); else s_plot(i) = NaN; end;
end;
hold on;       
set(gca,'XLim',[10 max(f)]);
XL = get(gca,'XLim');
YL = get(gca,'YLim');
lerrb = ((log10(XL(2))-log10(XL(1)))/(2*nbin)); 
for i = 1 : length(s_plot)       %errorbars
    if s_plot(i) < Pb_plot(i)
        h = plot([fb_plot(i) fb_plot(i)],[Pb_plot(i)-s_plot(i) Pb_plot(i)+s_plot(i)],'k-'); 
        if ~isempty(color), set(h,'Color',color); end;
        h = plot([exp(log(fb_plot(i))-lerrb) exp(log(fb_plot(i))+lerrb)],[Pb_plot(i)-s_plot(i) Pb_plot(i)-s_plot(i)],'k-');
        if ~isempty(color), set(h,'Color',color); end;
        h = plot([exp(log(fb_plot(i))-lerrb) exp(log(fb_plot(i))+lerrb)],[Pb_plot(i)+s_plot(i) Pb_plot(i)+s_plot(i)],'k-');
        if ~isempty(color), set(h,'Color',color); end;
    end;
end;

hd = plot(fb_plot, Pb_plot, 'k.'); set(hd,'MarkerSize',10,'MarkerFaceColor',color,'MarkerEdgeColor',color);

h = xlabel('Frequency (Hz)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
h = ylabel(['P(f) (arbitrary units) ']); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
set(gca,'XScale','log','YScale','log','FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
