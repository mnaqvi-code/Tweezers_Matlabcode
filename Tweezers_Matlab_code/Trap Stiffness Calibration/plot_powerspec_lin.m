function hf = plot_powerspec_lin(f,P,fNyq,color)

if nargin < 4, color = 'k'; end;
nbin = fNyq;
nbplot = floor(length(f)/nbin);
[count, centers] = hist(f,nbin);
delta   = centers(2) - centers(1);

clear fb_plot Pb_plot s_plot
ind = quickfind(f,nbin);
for i = 1 : nbin
        %ind         = find(f >= centers(i)-delta/2 & f < centers(i)+delta/2);
        fb_plot(i)  = mean_nan(f(ind(i,:)));
        Pb_plot(i)  = mean_nan(P(ind(i,:)));
        if sum(isfinite(P(ind(i,:)))) > 0, s_plot(i) = Pb_plot(i)/sqrt(sum(isfinite(P(ind(i,:))))); else s_plot(i) = NaN; end;
    end;

hold on;       

    %set(gca,'XLim',[10 max(f)]);
    %XL = get(gca,'XLim');
    %YL = get(gca,'YLim');
    %lerrb = ((log10(XL(2))-log10(XL(1)))/(2*nbin)); 
    %for i = 1 : length(s_plot),       %errorbars
    %    if s_plot(i) < Pb_plot(i),
    %        h = plot([fb_plot(i) fb_plot(i)],[Pb_plot(i)-s_plot(i) Pb_plot(i)+s_plot(i)],'k-'); 
    %            if ~isempty(color), set(h,'Color',color); end;
    %        h = plot([exp(log(fb_plot(i))-lerrb) exp(log(fb_plot(i))+lerrb)],[Pb_plot(i)-s_plot(i) Pb_plot(i)-s_plot(i)],'k-');
    %            if ~isempty(color), set(h,'Color',color); end;
    %        h = plot([exp(log(fb_plot(i))-lerrb) exp(log(fb_plot(i))+lerrb)],[Pb_plot(i)+s_plot(i) Pb_plot(i)+s_plot(i)],'k-');
    %            if ~isempty(color), set(h,'Color',color); end;
    %   end;
    %end;
    
hd = plot(fb_plot, Pb_plot, '-k'); set(hd,'LineWidth',0.7,'Color',color);

h = xlabel('Frequency (Hz)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
h = ylabel(['P(f) (arbitrary units) ']); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
set(gca,'XScale','log','YScale','log','FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
%STR(1) = {['Number of data points per block = ' num2str(nbplot,'%8.0f')]};
%caption(STR);

