function caption(STR,pos_left,pos_up)

if nargin == 1
    pos_left = 0.1;
    pos_up = 0;
end;

hcap = axes('Position',[0 0 1 1],'Visible','off');
set(gcf,'CurrentAxes',hcap);
text(pos_left+0.02,0.1+pos_up,STR,'FontUnits','normalized','FontSize',0.021,'fontweight','bold');
