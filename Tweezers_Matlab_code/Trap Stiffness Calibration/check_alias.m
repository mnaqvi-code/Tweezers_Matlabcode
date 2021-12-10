%this mfile creates a popup menu which gives the user the possibility to 
%increase the number of terms to include in the anti aliasing sum
alias_tol = 0.01/sqrt(nblock);
[P, alias_check] = P_theor(scal_fit,parameters0,xfin,1);
if (alias_check > alias_tol)
    set(0,'Units','pixels') 
    scnsize = get(0,'ScreenSize');
    h17 = figure(17);
    set(h17,'MenuBar','none','Color','k','Name','Message','NumberTitle','off','Position',[scnsize(3)/10 scnsize(4)/10 2*scnsize(3)/3 2*scnsize(4)/4])
    uicontrol('style','frame','Position',[scnsize(3)/100+10 scnsize(4)/100+10 2*scnsize(3)/3-40 2*scnsize(4)/4-40],'BackgroundColor',[.7 .7 .7]);
    uicontrol('style','text','Position',[scnsize(3)/10 scnsize(4)/4 2*scnsize(3)/4 2*scnsize(4)/10],...
    'Fontsize',11,'horizontalalignment','left','string','The precision of the fit you just performed can be improved if you increase the number of terms n in the aliasing sum in the mfile P_theor.m. The problem with your fit is probably that the corner frequency is too close to the Nyquist frequency, f_Nyq = f_sample/2. Try increasing your sampling frequency in a new experiment or increase the number of terms in the aliasing sum and press the OK button below'...
    ,'BackgroundColor',color1,'FontWeight','bold');
    alias_corr = uicontrol('style','popupmenu','Position',[scnsize(3)/6 1 scnsize(3)/3 scnsize(4)/5],...
    'string','Increase the number of aliasing terms to 20|Increase the number of aliasing terms to 30|Increase the number of aliasing terms to 40|Satisfied with the fit n = 10'...
    ,'BackgroundColor',color1,'FontWeight','bold','Callback','change_alias');
    new_fit = uicontrol('style','pushbutton','Position',[scnsize(3)/5+scnsize(3)/30 scnsize(4)/15 scnsize(3)/5 scnsize(4)/15],...
    'string','OK','FontWeight','bold'...
    ,'BackgroundColor',[1 0 0],'Callback','fit_powerspectrum;close(h17)');
end