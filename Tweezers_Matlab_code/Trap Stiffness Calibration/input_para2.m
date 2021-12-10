% INPUT_PARA2 reads NBLOCK (number of data points in a block), TOLX
% (termination tolerance on X in the fit; see file LSQNONLIN in the Optimization 
% Toolobox in /toolbox/optim and OPTIMSET in /toolbox/matlab/funfun),
% LFIT_START and LFIT_END (range of the Lorentzian fit), FFIT_START and FFIT_END 
% (range of the final fit), PLOT_START and PLOT_END (plotting range),
% and WANT_ALPHA (choice of fitting the diode with 1 or 2 parameters).

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if prgVr.setupStt == 0 %% Old Tweezers
    % Default values:
    nblock      =   150;   % number of data points in a block
    tolx        =   1e-7;   % tolerance in the fit
    Lfit_start  =   100;    % range of the Lorentzian fit (Hz)
    Lfit_end    =   20000;   %  -"-
    Ffit_start  =   100;    % range of the final fit (Hz)
    Ffit_end    =   20000;   %  -"-
    Plot_start  =   100;     % plotting range (Hz)
    Plot_end    =   20000;   %  -"-
    want_alias  =   1;      % want_alias = 1 for accounting for aliasing;
    want_diode  =   1;      % want_diode = 0 for no filtering by diode;
    want_alpha  =   0;      % want_alpha = 0 for filtering by diode with 1 parameter (f3dBeff); want_alpha = 1 for fitting filtering by diode with 2 parameters (alpha and f3dB); 

elseif prgVr.setupStt == 1 %% Foldometer  
        % Default values:
    nblock      =   150;   % number of data points in a block
    tolx        =   1e-7;   % tolerance in the fit
    Lfit_start  =   100;    % range of the Lorentzian fit (Hz)
    Lfit_end    =   20000;   %  -"-
    Ffit_start  =   100;    % range of the final fit (Hz)
    Ffit_end    =   20000;   %  -"-
    Plot_start  =   100;     % plotting range (Hz)
    Plot_end    =   20000;   %  -"-
    want_alias  =   1;      % want_alias = 1 for accounting for aliasing;
    want_diode  =   1;      % want_diode = 0 for no filtering by diode;
    want_alpha  =   0;      % want_alpha = 0 for filtering by diode with 1 parameter (f3dBeff); want_alpha = 1 for fitting filtering by diode with 2 parameters (alpha and f3dB); 

elseif prgVr.setupStt == 2 %% CTRAP
    % Default values:
    nblock      =   150;   % number of data points in a block
    tolx        =   1e-7;   % tolerance in the fit
    Lfit_start  =   1000;    % range of the Lorentzian fit (Hz)
    Lfit_end    =   23000;   %  -"-
    Ffit_start  =   1000;    % range of the final fit (Hz)
    Ffit_end    =   23000;   %  -"-
    Plot_start  =   2;     % plotting range (Hz)
    Plot_end    =   25000;   %  -"-
    want_alias  =   2;      % want_alias = 1 for accounting for aliasing;
    want_diode  =   0;      % want_diode = 0 for no filtering by diode;
    want_alpha  =   0;      % want_alpha = 0 for filtering by diode with 1 parameter (f3dBeff); want_alpha = 1 for fitting filtering by diode with 2 parameters (alpha and f3dB); 

end


ht  = uicontrol('Style','text', ...
    'Position',[left1 bbot7-0.5*bbot13 wid1 height], ...
    'String','Number of data points in a block',...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hnblock = uicontrol('Style','edit', ...
    'Position',[left2 bbot7-0.5*bbot13 wid2 height], ...
    'String',num2str(nblock),...
    'BackgroundColor',color2,...
    'Callback','nblock = str2num(get(hnblock,''String''));');

ht  = uicontrol('Style','text', ...
    'Position',[left1 bbot8-0.5*bbot13 wid1 height], ...
    'String','Termination tolerance in fit',...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
htolx = uicontrol('Style','edit', ...
    'Position',[left2 bbot8-0.5*bbot13 wid2 height], ...
    'String',num2str(tolx,'%5.1e'),...
    'BackgroundColor',color2,...
    'Callback','tolx = str2num(get(htolx,''String''));');

w = 0.5*wid2; left2_2 = left2+0.6*wid2;
ht  = uicontrol('Style','text', ...
    'Position',[left1 bbot9-0.5*bbot13 wid1 height], ...
    'String','Fitting range for the Lorentzian fit (Hz)',...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hLfit_start = uicontrol('Style','edit', ...
    'Position',[left2 bbot9-0.5*bbot13 w height], ...
    'String',num2str(Lfit_start),...
    'BackgroundColor',color2,...
    'Callback','Lfit_start = str2num(get(hLfit_start,''String''));');
hLfit_end = uicontrol('Style','edit', ...
    'Position',[left2_2 bbot9-0.5*bbot13 w height], ...
    'String',num2str(Lfit_end),...
    'BackgroundColor',color2,...
    'Callback','Lfit_end = str2num(get(hLfit_end,''String''));');

ht  = uicontrol('Style','text', ...
    'Position',[left1 bbot10-0.5*bbot13 wid1 height], ...
    'String','Fitting range for the final fit (Hz)',...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hFfit_start = uicontrol('Style','edit', ...
    'Position',[left2 bbot10-0.5*bbot13 w height], ...
    'String',num2str(Ffit_start),...
    'BackgroundColor',color2,...
    'Callback','Ffit_start = str2num(get(hFfit_start,''String''));');
hFfit_end = uicontrol('Style','edit', ...
    'Position',[left2_2 bbot10-0.5*bbot13 w height], ...
    'String',num2str(Ffit_end),...
    'BackgroundColor',color2,...
    'Callback','Ffit_end = str2num(get(hFfit_end,''String''));');

ht  = uicontrol('Style','text', ...
    'Position',[left1 bbot11-0.5*bbot13 wid1 height], ...
    'String','Plotting range (Hz)',...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hPlot_start = uicontrol('Style','edit', ...
    'Position',[left2 bbot11-0.5*bbot13 w height], ...
    'String',num2str(Plot_start),...
    'BackgroundColor',color2,...
    'Callback','Plot_start = str2num(get(hPlot_start,''String''));');
hPlot_end = uicontrol('Style','edit', ...
    'Position',[left2_2 bbot11-0.5*bbot13 w height], ...
    'String',num2str(Plot_end),...
    'BackgroundColor',color2,...
    'Callback','Plot_end = str2num(get(hPlot_end,''String''));');

halias  = uicontrol('Style','popup', ...
    'Position',[left1 0.8*bbot12-0.4*bbot13 wid2+left2-left1 height], ...
    'String',...
    ' Account for aliasing | Do not account for aliasing ', ...
    'Callback','alias_decision', ...
    'BackgroundColor',color1);
set(halias,'Value',want_alias);

hdiode  = uicontrol('Style','popup', ...
    'Position',[left1 0.6*bot12+0.3*bot13 wid2+left2-left1 height], ...
    'String',...
    ' Position detector not treated as virtual low-pass filter | Detector treated as filter with 1 parameter (fdiode-eff) | Detector treated as filter with 2 parameters (alpha and fdiode)', ...
    'Callback','diode_decision',...
    'BackgroundColor',color1);
set(hdiode,'Value',want_diode+want_alpha+1);