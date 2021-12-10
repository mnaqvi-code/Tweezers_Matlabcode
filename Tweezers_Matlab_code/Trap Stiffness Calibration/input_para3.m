% INPUT_PARA3 reads WANT_HYDRO (choice of including hydrodynamic corrections), 
% R (radius of the bead), L (distance to surface), NU (kinematic
% viscosity), and RHO (density of the bead), and a choice of characteristic filter functions. 
% Finally, FIT_POWERSPECTRUM button is shown. 

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
want_hydro  =   0;              % want_hydro = 1 for including hydrodynamic corrections; otherwise 0 
R           =   0.5e-6;         % radius of the bead (m)
l           =   5e-6;           % distance to surface (m)
nu          =   1e-12/1e-6;     % kinematic viscosity (m2/s)
rho         =   1e-3/1e-6;      % density of the bead (kg/m3)

ht  = uicontrol('Style','text', ...
    'Position',[left3-edge bot5-edge wid1+wid2+3*edge bot_title+height-(bot5-edge)+0.02], ...
    'String',' ',...
    'BackgroundColor',color);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
ht  = uicontrol('Style','text', ...
    'FontSize',titlefs, ...
    'Position',[left3 bot_title wid1+wid2+edge height], ...
    'String','Frequency-dependent hydrodynamic friction',...
    'FontSize',titlefs+0.2,...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hhydro  = uicontrol('Style','popup', ...
    'Position',[left3 bot1 wid1+wid2+edge height], ...
    'String',' Do not include frequency-dependence of friction | Include frequency-dependence of friction', ...
    'Callback','input_para3_hydro; hydro_decision;', ...
    'BackgroundColor',color1);
input_para3_hydro;

hydro_decision;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Filters

ht  = uicontrol('Style','text', ...
    'Position',[left3-edge 0.02 wid1+wid2+3*edge bot5-0.02-2*edge], ...
    'String',' ',...
    'BackgroundColor',color);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
ht  = uicontrol('Style','text', ...
    'FontSize',titlefs, ...
    'Position',[left3 bot7 wid1+wid2+edge height], ...
    'String','Filter corrections',...
    'FontSize',titlefs+0.2,...
    'BackgroundColor',color1);
jh = findjobj(ht);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hfilter  = uicontrol('Style','popup', ...
    'Position',[left3 bot8 wid1+wid2+edge height], ...
    'String',' Type your own characteristic function of the filter | 1 first-order filter | 2 first-order filters | 3 first-order filters | 4 first-order filters | 5 first-order filters | No filter', ...
    'Callback','filter_decision', ...
    'BackgroundColor',color1);
filter_decision;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Fitting

hfit_powerspectrum  = uicontrol('Style','Pushbutton', ...
    'Position',[left3 bot11 1.2*wid1 height], ...
    'String','Fit power spectrum/spectra',...
    'FontSize',titlefs+0.2,...
    'BackgroundColor',[1 0.1 0.1],...
    'Callback','fit_powerspectrum');
jh = findjobj(hfit_powerspectrum);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
hphist = uicontrol('Style','Pushbutton', ...
    'Position',[left3 bot12 1.2*wid1 height], ...
    'String','Show distribution of Pexp/Pfit',... 
    'FontSize',titlefs+0.1,...
    'BackgroundColor',[1 0.1 0.1],...
    'Callback','spectrum_histogram');
jh = findjobj(hphist);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );
close_windows = uicontrol('Style','Pushbutton', ...
    'Position',[left3+(wid1+wid2) 0 wid1/2 height], ...
    'String','Close all',...
    'BackgroundColor',[.5 .5 .5],...
    'Callback','ccc');
jh = findjobj(close_windows);
jh.setVerticalAlignment( javax.swing.AbstractButton.CENTER );

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx