%% TWEEZERSTART - Initializes the tweezers data analysis program.
% Variables used in the program:  
% fitPar: The value controls the fit.
% rulVal: This value controls the ruler used.
% sigVal : This value controls which signal is used for the extension: 1 is
% piezo signal, 2 is the camera tracking signal and 3 is the combined
% signal as described in the paper "Noise reduction by signal combination
% in Fourier space applied to drift correction in optical tweezers."%
% forVal: This value contains the trap stiffness value [pN/um].%
% forState: This value is 1 for a spline sensitivity fit and 2 for a linear sensitivity fit. %
% fitVal: This is a 1x1 or 2x1 vector containing a fit value or fit
% interval.
%% INITIALIZATION Even though not programmatically nice, a lot of variables are defined as globals in the program, since not much interaction is done outside the program. Functions will call the global variables that they use at the start of their script.
MSGID = 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame';
warning('off', MSGID)
if exist('cPos')
    clear cPos;
    close all;
else
    clear all;
    close all;
    clc;
    if(~isdeployed)
        tmp = matlab.desktop.editor.getActive;
        cd(fileparts(tmp.Filename));
    end
    addpath(genpath(pwd))  
    global binaryD calData;
    calData = [];
    binaryD = 0; % Toggle for using binary format or not. 
end
PrgInit = 1;
global p_CCD p_QPD_piezo data fluodata pltVr prgVr datVr datWndVr datWndLabVr f g fluof;
data = [];
fluodata = [];


prgVr.forStt = 2; %1 = Fsmooth, 2 = Flinear
prgVr.setupStt = 1; %Toggle between Older Tweezers Setup, Foldometer Setup and Ctrap.
%% Data Variables

datWndLabVr.wlcCst = [];
datVr.ver = 'v1.16';
datVr.fitPar = 0.9; %Piece-wise Curve Fitting Paramater
datVr.sigVal = 1; % This value controls which signal is used for the extension: 1 is
datVr.forVal = 0; %This value contains the trap stiffness value [pN/um].
datVr.tePoints = []; %This value holds the unfolding event regions in the Time-Extension Trace.
datVr.fePoints = []; %This value holds the unfolding event regions in the Force-Extension Trace.
datVr.WLCacc = 0; %This value sets if calculated WLC fits and contour lengths are calculated analyitically or approximated (improving calculation times drastically).
datVr.LpCalc = []; %This value holds any fitted WLC contour lengths for proteins.
datVr.CustRul = []; %This value holds a custom loaded ruler.
datVr.WLCRuler = []; %This value holds the WLC Force-extension ruler to calculate protein contour length.
datVr.contourL = []; %This value holds the protein contour length.
datVr.forceSort = []; %This value holds sorted forces, in a descending order.
datVr.contourSort = []; %This value holds contour lengths sorted after descending force.
datVr.FBar = [];
datVr.contBar = [];
datVr.errBar = [];
datWndVr.contourL = []; 
datWndLabVr.contourL = []; 
datVr.basepairs = 2666;%558*2; %This value holds the custom amount of basepairs for the DNA eWLC length.
datWndVr.basepairs = [];
datWndLabVr.basepairs = [];
datVr.baseL = 0.34; %This value holds the size of a basepair in nm.
datWndVr.baseL = [];
datWndLabVr.baseL = [];
datVr.pDNA = 45; %This value holds the persistence length of the dsDNA (nm).
datWndVr.pDNA = [];
datWndLabVr.pDNA = [];
datVr.KDNA = 1.2E-9; % This value holds the stiffness of the dsDNA.
datWndVr.KDNA = [];
datWndLabVr.KDNA = [];
datVr.residues = 396;%184;%365; %This value holds the custom amount of resideus for the protein eWLC length.
datWndVr.residues = [];
datWndLabVr.residues = [];
datVr.resL = 134.64/396;%0.36; %0.33 %This value holds the length of a residue in nm.
datWndVr.resL = [];
datWndLabVr.resL = [];
datVr.pPROT = 1.5E-9; %This value holds the persistence length of the protein (nm).
datWndVr.pPROT = [];
datWndLabVr.pPROT = [];
datVr.KPROT = 2.0E-9; %This value holds the stiffness of the protein.
datWndVr.KPROT = [];
datWndLabVr.KPROT = [];

datWndLabVr.rulB = [];
datVr.xBias = []; %Will contain X-coordinates to calculate bias correction.
datVr.yBias = []; %Will contain Y-coordinates to calculate bias correction.
datVr.decFact = 1; %Data decimation factor.
datVr.Nsplit = 8; %Calibration data number of windows to average over.
%% Dual trapping constants
datWndLabVr.physCst = [];
datVr.dualCst.kB = 1.3806488*10^(-23); %J/K
datWndVr.dualCst.kB = [];
datWndLabVr.dualCst.kB = [];
datVr.dualCst.T = 294.15; %K (26 C )
datWndVr.dualCst.T = [];
datWndLabVr.dualCst.T = [];
datVr.dualCst.TC = datVr.dualCst.T-273.15; %C
datWndVr.dualCst.TC = [];
datWndLabVr.dualCst.TC = [];
datVr.dualCst.m = 250E-3; %M NaCl
datWndVr.dualCst.m = [];
datWndLabVr.dualCst.m = [];
datVr.dualCst.m2 = 0; %M +2 Ions
datWndVr.dualCst.m2 = [];
datWndLabVr.dualCst.m2 = [];
datVr.dualCst.viscC = [0.1256735; 1.265347; -1.105369; 0.2044679; 1.308779];
datWndVr.dualCst.viscC = [];
datWndLabVr.dualCst.viscC = [];
datVr.dualCst.alphaC = [-0.04296718; 0.3710073; 0.4230889; -0.3259828];
datWndVr.dualCst.alphaC = [];
datWndLabVr.dualCst.alphaC = [];
datVr.dualCst.eta = (datVr.dualCst.viscC(1) +  datVr.dualCst.viscC(2)*exp(datVr.dualCst.alphaC(1).*datVr.dualCst.TC) + datVr.dualCst.viscC(3)*exp(datVr.dualCst.alphaC(2).*datVr.dualCst.m) + datVr.dualCst.viscC(4)*exp(datVr.dualCst.alphaC(3).*(0.01.*datVr.dualCst.TC + datVr.dualCst.m)) + datVr.dualCst.viscC(5)*exp(datVr.dualCst.alphaC(4).*(0.01.*datVr.dualCst.TC - datVr.dualCst.m)))*1E-3;%0.000891937; %Pa*s (measured for 37 C w rheometer)
datWndVr.dualCst.eta = [];
datWndLabVr.dualCst.eta = [];

%Fit Parameters Trap 1 - X
datWndLabVr.trap1Cst = [];
datVr.dualCst.Dex1 =  538049.7;%3388577.1;%718007.4;%858644.2;%3388577.1;%858644.2;%3388577.1;%3265273.4;%3400077.3;%3593096.1;%3164569.0;%2096250.5;%3164569.0;%2988709.9;%2049168.7;%3960212.5;%3176996.9;%3510159.2;%2851943.5;%1899175.3;%4017662;%1981958.8;%3102422.8;%5994570.1;%5588493.9;%3605859.8;%5588493.9;%4616429.6;%4842977.2;%3258786.9; %3380695.9;%3441669.6;%3706348.3; %4272991.0; %4392557.4; %4057743.8; %3990864.4; 4133423.7; %4712070.3; %a.u.^2/s
datWndVr.dualCst.Dex1 = [];
datVr.dualCst.cFreq1 =  3320.3;%3191.5;%3117.4;%3191.5;%3981.5;%3191.5;%3238.2;%3401.7;%3198.6;%3556.6;%3084.5;%3247.9;%3362.7;%4609.6;%3197.3;% 3932.3;%3208;%3179.7;%3257.3;%2824.4;%4683.7;%4386.4;%4307.3;%4386.4;%4385.6;%4108.2; %3245;%3301.5;%4214.4; %3520.1; %3582.9; %2050; %2139.5; 2135.1; 2148.9; %2237.6; %2305; %2182.8; %Hz
datWndVr.dualCst.cFreq1 = [];
datVr.dualCst.R1 = 0.0000011; %m
datWndVr.dualCst.R1 = [];
datVr.dualCst.gam1 = 6*pi*datVr.dualCst.eta*datVr.dualCst.R1; %N*s/m
datWndVr.dualCst.gam1 = [];
datVr.dualCst.k1 = 2*pi*datVr.dualCst.gam1*datVr.dualCst.cFreq1*10^6; %pN/um
datWndVr.dualCst.k1 = [];
datVr.dualCst.Dth1 = datVr.dualCst.kB*datVr.dualCst.T/datVr.dualCst.gam1; %m^2/s
datWndVr.dualCst.Dth1 = [];
datVr.dualCst.alpha1 = sqrt(datVr.dualCst.Dth1/datVr.dualCst.Dex1)*10^6; %um/a.u. %0.00047568
datWndVr.dualCst.alpha1 = []; 

%Fit Parameters Trap 1 - Y
datWndLabVr.trap3Cst = [];
datVr.dualCst.Dex3 =  251373.4;
datWndVr.dualCst.Dex3 =  [];
datVr.dualCst.cFreq3 = 2376.0; %Hz
datWndVr.dualCst.cFreq3 = [];
datVr.dualCst.R3 = 0.0000011; %m
datWndVr.dualCst.R3 = [];
datVr.dualCst.gam3 = 6*pi*datVr.dualCst.eta*datVr.dualCst.R3; %N*s/m
datWndVr.dualCst.gam3 = [];
datVr.dualCst.k3 = 2*pi*datVr.dualCst.gam3*datVr.dualCst.cFreq3*10^6; %pN/um
datWndVr.dualCst.k3 = [];
datVr.dualCst.Dth3 =  datVr.dualCst.kB*datVr.dualCst.T/datVr.dualCst.gam3; %m^2/s
datWndVr.dualCst.Dth3 = [];
datVr.dualCst.alpha3 = sqrt(datVr.dualCst.Dth3/datVr.dualCst.Dex3)*10^6; %um/a.u.
datWndVr.dualCst.alpha3 = [];

%Fit Parameters Trap 2 - X
datWndLabVr.trap2Cst = [];
datVr.dualCst.Dex2  =  2012759.2;%1939800.8;%2140639.6;%1939800.8;%2052545.5;%1939800.8;%1722659.0;%2052545.5;%2107779.4;%2323934.5;%2796250.5;%1306373.3;%2531644.7;%1601046.7;%2422523.2;%5440914.1;%6167735.0;%5108054.7;%3102422.8;%5916780.3;%5118046.2;%4588493.9;%1039781.2;%1719122.2;%4588493.9;%2024019.2;%4999973.4;%3627464.8; %4098389.4;%4660117.7;%2788353.3; %4437965.6;%3074985.3;%3449082.23;%2319588.6; %4712070.3;%2658208.7;%4712070.3;%2658208.7;%1919907; %3309766.6; %a.u.^2/s
datWndVr.dualCst.Dex2 = [];
datVr.dualCst.cFreq2 =  3179.5;%3425.3;%3134.7;%3272.2;%2434.7;%3134.7;%;%2949.0;%3109.5;%2941.6;%1753.4;%2804.7;%3127.3;%2882.4;%4619.6;%5000.6;%4402.9;%2824.4;%3649.2;%3777.2;%3233.9;%3777.2;%3319.6;%3408.6; %4791.2;%3341.9;%4159.6; %3451.6; %4061.2; %1500;%1749.6; %1376.9;%2182.8; %Hz
datWndVr.dualCst.cFreq2 = [];
datVr.dualCst.R2 = 0.0000011;%0.0000011; %m
datWndVr.dualCst.R2= [];
datVr.dualCst.gam2 = 6*pi*datVr.dualCst.eta*datVr.dualCst.R2; %N*s/m
datWndVr.dualCst.gam2 = [];
datVr.dualCst.k2 = 2*pi*datVr.dualCst.gam2*datVr.dualCst.cFreq2*10^6; %pN/um
datWndVr.dualCst.k2 = [];
datVr.dualCst.Dth2 = datVr.dualCst.kB*datVr.dualCst.T/datVr.dualCst.gam2; %m^2/s
datWndVr.dualCst.Dth2 = [];
datVr.dualCst.alpha2 = sqrt(datVr.dualCst.Dth2/datVr.dualCst.Dex2)*10^6; %um/a.u.
datWndVr.dualCst.alpha2 = [];

%Fit Parameters Trap 2 - Y
datWndLabVr.trap4Cst = [];
datVr.dualCst.Dex4 =  1344605.3;
datWndVr.dualCst.Dex4 =  [];
datVr.dualCst.cFreq4 = 2724.1; %Hz
datWndVr.dualCst.cFreq4 = [];
datVr.dualCst.R4 = 0.0000011; %m
datWndVr.dualCst.R4 = [];
datVr.dualCst.gam4 = 6*pi*datVr.dualCst.eta*datVr.dualCst.R4; %N*s/m
datWndVr.dualCst.gam4 = [];
datVr.dualCst.k4 = 2*pi*datVr.dualCst.gam4*datVr.dualCst.cFreq4*10^6; %pN/um
datWndVr.dualCst.k4 = [];
datVr.dualCst.Dth4 =  datVr.dualCst.kB*datVr.dualCst.T/datVr.dualCst.gam4; %m^2/s
datWndVr.dualCst.Dth4 = [];
datVr.dualCst.alpha4 = sqrt(datVr.dualCst.Dth4/datVr.dualCst.Dex4)*10^6; %um/a.u.
datWndVr.dualCst.alpha4 = [];
%% Plotting Variables


pltVr.legOn = 1; %Legend toggle 1/0;
pltVr.legArry = []; %Contains legend labels.
pltVr.rstAx = 1; %Axis reset toggle to default axis;
pltVr.ftAx = []; %Contains current Force-Time axis;
pltVr.teAx = []; %Contains current Time-Extension axis;
pltVr.feAx = [0.5 1.5 -10  85]; %Contains current Force-Extension axis;
pltVr.ssAx = []; %Contains current Sensitivity axis;
pltVr.grphAct = 0; %Currently active graph; 1 = Force-Time; 2 = Extension-Time; 3 = Force-Extension; 4 = Sensitivity;
pltVr.rulVal = 2; %This value controls the ruler used.
pltVr.fitVal = [0 1]; %Fitting bounds where 0 is the most left bound and 1 is the most right bound.
pltVr.lIdx = []; %Selected regions left boundaries.
pltVr.rIdx = []; %Selected regions right boundaries.
pltVr.unfDet = 0; %Unfolding detection toggle 1/0;
pltVr.WLCfit = 0; %WLC fitting enabled;
pltVr.LpPlot = 0; %Plots Contour Length;
pltVr.splitPulls = 1; %Toggles splitting Extension/Relaxation Curves.
pltVr.figs = []; %Stores figure handles.
%% SIGNAL COMBINATION CALIBRATION - Used for single trap tweezers setup to correct for drift.
doCal = false;
if doCal == true    

h = waitbar(0,'Performing Signal Combination Calibration');
waitbar(0,h);
% To use this program replace currentpath, with path to where the folder is
% saved. Then open 'values sx and start' in Matlab and run >> Signalcombination(values)
currentpath = pwd;
values = [3 4.80000000000000 1610;4 4.40000000000000 1350;5 4.10000000000000 1600;6 4.10000000000000 3100;7 4 1950;8 4.30000000000000 2700;9 4.80000000000000 2500;10 4.40000000000000 1900;11 5.20000000000000 2300;12 4.20000000000000 2400;13 3.80000000000000 1700;14 4.50000000000000 3200;15 5.30000000000000 1;17 4.50000000000000 1810;18 4 2400;19 4.10000000000000 2000;20 4.30000000000000 1400;21 5.20000000000000 2500;22 5.10000000000000 3600;23 4.20000000000000 3350];

% Program makes distributions based on recorded data, sampling frequency is
% currently 50 Hz, and uses these distributions for drift correction.

% first column is p_CCD_pipette, second column is p_QPD_pipette
lmax = 10050;
distributionsraw = zeros(lmax, 3);
% sum up all 20 traces and interpolate linearly, if trace is too short.
for i = 1:20
    
    n = values(i,1);
    if n>9
        n = ['0' num2str(n)];
    else
        n = ['00' num2str(n)];
    end
    filename = ['20110722' n '.dat'];
    pathname = fullfile(currentpath, 'Signalcombination/20110722', filename);
    fid = fopen(strrep(pathname, '\', filesep), 'r');
    
    dataTemp = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f','headerlines', 1);
    fclose(fid);
    xpz = dataTemp{2};
    xspt1 = dataTemp{5};
    xspt2 = dataTemp{7};
    vx = dataTemp{9};
    
    ls = length(xpz);
    
    % find sensitivity
    c = polyfit(xspt1(values(i,3):ls).*0.0615, vx(values(i,3):ls), 1);
    sx = c(1);
    
    % The voltage signal is comparatively noise free and we use it here as
    % true signal. + sign, since vx is already negative compared to video
    % tracking
    n1 = (xspt1(values(i,3):ls)).*0.0615 - vx(values(i,3):ls)./sx; % video tracking one bead
    n2 = vx(values(i,3):ls)./sx; % since piezo is not moved.
    n3 = (xspt1(values(i,3):ls) - xspt2(values(i,3):ls)).*0.0615; % noise in differential signal
    
    ln = length(n1);
    fn1 = fft(n1 - mean(n1));
    fn2 = fft(n2 - mean(n2));
    fn3 = fft(n3 - mean(n3));
    
    lp = floor(ln./2);
    
    %pn1 is p_CCD, pn2 is p_pipette
    pn1 =(abs(fn1(1:lp))).^2./ln;
    pn2 =(abs(fn2(1:lp))).^2./ln;
    pn3 =(abs(fn3(1:lp))).^2./ln;
    
    pnew1 = zeros(lmax,1);
    pnew2 = zeros(lmax,1);
    pnew3 = zeros(lmax,1);
    
    for j = 1:lp
        k = ceil(j/lp*lmax);
        pnew1(k) = pn1(j);
        pnew2(k) = pn2(j);
        pnew3(k) = pn3(j);
    end
    % now remove zeros in pnew1 by linear interpolation
    g = 0;
    if pnew1(1) == 0
        while(pnew1(1+g)==0)
            g = g+1;
        end
        for m = 1:g
            pnew1(m) = pnew1(g+1);
        end
    end
    
    for t = 1:lmax
        g = 0;
        if pnew1(t) == 0
            while(pnew1(t+g)==0)
                g = g+1;
                if t+g > lmax
                    g = g-1;
                    break
                end
            end
            for m = 0:g-1
                pnew1(t+m) = (pnew1(t+g)-pnew1(t-1))*(m)/(g)+pnew1(t-1);
            end
        end
    end
    
    if pnew1(lmax) == 0
        pnew1(lmax) = pnew1(lmax-1);
    end
    
    % now the same for pnew2
    g = 0;
    if pnew2(1) == 0
        while(pnew2(1+g)==0)
            g = g+1;
        end
        for m = 1:g
            pnew2(m) = pnew2(g+1);
        end
    end
    
    for t = 1:lmax
        g = 0;
        if pnew2(t) == 0
            while(pnew2(t+g)==0)
                g = g+1;
                if t+g > lmax
                    g = g-1;
                    break
                end
            end
            for m = 0:g-1
                pnew2(t+m) = (pnew2(t+g)-pnew2(t-1))*(m)/(g)+pnew2(t-1);
            end
        end
    end
    
    if pnew2(lmax) == 0
        pnew2(lmax) = pnew2(lmax-1);
    end
    
    % now the same for pnew3
    g = 0;
    if pnew3(1) == 0
        while(pnew3(1+g)==0)
            g = g+1;
        end
        for m = 1:g
            pnew3(m) = pnew3(g+1);
        end
    end
    
    for t = 1:lmax
        g = 0;
        if pnew3(t) == 0
            while(pnew3(t+g)==0)
                g = g+1;
                if t+g > lmax
                    g = g-1;
                    break
                end
            end
            for m = 0:g-1
                pnew3(t+m) = (pnew3(t+g)-pnew3(t-1))*(m)/(g)+pnew3(t-1);
            end
        end
    end
    
    if pnew3(lmax) == 0
        pnew3(lmax) = pnew3(lmax-1);
    end

    % now add new distritbutions
    
    distributionsraw(:,1) = distributionsraw(:,1) + pnew1;
    distributionsraw(:,2) = distributionsraw(:,2) + pnew2;
    distributionsraw(:,3) = distributionsraw(:,3) + pnew3;
    waitbar(i/26,h);
end


% block data to 25./q Hz resolution
q = 250;
p_QPD_pipette = zeros(q);
for i = 1:q
    p_QPD_pipette(i) = mean(distributionsraw((i-1)*floor(lmax./q)+1:floor(lmax./q)*i,2)./20);
end
p1 = distributionsraw(:,2)./20;
waitbar(21/26,h);
p_CCD_pipette = zeros(q);
for i = 1:q
    p_CCD_pipette(i) = mean(distributionsraw((i-1)*floor(lmax./q)+1:floor(lmax./q)*i,1)./20);
end
p2 = distributionsraw(:,1)./20;
waitbar(22/26,h);
p_differential = zeros(q);
for i = 1:q
    p_differential(i) = mean(distributionsraw((i-1)*floor(lmax./q)+1:floor(lmax./q)*i,3)./20);
end
p0 = distributionsraw(:,3)./20;
waitbar(23/26,h);
% Now make distribution of p_QPD without bead

filename = ['20110722' '000' '.dat'];
pathname = fullfile(currentpath, 'Signalcombination/20110722', filename);

fid = fopen(strrep(pathname, '\', filesep), 'r');

dataTemp = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f','headerlines', 1);
fclose(fid);

vx = dataTemp{9};
l2 = length(vx);

s2 = vx./4.5; % sensitivity is on average 4.5
fs2 = fft(s2 - mean(s2));
n2 = floor(l2./2);
p3 = abs(fs2(1:n2)).^2./l2;
p_QPD_nobead = zeros(q);
for i = 1:q
    p_QPD_nobead(i) = mean(p3((i-1)*floor(n2./q)+1: i*floor(n2./q)));
end
waitbar(24/26,h);
% Now make distribution of p_QPD with bead

filename = ['20110725' '000' '.dat'];
pathname = fullfile(currentpath, 'Signalcombination/20110725', filename);

fid = fopen(strrep(pathname, '\', filesep), 'r');

dataTemp = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f','headerlines', 1);
fclose(fid);

vx = dataTemp{9};
xvideo1 = dataTemp{5};
s3 = vx(1:57000)./4.5; % sensitivity is on average 4.5
l3 = length(s3);
fs3 = fft(s3 - mean(s3));
n3 = floor(l3./2);
p4 = abs(fs3(1:n3)).^2./l3;
p_QPD_bead = zeros(q);
for i = 1:q
    p_QPD_bead(i) = mean(p4((i-1)*floor(n3./q)+1: i*floor(n3./q)));
end
waitbar(25/26,h);
% And distribution of video tracking noise for bead in trap
s4 = xvideo1(1:57000).*0.0615;
l4 = length(s3);
fs4 = fft(s4 - mean(s4));
n4 = floor(l4./2);
p5 = abs(fs4(1:n4)).^2./l4;
p_CCD_trap = zeros(q);
for i = 1:q
    p_CCD_trap(i) = mean(p5((i-1)*floor(n4./q)+1: i*floor(n4./q)));
end
waitbar(26/26,h);
% now make p_CCD
p_CCD = p_CCD_pipette + p_CCD_trap;

%now make p_QPD_piezo
p_QPD_piezo = p_QPD_pipette + p_QPD_nobead;

close(h);
else
   
    load('calVal.mat');
    
end
%% DRAW ANALYSIS WINDOW - Various GUI Elements are created both for the single as the dual trap tweezers, such as buttons, sliders, toggles etc.

if exist('progState.mat', 'file') == 2    
    load('progState.mat','prgVr');    
else
    save('progState.mat','PrgInit');   
end
global wndVr;

wndVr.monPos = get(0, 'MonitorPositions');

if prgVr.setupStt == 0
       
    wndVr.hTrig = 0;
    wndVr.vTrig = 0;
    wndVr.vSlider = [];
    wndVr.hSlider = [];
    wndVr.sampling_shift = [];
    wndVr.sampling_xlim = [];
    wndVr.sampling_ylim = [];
    wndVr.hpc = [];
    wndVr.botPanel = [];
    wndVr.axesPos = [];
    wndVr.forBut = [];
    wndVr.sampling_for = [];
    wndVr.sensBut = [];
    wndVr.sampling_sens = [];
    wndVr.hruler = [];
    wndVr.hsignal = [];
    wndVr.pt = [];
    wndVr.ht = [];
    wndVr.vt = [];
    wndVr.st = [];
    wndVr.sampling_pol = [];
    wndVr.a = [];
    wndVr.legendBut = [];
    wndVr.resetAxBut = [];
    wndVr.unfBut = [];
    wndVr.WLCBut = [];
    wndVr.contourPlotBut = [];
    wndVr.biasCalcBut = [];
    wndVr.fluoPanel = [];
    wndVr.axesPosFluo = [];
    wndVr.b = [];
    drawTweezersWindow()
    
elseif prgVr.setupStt == 1
    wndVr.hTrig = 0;
    wndVr.vTrig = 0;
    wndVr.vSlider = [];
    wndVr.hSlider = [];
    wndVr.sampling_shift = [];
    wndVr.sampling_xlim = [];
    wndVr.sampling_ylim = [];
    wndVr.hpc = [];
    wndVr.botPanel = [];
    wndVr.axesPos = [];
    wndVr.forBut = [];
    wndVr.sampling_for = [];
    wndVr.sampling_const = [];
    wndVr.sensBut = [];
    wndVr.sampling_sens = [];
    wndVr.hruler = [];
    wndVr.hsignal = [];
    wndVr.pt = [];
    wndVr.ht = [];
    wndVr.vt = [];
    wndVr.st = [];
    wndVr.sampling_pol = [];
    wndVr.a = [];
    wndVr.legendBut = [];
    wndVr.resetAxBut = [];
    wndVr.unfBut = [];
    wndVr.WLCBut = [];
    wndVr.contourPlotBut = [];
    wndVr.biasCalcBut = [];
    wndVr.fluoPanel = [];
    wndVr.axesPosFluo = [];
    wndVr.b = [];
    drawFoldometerWindow()
    
elseif prgVr.setupStt == 2
    
    wndVr.hTrig = 0;
    wndVr.vTrig = 0;
    wndVr.vSlider = [];
    wndVr.hSlider = [];
    wndVr.sampling_shift = [];
    wndVr.sampling_xlim = [];
    wndVr.sampling_ylim = [];
    wndVr.hpc = [];
    wndVr.botPanel = [];
    wndVr.axesPos = [];
    wndVr.forBut = [];
    wndVr.sampling_for = [];
    wndVr.sampling_const = [];
    wndVr.sensBut = [];
    wndVr.sampling_sens = [];
    wndVr.hruler = [];
    wndVr.hsignal = [];
    wndVr.pt = [];
    wndVr.ht = [];
    wndVr.vt = [];
    wndVr.st = [];
    wndVr.sampling_pol = [];
    wndVr.a = [];
    wndVr.legendBut = [];
    wndVr.resetAxBut = [];
    wndVr.unfBut = [];
    wndVr.WLCBut = [];
    wndVr.contourPlotBut = [];
    wndVr.biasCalcBut = [];
    wndVr.fluoPanel = [];
    wndVr.axesPosFluo = [];
    wndVr.b = [];
    drawCTrapWindow()
    
end   