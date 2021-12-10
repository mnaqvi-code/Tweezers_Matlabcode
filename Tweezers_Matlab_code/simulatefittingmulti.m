clear all;
close all;
clc;

global want_hydro want_alias want_filter elec_filters fNyq delta_t sampling_f

want_hydro = 0;
want_alias = 0;
want_filter = 0;
elec_filters = [];
% Analysis

fSampAn = [50,100 150 200 250];
DAn = [2.050899925429805e-13 2.050899925429805e-13 2.050899925429805e-13 2.050899925429805e-13 2.050899925429805e-13];
DdetAn = [4.411253215116754e-13 2.422462969832773e-13 2.209831644483197e-13 2.145519004958776e-13 2.090865309226233e-13];
NwindowAn = [32 32 23 23 23];
NSampAn = [800 800 800 800 800];
NbinAn = [512 512 512 512 512];

figure;
plot(fSampAn,DAn,fSampAn,DdetAn);
legend('Real Diffusion Coeff','Estimated Diffusion Coeff')
xlabel('Frequency (kHz)');
ylabel('Diffusion Coeff');
%Papers used

% [1]
% Power spectrum analysis with least-squares fitting: Amplitude bias and
% its elimination, with application to optical tweezers and atomic force
% microscope cantilevers. Nørrelykke et al.
% Review of Scientific Instruments 81, 075103 (2010); doi: 10.1063/1.3455217
% [2]
% Power spectrum analysis for optical tweezers. Berg-Sørensen et al.
% Review of Scientific Instruments 75, 594 (2004); doi: 10.1063/1.1645654

%Important Parameters
fSamp = 150000;%16384*16;%16384*16; %Hz Sampling frequency
Nsamp = 800000;%262144*4;%262144*4; %N Number of samples
Nsplit = 32; %N Number of powerspectra to split into

%Constants

dT = 1/fSamp; %s
kB = 1.38064852*10^(-23); %J/K

kapList = [100:100:1500].*1E-6;
rPlist = [0.5 : 0.1: 2.0].*1E-6;

T = 294.15; %K
for kk = 1:length(kapList)
    for gg = 1:length((rPlist))
        close all;
        clear xv deltXV Tmax f P FT
        kap = kapList(kk);%1000E-6; %N/m;
        rP = rPlist(gg);%1.05E-6; %m
        gam = 6*pi*0.0010005*rP;%
        D = kB*T/gam;
        vP = (4/3)*pi*rP^3; %m^3
        rhoP = 1.04; %kg/m^3
        mP = vP*rhoP; %kg
        fc = kap/(2*pi*gam);
        
        %Implementing Appendix A from [1]
        
        lP = (gam/(2*mP)) + sqrt(gam^2/(4*mP^2)-kap/mP); %Lambda+
        lM = (gam/(2*mP)) - sqrt(gam^2/(4*mP^2)-kap/mP); %Lambda-
        cP = exp(-lP*dT); %c+
        cM = exp(-lM*dT); %c-
        
        expMT = 1/(lP-lM).*[-lM*cP+lP*cM, -cP+cM; lP*lM*(cP-cM), lP*cP-lM*cM ]; %Exp(-M*deltaT)
        
        aP = (lP+lM)/(lP-lM)*sqrt(((1-cP^2)*D)/(2*lP)); %A+
        aM = (lP+lM)/(lP-lM)*sqrt(((1-cM^2)*D)/(2*lM)); %A-
        alp = 2*sqrt(lP*lM)/(lP+lM)*(1-cP*cM)/sqrt((1-cP^2)*(1-cM^2)); %alpha
        
        %xv(1,1) = 0;
        %xv(2,1) = 0;
        
        xv = zeros(2,Nsamp);
        %
        % h = waitbar(0,'Initializing waitbar...');
        % set(h,'Name','Simulating Brownian Motion');
        % tTot = 0;
        % tCount = 0;
        
        %Performing Monte Carlo Simulation of the Einstein-Ornstein-UhlenBeck
        %Theory of Brownian Motion in a Harmonic Potential
        
        %tic;
        deltXV = (aP.*[-1; lP]+aM.*[1; -lM]).*sqrt(1+alp).*normrnd(0,1,[1,length(2:Nsamp)])+(aP.*[-1; lP]-aM.*[1; -lM]).*sqrt(1-alp).*normrnd(0,1,[1,length(2:Nsamp)]);
        
        %xvtest(1,:) = filter(1,[1 -(expMT(1,1)+expMT(1,2))],xv(1,:),[0 deltXV(1,:)],1);
        %xvtest(1:2,:) = filter([1; 1],[[1; 1] [-(expMT(1,1)+expMT(1,2)); -(expMT(2,1)+expMT(2,2))]],xv(1:2,:),[[0; 0] deltXV(1:2,:)]);
        for ii = 2:Nsamp
            %deltXV = (aP.*[-1; lP]+aM.*[1; -lM]).*sqrt(1+alp).*normrnd(0,1)+(aP.*[-1; lP]-aM.*[1; -lM]).*sqrt(1-alp).*normrnd(0,1);
            xv(1:2,ii) = expMT*xv(1:2,ii-1)+deltXV(:,ii-1);
            %     if mod(ii,floor(Nsamp/10)) == 0
            %         t = toc;
            %         tCount = tCount + 1;
            %         tTot = t + tTot;
            %         hours = floor((tTot/tCount)*(length(2:Nsamp)/(Nsamp/1000)-tCount)/3600);
            %         mins = floor(((tTot/tCount)*(length(2:Nsamp)/(Nsamp/1000)-tCount)-hours*3600)/60);
            %         secs = (tTot/tCount)*(length(2:Nsamp)/(Nsamp/1000)-tCount)-hours*3600-mins*60;
            %         waitbar((ii)/(length(2:Nsamp)),h,[num2str(round(hours)) ':' num2str(round(mins)) ':' num2str(round(secs)) ' remaining.'])
            %         tic;
            %     end
        end
        
        
        %max(abs(xvtest(1,:)-xv(1,:)))
        % close(h)
        %keyboard;
        
        % load('beadData.mat');
        % xv = Data{1}';
        % fSamp = 50000;
        
        
        figure; %Plot simulated brownian motion.
        plot(dT:dT:length(xv(1,:))*dT,xv(1,:).*1E6);
        
        fNyq    =   fSamp / 2;
        delta_t =   dT;%1 / fSamp;
        xv(:,end-mod(length(xv(1,:)),Nsplit)+1:end) = [];
        if mod(length(xv(1,:))/Nsplit,2) ~= 0
            xv(:,end-Nsplit+1:end) = [];
        end
        for ii = 1:Nsplit
            
            time{ii}    =   [ceil(length(xv(1,:))*(ii-1)/Nsplit)*delta_t : delta_t : (ceil(length(xv(1,:))*ii/Nsplit)-1)*delta_t]';
            Tmax{ii}       =   max(time{ii})-min(time{ii});
            f{ii}       =   (([ceil(length(xv(1,:))*(ii-1)/Nsplit) : floor(length(xv(1,:))*ii/Nsplit)+1] - ceil(length(xv(1,:))*(ii-1)/Nsplit)) / (Tmax{ii}))';
            
            FT{ii}      =   delta_t*fft(xv(1,ceil(length(xv(1,:))*(ii-1)/Nsplit)+1 : floor(length(xv(1,:))*ii/Nsplit)));
            P{ii}       =   FT{ii} .* conj(FT{ii}) / Tmax{ii};
            
            ind     =   find(f{ii} <= fNyq); % only to the Nyquist f
            f{ii}       =   f{ii}(ind);
            P{ii}       =   P{ii}(ind)';
            
        end
        
        
        ftot = mean([f{:}],2);
        Ptot = mean([P{:}],2);
        
        
        
        %Ptot = Ptot .*(1 ./ (1 + (ftot./(fSamp./2.0833)).^128));
        
        logf=log(ftot);
        
        nbin = fSamp/(64*8);
        turn = 1;
        color_x     =   [0.6 0.1 0.8];
        color_y     =   [0.2 0.7 0.3];
        color = color_x;
        
        hold on;
        cP = exp(-pi*fc/fNyq);
        deltXX = sqrt((1-cP^2)*D/(2*pi*fc));
        figure; % plot Aliased lorentzian (see formula 23 in [2]).
        plot(1:fSamp,deltXX.^2.*dT./(1+cP.^2-2.*cP.*cos(2.*pi.*(1:fSamp)./fSamp)));
        hold on;
        plot(ftot,(D./(2.*pi.^2))./(fc^2+ftot.^2));
        
        [count, centers] = hist(logf,nbin);
        delta   = centers(2) - centers(1);
        clear fb_plot Pb_plot s_plot
        if  (turn == 1)
            for i = 1 : nbin
                indp(i).indx  = find(logf >= centers(i)-delta/2 & logf < centers(i)+delta/2);
            end
        end
        
        for i = 1:nbin
            fb_plot(i)  = mean_nan(ftot(indp(i).indx));
            Pb_plot(i)  = mean_nan(Ptot(indp(i).indx));
            if sum(isfinite(Ptot(indp(i).indx))) > 0, s_plot(i) = Pb_plot(i)/sqrt(sum(isfinite(Ptot(indp(i).indx)))); else s_plot(i) = NaN; end;
        end;
        
        set(gca,'XLim',[10 25000]);%max(ftot)]);
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
        
        % Plot binned powerspectrum.
        
        
        
        hd = plot(fb_plot, Pb_plot, 'k.'); set(hd,'MarkerSize',10,'MarkerFaceColor',color,'MarkerEdgeColor',color);
        hold off;
        h = xlabel('Frequency (Hz)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
        h = ylabel(['P(f) (arbitrary units) ']); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
        set(gca,'XScale','log','YScale','log','FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
        title('Binned power spectrum with binsize = 512, Fsamp = 50 kHz, NSamp = 800k, NWindow = 32 2x Bead')
        legend('Aliased Lorentzian','Lorentzian','Binned Data');
        ind = find(ftot <= 25000);
        ftotdet = ftot(ind);
        Ptotdet = Ptot(ind);
        
        for p = 0 : 2
            for q = 0 : 2
                eval(['s' num2str(p) num2str(q) ' = sum ((ftotdet .^ (2*' num2str(p) ')) .* (Ptotdet .^ (' num2str(q) ')));']);
            end;
        end;
        
        fcDet = sqrt((s01*s22-s11*s12)/(s11*s02-s01*s12));
        DDet = (2*pi^2*Nsplit)/(Nsplit+1)*(s02*s22-s12^2)/(s11*s02-s01*s12);
        
        figure;
        plot(ftot,Ptot);
        hold on;
        plot(1:fSamp,deltXX.^2.*dT./(1+cP.^2-2.*cP.*cos(2.*pi.*(1:fSamp)./fSamp)),'LineWidth',3);
        set(gca,'XLim',[10 25000]);%max(ftot)]);
        h = xlabel('Frequency (Hz)'); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
        h = ylabel(['P(f) (arbitrary units) ']); set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
        set(gca,'XScale','log','YScale','log','FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        % S.D.
        
        x   = min(ftot) / fc;
        y   = max(ftot) / fc;
        s   = sqrt(pi) * ( (2*y) / (1 + y^2) - (2*x) / (1 + x^2) + 2 * atan((y - x) / (1 + x*y)) - ...
            (4/(y - x)) * (atan((y - x) / (1 + x*y)))^2) ^ (-1/2);
        
        sfc = s * fc / sqrt(pi * fc * T);
        
        g   = sqrt( ((2*y)/(1 + y^2)-(2*x)/(1 + x^2) + 2*atan((y - x) / (1 + x*y)) )/((1 + pi/2)*(y - x)) );
        
        sD  = D * sqrt( (1 + pi/2) / (pi * fc * T) )*g*s;
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        % Bin the powerspectrum
        nblock = 100;
        nbin    =   floor(length(ftot)/nblock);
        fb = []; Pb = []; s = [];
        for i = 1 : nbin
            fb(i)   = mean_nan(ftot((i-1)*nblock+1 : i*nblock));
            Pb(i)   = mean_nan(Ptot((i-1)*nblock+1 : i*nblock));
            s(i)    = (1/Pb(i))/sqrt(sum(isfinite(Ptot((i-1)*nblock+1 : i*nblock))));
        end;
        fb = fb(:); Pb = Pb(:); s = s(:);
        a   = (s01 * s22 - s11 * s12) / (s02 * s22 - s12.^2);
        b   = (s11 * s02 - s01 * s12) / (s02 * s22 - s12.^2);
        
        %fc  = sqrt(a/b);
        %D   = (1/b) * 2 * (pi.^2) * Nsplit/(Nsplit+1);
        
        Pfit = 1 ./ (a + b .* ftot.^2);
        
        parameters0 = [fcDet,DDet,[],[]];
        scal_fit = ones(1,length(parameters0));         %Scaled fitting parameters
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  Choose data to be fit by a Lorentzian
        Ffit_start = 2;
        Ffit_end = 20000;
        ind     = find(fb > Ffit_start & fb <= Ffit_end);
        xfin    = fb(ind);
        yfin    = Pb(ind);
        sfin    = s(ind);
        tolx = 10E-99;
        check_flag = 0;
        sampling_f = fSamp;
        
        weighted_difference = @funn; % Creating a function handle better than defining a string (see matlab help on feval)
        
        sampling_f = sampling_f/1000;
        [scal_fit,RESNORM,RESIDUAL,JACOBIAN] = fit_nonl(weighted_difference,scal_fit,tolx,1500,parameters0,xfin,yfin,sfin,check_flag);
        sampling_f = sampling_f*1000;
        scal_fit = abs(scal_fit);               % The function P_theor is symmetric in alpha and fdiode
        parameters = scal_fit.*parameters0;
        color_now = color_x;
        load('testSpectrum.mat')
        %ftot = f;
        %Ptot = P/10^11.1;
        P = P./(1 ./ (1 + (f./(15000)).^2));
        %plot(f,P/10^11)
        plot(ftot,(D./(2.*pi.^2))./(fc^2+ftot.^2),'LineWidth',3);
        plot(ftot,((DDet*scal_fit(2))./(2.*pi.^2))./((fcDet*scal_fit(1))^2+ftot.^2),'LineWidth',2);
        plot(ftot,(1-(exp(-pi*fcDet*scal_fit(1)/fNyq))^2)*DDet*scal_fit(2)./(2*pi*fcDet*scal_fit(1)).*delta_t./(1+(exp(-pi*fcDet*scal_fit(1)/fNyq)).^2-2.*(exp(-pi*fcDet*scal_fit(1)/fNyq)).*cos(2.*pi.*(ftot)./(sampling_f))),'LineWidth',2)
        
        
        legend('Windowed Data','Aliased Lorentzian','Real Lorentzian','Fitted Lorentzian','Fitted Aliased Lorentzian');
        
        figure;
        loglog(ftot,Ptot./decimate(P,2));
        figure;
        hist(Ptot./decimate(P,2));
        
        Dreal(kk,gg) = D;
        Dest(kk,gg) = DDet;
        Dfit(kk,gg) = DDet*scal_fit(2);
        kk
        gg
        
    end
end

figure; imagesc(rPlist.*1E6,kapList.*1E6,Dfit./Dreal.*100); colormap hot; colorbar
set(gca,'Fontsize',16,'FontWeight','bold')
xlabel('Bead Diameter (\mum)'); %set(h,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
ylabel('Trap Stiffness (pN \mum)');
axis square;