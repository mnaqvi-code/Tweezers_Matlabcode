
global pltVr wndVr

% FIT_POWERSPECTRUM fits the powerspectrum data vs. frequency.
% The function fit_nonl is used.

%%%%%%%%%%%%%%%%%%%%%%checks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform various checks to ensure that the input values are physically
%meaningful
check1(nblock,color1);
check(sampling_f,0,1000,'sampling frequency',color1);
check2(channel_x,1,size(File,2),'Channel X',color1);
check2(channel_y,0,size(File,2),'Channel Y',color1);
check2(channel_z,0,size(File,2),'Channel Z',color1);
check2(Lfit_start,0,fNyq,'Lorentzian fitting range',color1);
check(Lfit_end,Lfit_start,fNyq,'Lorentzian fitting range',color1);
check2(Ffit_start,0,fNyq,'Final fitting range',color1);
check(Ffit_end,Ffit_start,fNyq,'Final fitting range',color1);
check2(Plot_start,1,fNyq,'Plotting range',color1);
check(Plot_end,Plot_start,fNyq,'Plotting range',color1);

if (want_hydro == 1)
    check_hydro(R,0,'Bead radius',color1);
    check_hydro(l,0,'Distance to surface',color1);
    check_hydro(rho,0,'Density of bead',color1);
    check_hydro(rhol,0,'Density of fluid',color1);
    check_hydro(nu,0,'Kinematic viscosity',color1);
end

for jj = 1:length(Px)
    
    if isfield(pltVr.figs, 'fitCst') && jj == 1
        if iscell(pltVr.figs.fitCst)
            for kk = 1:length(pltVr.figs.fitCst)
                if ~isempty(findobj('name',['Consistency of Fit in Trap ', num2str(kk)]))
                    close(pltVr.figs.fitCst{kk})
                end
            end
        else
            if ~isempty(findobj('name',['Consistency of Fit in Trap ', num2str(kk)]))
                close(pltVr.figs.fitCst{jj})
            end
        end
        rmfield(pltVr.figs, 'fitCst');
    end
    
    if isfield(pltVr.figs, 'fitFin') && jj == 1
        if iscell(pltVr.figs.fitFin)
            for kk = 1:length(pltVr.figs.fitFin)
                if ~isempty(findobj('name',['Final Fit, Y-', num2str(kk/2)]))
                    close(pltVr.figs.fitFin{kk})
                end
                if ~isempty(findobj('name',['Final Fit, X-', num2str(kk/2+0.5)]))
                    close(pltVr.figs.fitFin{kk})
                end
            end
        else
            if ~isempty(findobj('name',['Final Fit, Y-', num2str(jj)]))
                close(pltVr.figs.fitFin{jj})
            end
            if ~isempty(findobj('name',['Final Fit, X-', num2str(jj)]))
                close(pltVr.figs.fitFin{jj})
            end
        end
        rmfield(pltVr.figs, 'fitFin');
    end
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    if ~exist('Px') %view_data; end;
        htview  = uicontrol('Style','text', ...
            'Position',[left3+(wid1+edge) bot11 wid2 height], ...
            'String','Activate ''View data'' before clicking this button',...
            'BackgroundColor',color1);
        return
    else
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        % Fit X and Y
        
        for ixy = IXY % IXY = [1 2] if there exist 2 channels; IXY = 1 if there is only 1 channel
            if ixy == 1
                if prgVr.setupStt == 0 %% Old Setup
                    P = Px{jj}; Title = 'X';
                elseif prgVr.setupStt == 1 %% Foldometer
                    P = Px{jj}; Title = 'X';
                elseif prgVr.setupStt == 2 %% CTRAP
                    %P = Px{jj}./(1 ./ (1 + (f/15000).^2)); Title = 'X';
                    P = Px{jj}; Title = 'X';
                end
            else
                if prgVr.setupStt == 0 %% Old Setup
                    P = Py{jj}; Title = 'Y';
                elseif prgVr.setupStt == 1 %% Foldometer
                    P = Py{jj}; Title = 'Y';
                elseif prgVr.setupStt == 2 %% CTRAP
                    P = Py{jj}; Title = 'Y';%./(1 ./ (1 + (f/15000).^2));
                end
            end;
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            % Bin the powerspectrum
            nbin    =   floor(length(f)/nblock);
            fb = []; Pb = []; s = [];
            for i = 1 : nbin
                fb(i)   = mean_nan(f((i-1)*nblock+1 : i*nblock));
                Pb(i)   = mean_nan(P((i-1)*nblock+1 : i*nblock));
                s(i)    = (1/Pb(i))/sqrt(sum(isfinite(P((i-1)*nblock+1 : i*nblock))));
            end;
            fb = fb(:); Pb = Pb(:); s = s(:);
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  Choose data to be fit by a Lorentzian
            ind     = find((fb > Lfit_start & fb <= Lfit_end));
            xfin    = fb(ind);
            yfin    = Pb(ind);
            sfin    = s(ind);
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  First guess for the fitting parameters FC (corner frequency) and
            %  D (diffusion coefficient)
            [fc0,D0,sfc,sD,Pfit] = lorentz_analyt(xfin,yfin,nblock);
            
            disp(' ');
            disp(['Lorentzian: fc = ' num2str(fc0,'%5.2f') ' +- ' num2str(sfc,'%5.2f') ', D = ' num2str(D0,'%10.3e') ' +- ' num2str(sD,'%10.3e')]);
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  First guess for F_DIODE (3dB frequency of the photodiode)
            %xxxxxxxxxxxxxxxxx
            %  First, does the diode filter the data, i.e., should we fit f_diode ?
            if want_diode == 1
                if want_alias == 1
                    P_aliasedNyq = sum((D0/(2*pi^2)) ./ ((fNyq + 2*[-10:10]*fNyq).^2 + fc0^2));
                else
                    P_aliasedNyq = (D0/(2*pi^2)) / ( fNyq^2 + fc0^2);
                end; %want_alias == 1
                
                if Pb(length(Pb)) < P_aliasedNyq
                    dif         =   Pb(length(Pb)) / P_aliasedNyq;
                    f_diode0    =   sqrt(dif * fNyq^2 / (1 - dif));
                else f_diode0   =   2*fNyq;
                end;
            else f_diode0  =  [];
            end; % want_diode == 1
            
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  Choose data to be fit with the final fit
            
            ind     = find(fb > Ffit_start & fb <= Ffit_end);
            xfin    = fb(ind);
            yfin    = Pb(ind);
            sfin    = s(ind);
            read_filter_function;
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  Fit
            check_flag =0;
            weighted_difference = @funn; % Creating a function handle better than defining a string (see matlab help on feval)
            if (want_alpha == 1)
                alpha0 = 0.3;
                a0 = sqrt(1/alpha0^2 - 1);      % Substitute a0 for alpha to ensure that alpha lies between 0 and 1
            else
                alpha0 = [];
                a0 =[];
            end;
            parameters0 = [fc0 D0 f_diode0 a0];             %Initial fitting parameters
            scal_fit = ones(1,length(parameters0));         %Scaled fitting parameters
            disp(' ');
            if (length(parameters0) == 2)
                disp(['Parameters0: fc = ' num2str(parameters0(1),'%5.1f') ', D = ' num2str(parameters0(2),'%10.3e')]);
            elseif (length(parameters0) == 3)
                disp(['Parameters0: fc = ' num2str(parameters0(1),'%5.1f') ', D = ' num2str(parameters0(2),'%10.3e') ',  fdiode = ' num2str(parameters0(3),'%10.3e')]);
            elseif (length(parameters0) == 4)
                disp(['Parameters0: fc = ' num2str(parameters0(1),'%5.1f') ', D = ' num2str(parameters0(2),'%10.3e') ',  fdiode = ' num2str(parameters0(3),'%10.3e') ',    alpha = ' num2str(1/sqrt(1+parameters0(4)^2),'%5.3f')]);
            end
            [scal_fit,RESNORM,RESIDUAL,JACOBIAN] = fit_nonl(weighted_difference,scal_fit,tolx,50,parameters0,xfin,yfin,sfin,check_flag);
            scal_fit = abs(scal_fit);               % The function P_theor is symmetric in alpha and fdiode
            parameters = scal_fit.*parameters0;
            if (length(parameters) == 2)
                disp(['Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.3e')]);
            elseif (length(parameters) == 3)
                disp(['Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.3e') ',  fdiode = ' num2str(parameters(3),'%10.3e')]);
            else
                disp(['Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.3e') ',  fdiode = ' num2str(parameters(3),'%10.3e') ',    alpha = ' num2str(parameters(4),'%5.3f')]);
            end
            disp(['chi^2 = ',num2str(RESNORM,'%6.2f')]);
            nfree = length(yfin) - length(parameters0);
            bbac = 1. - gammainc(RESNORM/2.,nfree/2.);     %Calculate backing of fit
            chi2 = RESNORM/nfree;
            disp(['chi^2 per degree of freedom = ',num2str(chi2,'%6.2f') ', n_{free} = ' num2str(nfree,'%6.2f')]);
            disp(['backing = ',num2str(bbac,'%10.3e')]);
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  If nblock is small, add extra minimization criteria
            %         if nblock < 200
            %             parameters0 = parameters;
            %             scal_fit = ones(1,length(parameters0));
            %             sumfun =  inline('sum(((1./P_theor(scal_fit,parameters0,xfin,check_flag) - 1./yfin).^2) ./ (sfin.^2) + 2*log(P_theor(scal_fit,parameters0,xfin,check_flag)))',...
            %                 'scal_fit','parameters0','xfin','yfin','sfin','check_flag');
            %             [scal_fit,RESNORM,EXITFLAG,OUTPUT] = fminsearch(sumfun,scal_fit,...
            %                 optimset('Display','iter','MaxIter',100,'MaxFunEvals',1000,'LargeScale','off'),parameters0,xfin,yfin,sfin,check_flag);
            %             scal_fit = abs(scal_fit); % The function P_theor is symmetric in alpha and fdiode
            %             parameters = scal_fit.*parameters0;
            %             disp(' ');
            %             disp(['Second fit. Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.1e')]);
            %             if length(parameters0) > 2, disp([',  fdiode = ' num2str(parameters(3),'%10.1e')]); end;
            %             if length(parameters0) > 3, disp(['            alpha = ' num2str(1/sqrt(1+parameters(4)^2),'%5.2f')]); end;
            %             disp(['Minimized sum = ',num2str(RESNORM,'%6.1f')]);
            %             %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %             %   Calculate the real chi^2 and backing:
            %             [scal_fit,RESNORM,RESIDUAL,JACOBIAN] = fit_nonl(weighted_difference,scal_fit,tolx,0,parameters0,xfin,yfin,sfin,check_flag);
            %             nfree = length(yfin) - length(parameters0);
            %             bbac = 1. - gammainc(RESNORM/2.,nfree/2.);
            %             chi2 = RESNORM/nfree;
            %             disp(['chi^2 per degree of freedom = ',num2str(chi2,'%6.2f') ', n_{free} = ' num2str(nfree,'%6.2f')]);
            %             disp(['backing = ',num2str(bbac,'%8.0f')]);
            %         end;
            
            if (length(parameters0) > 3)
                parameters(4) = 1/sqrt(1+parameters(4)^2);
                for i= 1: length(parameters0)-1
                    JACOBIAN(:,i)=JACOBIAN(:,i)/parameters0(i);     %Rescaling Jacobian back to unscaled parameters
                end
                JACOBIAN(:,4) = JACOBIAN(:,4)/(-1/sqrt(1/alpha0-1)*1/sqrt(1/parameters(4)-1)*1/parameters(4)^2);     %Rescaling Jacobian back to unscaled parameters
            else
                for i= 1: length(parameters0)
                    JACOBIAN(:,i)=JACOBIAN(:,i)/parameters0(i);     %Rescaling Jacobian back to unscaled parameters
                end
            end
            eval(['scal_fit' Title ' = scal_fit;']);
            eval(['parameters' Title ' = parameters;']);
            eval(['parameters0' Title ' = parameters0;']);
            eval(['bac' Title ' = bbac;']);
            eval(['chi2' Title ' = chi2;']);
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  Errors on parameters
            CURVATURE   = JACOBIAN' * JACOBIAN;
            cov = COVARIANCES(CURVATURE);       %This way of inverting the CURVATURE matrix prevents a nearly singlular matrix in many cases
            SIGMA_PAR   = [];
            sigma_par = sigmapar(CURVATURE,parameters);
            eval(['sigma_par' Title ' = sigma_par;']);
            eval(['cov' Title ' = cov;']);
            disp('                                                                                              ')
            disp(['cov(fc,D)            = ',num2str(cov(1)/sqrt(parameters(1)*parameters(2)),'%8.3f')]);
            if  length(parameters0) > 2
                disp(['cov(fc,fdiode)       = ' num2str(cov(2)/sqrt(parameters(1)*parameters(3)),'%8.3f')]);
                disp(['cov(D,fdiode)        = ' num2str(cov(3)/sqrt(parameters(2)*parameters(3)),'%8.3f')]);
            end;
            if length(parameters0) > 3
                disp(['cov(fc,alpha)        = ' num2str(cov(4)/sqrt(parameters(1)*parameters(4)),'%8.3f')]);
                disp(['cov(D,alpha)         = ' num2str(cov(5)/sqrt(parameters(2)*parameters(4)),'%8.3f')]);
                disp(['cov(fdiode,alpha)    = ' num2str(cov(6)/sqrt(parameters(3)*parameters(4)),'%8.3f')]);
            end;
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %  Plot the powerspectrum
            
            if(size(wndVr.monPos,1) > 1)
                lef = wndVr.monPos(size(wndVr.monPos,1)-1,1); bot = wndVr.monPos(size(wndVr.monPos,1)-1,2); wid = wndVr.monPos(size(wndVr.monPos,1)-1,3); hei = wndVr.monPos(size(wndVr.monPos,1)-1,4);
            else
                screensize = get(0,'ScreenSize');
                P2 = get(gcf,'pos');
                lef = 0; bot = P2(2)+5; wid = screensize(3); hei = screensize(4) - P2(2)*2-20;
            end;
            
            if ixy == 1
                color_now = color_x;
                pltVr.figs.fitCst{jj} = figure; clf; subplot('position',[0.1 .6 .35 .3]); set(gcf,'Numbertitle','off','Name',['Consistency of Fit in Trap ', num2str(jj)]); hold on; h = title('Final Fit, X','FontWeight','bold'); ...
                    set(h,'FontUnits','Normalized','FontSize',0.07); plot_P_cos(f,P,scal_fit,parameters0,color_now);
                set(0, 'currentfigure', pltVr.figs.fitCst{jj}); shg; subplot('position',[0.1 .15 .35 .3]); hold on; h = title('Final Fit, X','FontWeight','bold');ind = find(f > Ffit_start & f < Ffit_end);...
                    set(h,'FontUnits','Normalized','FontSize',0.07);
                drawnow;
                pause(0.1);
                jFrame = get(handle(gcf),'JavaFrame');
                jFrame.setMaximized(true);
                plot_data_div_fit(f(ind),P(ind),scal_fit,parameters0,color_now);
                pltVr.figs.fitFin{jj*2-1} = figure; clf; set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name',['Final Fit, X-', num2str(jj)]); hold on; h = title('Final Fit, X','FontWeight','bold');ind = find(f > Plot_start & f < Plot_end);...
                    plot_fit(f(ind),P(ind),scal_fit,parameters0,color_now,1);
                drawnow;
                pause(0.1);
                jFrame = get(handle(gcf),'JavaFrame');
                jFrame.setMaximized(true);
            else
                color_now = color_y;
                set(0, 'currentfigure', pltVr.figs.fitCst{jj}); shg; subplot('position',[0.6 .6 .35 .3]); hold on; h = title('Final Fit, Y','FontWeight','bold');plot_P_cos(f,P,scal_fit,parameters0,color_now);
                set(h,'FontUnits','normalized','FontSize',0.07);
                set(0, 'currentfigure', pltVr.figs.fitCst{jj}); shg; subplot('position',[0.6 .15 .35 .3]); hold on; h = title('Final Fit, Y','FontWeight','bold');ind = find(f > Ffit_start & f < Ffit_end);...
                    set(h,'FontUnits','normalized','FontSize',0.07);
                drawnow;
                pause(0.1);
                jFrame = get(handle(gcf),'JavaFrame');
                jFrame.setMaximized(true);
                plot_data_div_fit(f(ind),P(ind),scal_fit,parameters0,color_now);
                pltVr.figs.fitFin{jj*2} = figure; clf; set(gcf,'Position',[lef bot wid hei],'Numbertitle','off','Name',['Final Fit, Y-', num2str(jj)]); hold on; h = title('Final Fit, Y','FontWeight','bold');ind = find(f > Plot_start & f < Plot_end);...
                    plot_fit(f(ind),P(ind),scal_fit,parameters0,color_now,2);
                drawnow;
                pause(0.1);
                jFrame = get(handle(gcf),'JavaFrame');
                jFrame.setMaximized(true);
            end;
            
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            
        end; %(for ixy = [1 2])
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %Create caption for figures
        if jj == 1
            clear writeVal;
            writeVal(1) = FileInfo.datenum;
        end
        set(0, 'currentfigure', pltVr.figs.fitFin{jj*2-1}); shg;
        STR ={};
        STR2 = {};
        STR(1) = {'               '};        
        if (length(parameters) == 2),    STR2(1) = {['Cov(f_c,D)    = ',num2str(covX(1)/sqrt(parametersX(1)*parametersX(2)),'%8.3f')]};end;
        STR(2) = {['f_c (Hz) = ' num2str(parametersX(1),'%0.1f') ' \pm ' num2str(sigma_parX(1),'%0.1f')]};        
        writeVal(4*jj-3+1) = parametersX(1);
        STR(3) = {['D (arb. units)^2/s = ' num2str(parametersX(2),'%0.1f') ' \pm ' num2str(sigma_parX(2),'%0.1f')]};
        writeVal(4*jj-2+1) = parametersX(2);
        if length(parameters) > 2
            STR(4) = {['f_{3dB,diode} (Hz) = ' num2str(parametersX(3),'%0.1f') ' \pm ' num2str(sigma_parX(3),'%0.1f')]};
            
            if length(parameters) > 3
                STR(5) = {['\alpha_{diode}   = ' num2str(parametersX(4),'%0.1f') ' \pm ' num2str(sigma_parX(4),'%0.1f')]};
                STR(6) = {['\chi^2 per degree of freedom = ' num2str(chi2X,'%8.2f')]};
                STR(7) = {['Backing  = ' num2str(100*bacX,'%8.0f') ' %']};
                STR2(2) = {['Cov(f_c,D)             = ',num2str(covX(1)/sqrt(parametersX(1)*parametersX(2)),'%0.3f')]};
                STR2(3) = {['Cov(f_c,f_{diode})        = ' num2str(covX(2)/sqrt(parametersX(1)*parametersX(3)),'%0.3f')]};
                STR2(4) = {['Cov(fc,\alpha)              = ' num2str(covX(4)/sqrt(parametersX(1)*parametersX(4)),'%0.3f')]};
                STR2(5) = {['Cov(D,f_{diode})        = ' num2str(covX(3)/sqrt(parametersX(2)*parametersX(3)),'%0.3f')]};
                STR2(6) = {['Cov(D,\alpha)             = ' num2str(covX(5)/sqrt(parametersX(2)*parametersX(4)),'%0.3f')]};
                STR2(7) = {['Cov(f_{diode},\alpha)        = ' num2str(covX(6)/sqrt(parametersX(3)*parametersX(4)),'%0.3f')]};
            else
                STR(5) = {['\chi^2 per degree of freedom = ' num2str(chi2X,'%8.2f')]};
                STR(6) = {['Backing  = ' num2str(100*bacX,'%8.0f') ' %']};
                STR2(2) = {['Cov(f_c,D)         = ',num2str(covX(1)/sqrt(parametersX(1)*parametersX(2)),'%0.3f')]};
                STR2(3) = {['Cov(f_c,f_{diode}) = ' num2str(covX(2)/sqrt(parametersX(1)*parametersX(3)),'%0.3f')]};
                STR2(4) = {['Cov(D,f_{diode})   = ' num2str(covX(3)/sqrt(parametersX(2)*parametersX(3)),'%0.3f')]};
            end;
        else
            STR(4) = {['\chi^2 per degree of freedom = ' num2str(chi2X,'%8.2f')]};
            STR(5) = {['Backing = ' num2str(100*bacX,'%8.0f') ' %']};
        end;
        caption(STR,0.1,0.15);
        caption(STR2,0.65,0.15);
        writeVal(4*jj-1+1) = 1;
        writeVal(4*jj+1) = 1;
        if (length(IXY) >1)
            
            
            set(0, 'currentfigure', pltVr.figs.fitFin{jj*2}); shg;
            STR ={};
            STR2 = {};
            if (length(parameters) == 2), STR2(1) = {['cov(f_c,D)           = ',num2str(covY(1)/sqrt(parametersY(1)*parametersY(2)),'%8.3f')]}; end;
            STR(1) = {'               '};
            STR(2) = {['f_c (Hz) = ' num2str(parametersY(1),'%0.1f') ' \pm ' num2str(sigma_parY(1),'%0.1f')]};
            writeVal(4*jj-1+1) = parametersY(1);
            STR(3) = {['D (arb. units)^2/s = ' num2str(parametersY(2),'%0.1f') ' \pm ' num2str(sigma_parY(2),'%0.1f')]};
            writeVal(4*jj+1) = parametersY(2);
            if length(parameters) > 2
                STR(4) = {['f_{3dB,diode} (Hz) = ' num2str(parametersY(3),'%0.1f') ' \pm ' num2str(sigma_parY(3),'%0.1f')]};
                if length(parameters) > 3
                    STR(5) = {['\alpha_{diode} = ' num2str(parametersY(4),'%0.1f') ' \pm ' num2str(sigma_parY(4),'%0.1f')]};
                    STR(6) = {['\chi^2 per degree of freedom = ' num2str(chi2Y,'%8.2f')]};
                    STR(7) = {['Backing  = ' num2str(100*bacY,'%8.0f') ' %']};
                    STR2(2) = {['Cov(f_c,D)               =  ',num2str(covY(1)/sqrt(parameters(1)*parametersY(2)),'%0.3f')]};
                    STR2(3) = {['Cov(f_c,f_{diode})         =  ' num2str(covY(2)/sqrt(parametersY(1)*parametersY(3)),'%0.3f')]};
                    STR2(4) = {['Cov(f_c,\alpha)                =  ' num2str(covY(4)/sqrt(parametersY(1)*parametersY(4)),'%0.3f')]};
                    STR2(5) = {['Cov(D,f_{diode})          =  ' num2str(covY(3)/sqrt(parametersY(2)*parametersY(3)),'%0.3f')]};
                    STR2(6) = {['Cov(D,\alpha)                =  ' num2str(covY(5)/sqrt(parametersY(2)*parametersY(4)),'%0.3f')]};
                    STR2(7) = {['Cov(f_{diode},\alpha)          =  ' num2str(covY(6)/sqrt(parametersY(3)*parametersY(4)),'%0.3f')]};
                else
                    STR(5) = {['chi^2 per degree of freedom = ' num2str(chi2Y,'%0.2f')]};
                    STR(6) = {['Backing  = ' num2str(100*bacY,'%8.0f') ' %']};
                    STR2(2) = {['Cov(f_c,D)             = ',num2str(covY(1)/sqrt(parametersY(1)*parametersY(2)),'%0.3f')]};
                    STR2(3) = {['Cov(f_c,f_{diode})     = ' num2str(covY(2)/sqrt(parametersY(1)*parametersY(3)),'%0.3f')]};
                    STR2(4) = {['Cov(D,f_{diode})       = ' num2str(covY(3)/sqrt(parametersY(2)*parametersY(3)),'%0.3f')]};
                end;
            else
                STR(4) = {['\chi^2 per degree of freedom = ' num2str(chi2Y,'%8.2f')]};
                STR(5) = {['Backing = ' num2str(100*bacY,'%8.0f') ' %']};
            end;
            caption(STR,0.1,0.15);
            caption(STR2,0.65,0.15);
        end;
        
    end; %if ~exist('Px')
    
    check_alias;    %Checks if the precision of the fit could become better by including more aliasing terms
    
end

%%Update calibration values in stored file.
if prgVr.setupStt == 1
    if exist('DualTrapCalValues.txt', 'file') == 2
        %fid = fopen( 'DualTrapCalValues.txt', 'r' );
        writeValAppend = dlmread('DualTrapCalValues.txt');
        if sum(sum(eq(writeValAppend(:,1),writeVal),2)) == 1
            writeValAppend(find(writeValAppend(:,1) == writeVal),:) = writeVal;
            dlmwrite('DualTrapCalValues.txt',writeValAppend,'precision',16);
        else
            writeVal = [writeVal; writeValAppend];
            writeVal = sortrows(writeVal,1);
            dlmwrite('DualTrapCalValues.txt',writeVal,'precision',16);
        end
        %hisRead = textscan(fid, '%s %s %s %s %s %s %s\n', 'delimiter', '|','collectoutput',true);
       % fclose(fid);
        
    else        
        dlmwrite('DualTrapCalValues.txt',writeVal,'precision',16);
    end
    
elseif prgVr.setupStt == 2
    if exist('CTRAPCalValues.txt', 'file') == 2
        %fid = fopen( 'DualTrapCalValues.txt', 'r' );
        writeValAppend = dlmread('CTRAPCalValues.txt');
        if sum(sum(eq(writeValAppend(:,1),writeVal),2)) == 1
            writeValAppend(find(writeValAppend(:,1) == writeVal),:) = writeVal;
            dlmwrite('CTRAPCalValues.txt',writeValAppend,'precision',16);
        else
            writeVal = [writeVal; writeValAppend];
            writeVal = sortrows(writeVal,1);
            dlmwrite('CTRAPCalValues.txt',writeVal,'precision',16);
        end
        %hisRead = textscan(fid, '%s %s %s %s %s %s %s\n', 'delimiter', '|','collectoutput',true);
       % fclose(fid);
        
    else        
        dlmwrite('CTRAPCalValues.txt',writeVal,'precision',16);
    end
    
end
