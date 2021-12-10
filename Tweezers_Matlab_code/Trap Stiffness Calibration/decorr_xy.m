global pltVr



%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FORM Pxy
delta_t = 1./(2.*fNyq);

if nCal == 0
    clear PxSs PySs Pxy
    for ii = 1:length(Y{jj})/(1000*sampling_f)
        PxSs{ii} = delta_t*fft(X{jj}(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)));
        PySs{ii} = delta_t*fft(Y{jj}(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)));
        %PxyS{ii} = real(delta_t*fft(X{jj}(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii))) .* conj(delta_t*fft(Y{jj}(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii))))) / T;
        %PxyS{ii} = PxyS{ii}(find(f <= fNyq));
    end
    Pxy = real(mean([PxSs{:}],2).*conj(mean([PySs{:}],2))) / T;
    Pxy = Pxy(find(f <= fNyq));
    %Pxy = mean([PxyS{:}],2);
else
    Pxy = real(delta_t*fft(X{jj}) .* conj(delta_t*fft(Y{jj}))) / T;
    Pxy = Pxy(find(f <= fNyq));
end

% %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% % FORM Pxy
% delta_t = 1./(2.*fNyq);
%
% Pxy = real(delta_t*fft(X{jj}) .* conj(delta_t*fft(Y{jj}))) / T;
% Pxy = Pxy(find(f <= fNyq));
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%  Choose data to be fitted
ind     = find(f > 50 & f <= 20000);
fd      = f(ind);
Pxd     = Px{jj}(ind);
Pyd     = Py{jj}(ind);
Pxyd    = Pxy(ind);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% BIN Px, Pxy, Pxy
nbin    =   floor(length(fd)/nblock);
fb = []; Pxb = []; Pyb = []; Pxyb = [];
for i = 1 : nbin
    fb(i)   = mean_nan(fd((i-1)*nblock+1 : i*nblock));
    Pxb(i)  = mean_nan(Pxd((i-1)*nblock+1 : i*nblock));
    Pyb(i)  = mean_nan(Pyd((i-1)*nblock+1 : i*nblock));
    Pxyb(i) = mean_nan(Pxyd((i-1)*nblock+1 : i*nblock));
end
fb = fb(:); Pxb = Pxb(:); Pyb = Pyb(:); Pxyb = Pxyb(:);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FIT b,c

[parameters,RESNORM,RESIDUAL,JACOBIAN]=fit_nonl(@min_corr,[0.1 0.1],tolx,50,Pxb,Pyb,Pxyb);

% [parameters,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@min_corr,[0 0],[],[],optimset('TolX',tolx,'Display','iter','MaxIter',100,'MaxFunEvals',1000,'LargeScale','off'),Pxb,Pyb,Pxyb);

b = parameters(1);
c = parameters(2);
M =  ((1 + b*c) .* Pxyb + c*Pxb + b*Pyb) ./ sqrt((Pxb + 2*b*Pxyb + b^2*Pyb) .* (Pyb + 2*c*Pxyb + c^2*Pxb));

disp(['b =  ' num2str(b,'%5.2f') ',   c =  ' num2str(c,'%5.2f')]);

Px1     = Px{jj} + 2*b*Pxy + b^2*Py{jj};
Py1     = Py{jj} + 2*c*Pxy + c^2*Px{jj};
Px1y1   = (1 + b*c) .* Pxy + c*Px{jj} + b*Py{jj};
for i = 1 : nbin
    Px1b(i)   = mean_nan(Px1((i-1)*nblock+1 : i*nblock));
    Py1b(i)   = mean_nan(Py1((i-1)*nblock+1 : i*nblock));
    Px1y1b(i) = mean_nan(Px1y1((i-1)*nblock+1 : i*nblock));
end
Px1b = Px1b(:); Py1b = Py1b(:); Px1y1b = Px1y1b(:);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% PLOT

if jj == 1
    pltVr.figs.PSDNorm = figure;
    set(gcf,'Numbertitle','off','Name','PSD With No Crosstalk');
    hf2 = suptitle('Power spectra, after cross-talk has been eliminated');
    set(hf2,'FontSize',20,'FontWeight','Bold')
elseif jj > 1
    set(0, 'currentfigure', pltVr.figs.PSDNorm);
end

hold on;
set(0, 'currentfigure', pltVr.figs.PSDNorm); shg; subplot(1,cIter,jj); set(gcf,'Numbertitle','off','Name','PSD With No Crosstalk');
set(gca,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
%hf7 = title('Power spectra, after cross-talk has been eliminated'); set(hf7,'Fontweight','Bold');
plot_powerspectrum(f,Px1,color_x,2);
plot_powerspectrum(f,Py1,color_y,2);
plothx = findobj('color',color_x);
plothy = findobj('color',color_y);
legend([plothx(1) plothy(1)],{'X''','Y'''});
hold off;

if jj == 1
    pltVr.figs.elCross = figure;
    set(gcf,'Numbertitle','off','Name','Elimination of crosstalk');
    hf2 = suptitle('Eliminating cross-talk between channels');
    set(hf2,'FontSize',20,'FontWeight','Bold')
elseif jj > 1
    set(0, 'currentfigure', pltVr.figs.elCross);
end

set(0, 'currentfigure', pltVr.figs.elCross); shg; subplot(cIter,1,jj); set(gcf,'Numbertitle','off','Name','Elimination of crosstalk');
set(gca,'FontUnits','normalized','FontSize',0.04,'FontWeight','bold');
%hf3 = title('Eliminating cross-talk between channels'); set(hf3,'Fontweight','Bold');
hold on;
plot_corrxy(fb,Pxb,Pyb,Pxyb,'k',1);
plot_corrxy(fb,Px1b,Py1b,Px1y1b,'r',2);
plot_x = findobj('marker','o');
plot_y = findobj('marker','s');
legend([plot_x(1) plot_y(1)],'Cross-talk not eliminated','Cross-talk eliminated');
STR ={};
STR(1) = {['b = ' num2str(b,'%5.2f') ', c = ' num2str(c,'%5.2f')]};
caption(STR,0.1,0.1+0.5*(cIter-1)*(jj-1));
drawnow;
pause(0.1);
jFrame = get(handle(gcf),'JavaFrame');
jFrame.setMaximized(true);
hold off;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Px{jj} = Px1;
Py{jj} = Py1;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Also transform data for position:

X1{jj} =  X{jj} + b*Y{jj};
Y1{jj} =  Y{jj} + c*X{jj};
X{jj} = X1{jj};
Y{jj} = Y1{jj};

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
