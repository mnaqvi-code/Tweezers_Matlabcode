global pltVr

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% CALCULATE power spectra
check1(nblock,color1);
check(sampling_f,0,1000,'sampling frequency',color1);
check2(channel(1),1,size(File,2),'Channel X-1',color1);
check2(channel(2),0,size(File,2),'Channel Y-1',color1);
if prgVr.setupStt == 2
    if channel(3) ~= 0 && channel(4) ~= 0
        check2(channel(3),0,size(File,2),'Channel X-2',color1);
        check2(channel(4),0,size(File,2),'Channel Y-2',color1);
        cIter = 2;
    else
        cIter = 1;
    end
    
else
    cIter = 1;
end

% if channel_z ~= 0,
%   X       = File(:,channel_x)./File(:,channel_z);
%   Y       = File(:,channel_y)./File(:,channel_z);
% else
%   X       = File(:,channel_x);
%   Y       = File(:,channel_y);
% end;

%(if channel_z ~= 0)
% X       = X - mean(X);
% [f,Px]  = calc_powersp(X,1000*sampling_f);

for jj = 1:cIter
    
    X = File(:,channel(jj*2-1));
    Y = File(:,channel(jj*2));
    
    if nCal == 0
        tRemove = mod(length(X)*1/(1000*sampling_f),1/32);
        pRemove = tRemove * 1000*sampling_f;
        X = X(1:end-pRemove);
        X = X - mean(X);
        for ii = 1:length(X)/(1000*sampling_f)
            [fS{ii},PxS{ii},TS{ii}]  = calc_powersp(X(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)),1000*sampling_f);
        end
        f = mean([fS{:}],2);
        Px = mean([PxS{:}],2);
        T = mean([TS{:}],2);
    else
        X = X - mean(X);
        [f,Px,T]  = calc_powersp(X,1000*sampling_f);
    end
    if nCal == 0
        Y = Y(1:end-pRemove);
        Y = Y - mean(Y);
        for ii = 1:length(Y)/(1000*sampling_f)
            [fS{ii},PyS{ii},TS{ii}]  = calc_powersp(Y(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)),1000*sampling_f);
        end
        f = mean([fS{:}],2);
        Py = mean([PyS{:}],2);
        T = mean([TS{:}],2);
    else
        Y = Y - mean(Y);
        [f,Py,T]  = calc_powersp(Y,1000*sampling_f);
    end
    % Y       = Y - mean(Y);
    % [f,Py]  = calc_powersp(Y,1000*sampling_f);
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % FORM Pxy
    delta_t = 1./(2.*fNyq);
    if nCal == 0
        for ii = 1:length(Y)/(1000*sampling_f)
            PxSs{ii} = delta_t*fft(X(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)));
            PySs{ii} = delta_t*fft(Y(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii)));
            %         PxyS{ii} = real(delta_t*fft(X(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii))) .* conj(delta_t*fft(Y(1+(1000*sampling_f)*(ii-1):(1000*sampling_f)*(ii))))) / T;
            %         PxyS{ii} = PxyS{ii}(find(f <= fNyq));
        end
        Pxy = real(mean([PxSs{:}],2).*conj(mean([PySs{:}],2))) / T;
        Pxy = Pxy(find(f <= fNyq));
        %     Pxy = mean([PxyS{:}],2);
    else
        Pxy = real(delta_t*fft(X) .* conj(delta_t*fft(Y))) / T;
        Pxy = Pxy(find(f <= fNyq));
    end
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    %  Choose data to be fitted
    ind     = find(f > 50 & f <= fNyq);
    fd      = f(ind);
    Pxd     = Px(ind);
    Pyd     = Py(ind);
    Pxyd    = Pxy(ind);
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % BIN Px, Pxy, Pxy
    
    nbin    =   floor(length(fd)/nblock);
    fb = []; Pxb = [];Pyb = []; Pxyb = [];
    for i = 1 : nbin
        fb(i)   = mean_nan(fd((i-1)*nblock+1 : i*nblock));
        Pxb(i)  = mean_nan(Pxd((i-1)*nblock+1 : i*nblock));
        Pyb(i)  = mean_nan(Pyd((i-1)*nblock+1 : i*nblock));
        Pxyb(i) = mean_nan(Pxyd((i-1)*nblock+1 : i*nblock));
    end;
    fb = fb(:); Pxb = Pxb(:); Pyb = Pyb(:); Pxyb = Pxyb(:);
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % PLOT
    
    if isfield(pltVr.figs, 'decorr') && jj == 1
        if ~isempty(findobj('name','Crosstalk Check'))
            close(pltVr.figs.decorr)
        end
        rmfield(pltVr.figs, 'decorr');
    end
    
    if jj == 1
        pltVr.figs.decorr = figure;
        set(gcf,'Numbertitle','off','Name','Crosstalk Check');
        hf2 = title(['Checking for cross-talk between channels X-', num2str(jj), ' and Y-', num2str(jj)]);
        set(hf2,'Fontweight','Bold');
    elseif jj > 1
        set(0, 'currentfigure', pltVr.figs.decorr);
    end
    %pltVr.figs.decorr.(genvarname(['Plot', cIter])) = figure;
    
    subplot(1,cIter,jj); 
    hold on;      
    plot_corrxy(fb,Pxb,Pyb,Pxyb,'k',1);
    hold off;
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
end