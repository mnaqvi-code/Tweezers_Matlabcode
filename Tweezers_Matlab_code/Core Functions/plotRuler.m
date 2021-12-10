function [ ] = plotRuler( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global datVr pltVr;
if isempty(datVr.CustRul);

    if pltVr.rulVal == 1
        
    elseif pltVr.rulVal == 2
        
    load Luciferasewlc.mat %Luciferase ruler to get the unfolding length    
    kB = 1.3806488*10^(-23); %J/K
    T = 293.15; %K
    deltaF = 0.1;
    fVal = [0.001:deltaF:100];
    bDsMax = 920/1000; %nm/bp
    pDs = 50; %nm
    sDs = 1200; %pN
    bDs = bDsMax.*(1 - (1/2).*sqrt((kB*T*1E21)./(fVal.*pDs)) + fVal./sDs);
    plot(bDs,fVal,'b','linewidth',2);
    
    pSs = 0.75; %nm
    bSsMax = 2687*0.623/1000; %nm
    sSs = 800; %pN
    bSs = bSsMax.*(coth(2.*(fVal.*pSs)/(kB*T*1E21)) - 1./((2.*fVal.*pSs)/(kB*T*1E21))).*(1 + fVal./sSs); 
    %Fit constants
    
    a = 0.34E-9; %nm
    L = 920E-9; %2687*a; %um
    
    fRul = [0.001:deltaF:100];

    fRul = fRul*1E-12;
    F = a/(kB*T)*fRul;
    kapB = 147;
    kapS = 4;
    kapBS = kapS;

    gam = 1.829;
    mu = 4.1;
    J = 2.15;
    Eb = 1200*1E-12;
    Et = Eb/(kB*T).*a;
    L0 = 920E-9;%2687*a;

    alphaB = sqrt(kapB.*F+(F./2).^2);
    alphaS = sqrt(kapS.*F+(F./2).^2);
    mu0 = mu - log(gam) + F.*((1-gam)./2)+0.5.*log((kapS+F.*gam./2+alphaS)./(kapB+F./2+alphaB));
    J0 = J + 0.25.*log(((kapBS+F./2+alphaB).*(kapBS+F.*gam./2+alphaS))./((kapB+F./2+alphaB).*(kapS+F.*gam./2+alphaS)));
    sigNorm = sinh(mu0)./sqrt(sinh(mu0).^2+exp(-4.*J0));
    sigN1 = sigNorm.^2+(1-sigNorm.^2).*(cosh(mu0)-sqrt(sinh(mu0).^2+exp(-4.*J0)))./(cosh(mu0)+sqrt(sinh(mu0).^2+exp(-4.*J0)));
    phiB = (1+sigNorm)./2;
    phiS = (1-sigNorm)./2;
    zNorm = (1+F./Et-1./(2.*alphaB)).*phiB+gam.*(1-1./(2.*alphaS)).*phiS+((sigN1-1)./4).*(1./(2.*alphaB).*(kapB-kapBS)./(kapBS+F./2+alphaB)+gam./(2.*alphaS).*(kapS-kapBS)./(kapBS+gam.*F./2+alphaS));
     plot(zNorm*L0*1E6,fRul*1E12,'k','linewidth',1);
    
    kB = 1.3806488*10^(-23); %J/K
    T = 293.15; %K
    a = 0.665E-9; %nm
    L = 2553*a; %um
    z = 0.01:0.0001:1.6;
    z = z*1E-6;
    kap = 1.28;
    u = coth(kap) - 1/kap;
    F = z./L.*(3*(1-u)/(1+u)-1/sqrt(1+4*kap^2))+sqrt(1./(1-z./L).^2+4*kap^2)-sqrt(1+4*kap^2);
    fRul = F*kB*T/a;
    Unl = 1.172777.*fRul.*1E8-3.731836.*(fRul.*1E8).^2+4.118249.*(fRul.*1E8).^3;
    F = z.*(1+Unl)./L.*(3*(1-u)/(1+u)-1/sqrt(1+4*kap^2))+sqrt(1./(1-z.*(1+Unl)./L).^2+4*kap^2)-sqrt(1+4*kap^2);
    fRul = F*kB*T/a;
    
    plot(z*1E6,fRul*1E12,'g','linewidth',2);  
        
    elseif pltVr.rulVal == 3
        
    load sMBPrulerMinimal.mat %MBP ruler to get the unfolding length
    plot( rulerext,  rulerf, 'k')  
    
    load wormlikechains % wormlike chains
    plot( ext,  singleMBPForce,'r', ext, DNAforce, 'k')
       
    elseif pltVr.rulVal == 4
      
    load 4MBPRULErdual.mat
    plot(extension(:,1),force(:,1)*1E12,'k');
    plot(extension(:,2),force(:,1)*1E12,'k');
    plot(extension(:,3),force(:,1)*1E12,'k');
    plot(extension(:,4),force(:,1)*1E12,'k');
    plot(extension(:,5),force(:,1)*1E12,'k');
    plot(extension(:,6),force(:,1)*1E12,'k');
        
%     load 4MBPRULEr.mat    
%     plot(MBP4(1:151,1),MBP4(1:151,2),'k');
%     plot(MBP4(152:302,1),MBP4(152:302,2),'k');
%     plot(MBP4(303:453,1),MBP4(303:453,2),'k');
%     plot(MBP4(454:604,1),MBP4(454:604,2),'k');
%     plot(MBP4(605:755,1),MBP4(605:755,2),'k');
%     plot(MBP4(756:906,1),MBP4(756:906,2),'k');
    
    elseif pltVr.rulVal == 5
        
    load Luciferasewlc.mat %Luciferase ruler to get the unfolding length
    plot( intermediates(:,1), intermediates(:,2), 'k', ext, OneLuci,'r', ext, f0, 'k');
        
    elseif pltVr.rulVal == 6
        
    load 1MBPRULERdualhandles.mat
    %plot(extension(:,1),force(:,1),'k');
    %plot(extension(:,2),force(:,1),'k');
    %plot(extension(:,3),force(:,1),'k');    
    plot(extension(:,1)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,2)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,3)*1E6,force(:,1)*1E12,'k');
    
    elseif pltVr.rulVal == 7
    load aSynRuler.mat
    plot(x*1E6,force*1E12, 'k');
    plot((x+xprot1)*1E6,force*1E12, 'k');
    plot((x+xprot2)*1E6,force*1E12, 'k');
    
    elseif pltVr.rulVal == 8
    load NCS1RULER.mat
    plot(extension(:,1),force(:,1)*1E12,'k');
    plot(extension(:,2),force(:,1)*1E12,'k');
    plot(extension(:,3),force(:,1)*1E12,'k');
    plot(extension(:,4),force(:,1)*1E12,'k');
    
    elseif pltVr.rulVal == 9
    load NCS1RULERMISF.mat
    plot(extension(:,1),force(:,1)*1E12,'k');
    plot(extension(:,2),force(:,1)*1E12,'k');
    plot(extension(:,3),force(:,1)*1E12,'k');
    plot(extension(:,4),force(:,1)*1E12,'k');
    plot(extension(:,5),force(:,1)*1E12,'k');
    plot(extension(:,6),force(:,1)*1E12,'k');
    
    elseif pltVr.rulVal == 10
    load GRshortlong.mat
    plot(extension(:,1)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,2)*1E6,force(:,1)*1E12,'k');
        
    elseif pltVr.rulVal == 11
    load RHODRULER.mat
    plot(extension(:,1)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,2)*1E6,force(:,1)*1E12,'k');    
    
    elseif pltVr.rulVal == 12
    load 4SynRuler.mat
    plot(extension(:,1)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,2)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,3)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,4)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,5)*1E6,force(:,1)*1E12,'k');
    
    elseif pltVr.rulVal == 13
    load 1MBPRULERdualhandleslong.mat
    plot(extension(:,1)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,2)*1E6,force(:,1)*1E12,'k');
    plot(extension(:,3)*1E6,force(:,1)*1E12,'k');
    
    end
else
    for nn = 3:size(datVr.CustRul,2)-1
        plot(datVr.CustRul(:,nn)*1E6,datVr.CustRul(:,2)*1E12,'k');
    end
    
    
end

end

