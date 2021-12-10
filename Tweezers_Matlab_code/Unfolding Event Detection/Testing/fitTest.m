clear all;
close all;
clc;

%eWLCExpr = {'(kB*T)/(4*(1-x/Ld+F/K)^2)','(kB*T)/(4*pD)','(kB*T*x)/(Ld*pD)','(kB*T*F)/(K*pD)'}; %Extensible WLC Rief et al.
%wlcFit = fittype(@(kB,T,Ld,pD,K,x,F) (kB.*T)./(4.*(1-x./Ld+F./K).^2)-(kB.*T)./(4.*pD)+(kB.*T.*x)./(Ld.*pD)-(kB.*T.*F)./(K.*pD), 'independent', {'x', 'F'});

 load wormlikechains % wormlike chains
    plot( ext,  singleMBPForce,'r', ext, DNAforce, 'k')
   axis([0.5 1.5 0 85]);
kB = 1.3806488*10^(-23); %J/K
T = 293.15; %K
%Ld = 920E-9; %nm
pD = 50E-9; %nm
K = 1.2*10^(-9);
%wlceFit = @(kB,T,Ld,pD,K,x,F) (kB.*T)./(4.*(1-x./Ld+F./K).^2)-(kB.*T)./(4.*pD)+(kB.*T.*x)./(Ld.*pD)-(kB.*T.*F)./(K.*pD);
 
fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.

extension = ext(1:110).*1E-6;
force = DNAforce(1:110).*1E-12;

expFitVal = fminsearch(@(x) fitWLC( extension, force, kB, T, x(1), x(2), K ), [1000E-9 50E-9] , fMinOptions);

%plot(data(timeSpan,1),(data(timeSpan,protCol)-data(timeSpan(1),protCol)).*protConvFac,data(timeSpan,1),expFitVal(1)*(1-exp(-expFitVal(2).*(data(timeSpan,1)-data(timeSpan(1),1)))));

figure()
plot(extension*1E6,force*1E12,'linewidth',5);
 hold on;
 FEE = (kB.*T)./(4.*expFitVal(2).*(1-extension./expFitVal(1)+force./K).^2)-(kB.*T)./(4.*expFitVal(2))+(kB.*T.*extension)./(expFitVal(1).*expFitVal(2))-(kB.*T.*force)./(K.*expFitVal(2));
 plot(extension*1E6,FEE*1E12,'linewidth',3);
 axis([0.5 1.5 0 85]);
 %   wlceFit = @(kB,T,Ld,pD,x,F) (kB.*T)./(4.*(1-x./Ld).^2)-(kB.*T)./(4.*pD)+(kB.*T.*x)./(Ld.*pD);
% WLCFITTING = fit(ext, DNAforce, wlceFit)