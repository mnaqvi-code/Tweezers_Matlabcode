clear all;
close all;
clc;

load wormlikechains % wormlike chains
plot( ext,  singleMBPForce,'r', ext, DNAforce, 'k', ext,  DNAforce-singleMBPForce,'g')
axis([0.5 1.5 0 85]);
ext = ext;
DNAforce = DNAforce;
singleMBPForce = singleMBPForce;
extension = ext.*1E-6;
force = singleMBPForce.*1E-12;

kB = 1.3806488*10^(-23); %J/K
T = 293.15; %K
Ld = 920E-9; %nm
Lp = 120E-9; %nm
pP = 1.5E-9; %nm
pD = 45E-9; %nm
K = 1.2E-9;

%force = (DNAforce-singleMBPForce)*1E-12;


fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
expFitVal = fminsearch(@(x) fitMBPonly( extension, x(1), force, kB, T, Ld, Lp, pD, pP, K ), linspace(0, 120E-9, length(extension)) , fMinOptions);

%plot(data(timeSpan,1),(data(timeSpan,protCol)-data(timeSpan(1),protCol)).*protConvFac,data(timeSpan,1),expFitVal(1)*(1-exp(-expFitVal(2).*(data(timeSpan,1)-data(timeSpan(1),1)))));

figure()
plot(extension*1E6,force*1E12,'linewidth',5)
hold on;


FEE = (kB.*T)./(4.*pP.*(1-expFitVal./Lp).^2)-(kB.*T)./(4.*pP)+(kB.*T.*expFitVal)./(Lp.*pP);
plot(extension*1E6,FEE*1E12,'linewidth',3);
axis([0.5 1.5 0 85]);
%   wlceFit = @(kB,T,Ld,pD,x,F) (kB.*T)./(4.*(1-x./Ld).^2)-(kB.*T)./(4.*pD)+(kB.*T.*x)./(Ld.*pD);
% WLCFITTING = fit(ext, DNAforce, wlceFit)