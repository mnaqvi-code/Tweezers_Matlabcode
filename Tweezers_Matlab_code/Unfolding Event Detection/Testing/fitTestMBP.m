clear all;
close all;
clc;

kB = 1.3806488*10^(-23); %J/K
T = 293.15; %K
Ld = 920E-9; %nm
Lp = 120E-9; %nm
pP = 1.5E-9; %nm
pD = 45E-9; %nm
K = 1.2E-9;

load wormlikechains % wormlike chains
ext = ext(1:end);
DNAforce = DNAforce(1:end);
force = DNAforce.*1E-12;
singleMBPForce = singleMBPForce(1:end);
extension = ext.*1E-6;
DNAforce = (kB.*T)./(4.*pD.*(1-extension./Ld+force./K).^2)-(kB.*T)./(4.*pD)+(kB.*T.*extension)./(Ld.*pD)-(kB.*T.*force)./(K.*pD);
protForce = DNAforce-singleMBPForce*1E-12;
figure();
hold on;
plot(extension*1E6, protForce*1E12);
plot( ext,  singleMBPForce,'r', ext, DNAforce*1E12, 'k')
axis([0.5 1.5 0 85]);
hold off;
fMinOptions = optimset('MaxFunEvals', 10000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
expFitVal = fminsearch(@(x) fitWLCMBP( extension, x(1), DNAforce, singleMBPForce*1E-12, kB, T, Ld, Lp, pD, pP, K ), linspace(0, 0.9*Lp, length(extension)), fMinOptions); %linspace(0, 0.9*Lp, length(extension))

%plot(data(timeSpan,1),(data(timeSpan,protCol)-data(timeSpan(1),protCol)).*protConvFac,data(timeSpan,1),expFitVal(1)*(1-exp(-expFitVal(2).*(data(timeSpan,1)-data(timeSpan(1),1)))));

figure();
plot(extension*1E6,singleMBPForce,'linewidth',5)
hold on;
plot((extension+120E-9)*1E6, ((kB.*T)./(4.*pD.*(1-extension./Ld+DNAforce./K).^2)-(kB.*T)./(4.*pD)+(kB.*T.*extension)./(Ld.*pD)-(kB.*T.*DNAforce)./(K.*pD)+(kB.*T)./(4.*pP.*(1-(expFitVal)./Lp).^2)-(kB.*T)./(4.*pP)+(kB.*T.*(expFitVal))./(Lp.*pP)).*1E12);
hold off;
axis([0 1.7 0 85]);
figure;
FEE = (kB.*T)./(4.*pP.*(1-1020E-9./Lp).^2)-(kB.*T)./(4.*pP)+(kB.*T.*1020E-9)./(Lp.*pP);
plot(extension*1E6,FEE*1E12,'linewidth',3);
axis([0.5 1.5 0 85]);

figure;
pExtension = linspace(0, 0.9*Lp, length(extension))
plot((pExtension+Ld)*1E6,((kB.*T)./(4.*pP.*(1-pExtension./Lp).^2)-(kB.*T)./(4.*pP)+(kB.*T.*pExtension)./(Lp.*pP))*1E12, extension*1E6, DNAforce*1E12,extension*1E6, singleMBPForce)
axis([0 1.7 0 85]);
