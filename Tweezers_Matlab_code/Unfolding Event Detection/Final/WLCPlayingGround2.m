clear all;
close all;
clc;

load WLCfittingtest2.mat
fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
kB = 1.38064852*10^(-23); %J/K
T = 293.15; %K
Ld = 920E-9; %nm
Lp = 120E-9;%370*0.%120E-9; %nm
pP = 1.5E-9; %nm
pD = 45E-9; %nm
K = 1.2E-9;
P = pD;
L = Ld;

datI = 75:140;

figure;
hold on;
plot(data.scombined{1}+data.shiftX,data.force{1}+data.shiftY)
plot(data.scombined{1}(datI)+data.shiftX,data.force{1}(datI)+data.shiftY,'r')

datExt = data.scombined{1}(datI)+data.shiftX;
datForce = data.force{1}(datI)+data.shiftY;


[sdatForce, i]  = sort(datForce);
sdatExt = datExt(i);

%plot(ext,DNAforce,'linewidth',2)
axis([0.5 1.1 0 35]);
Forcetest = data.force{1}(datI)+data.shiftY;
for ii = 1:length(datForce)
    xDNA2(ii) = fminsearch(@(x) fitEWLC( Forcetest(ii)*1E-12, x(1), kB, T, P, L, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end
% for ii = 1:length(datForce)
%     xDNA2(ii) = fminsearch(@(x) fitEWLC( datForce(ii)*1E-12, x(1), kB, T, P, L, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))
% end

plot(xDNA2*1E6,Forcetest,'linewidth',2)
xprot = (data.scombined{1}(datI)+data.shiftX)*1E-6-xDNA2';
xprot = xprot;
dnaCalcforce2 = calcEWLC( xDNA2, kB, T, P, L, K  );
K = 1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
P = 1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
L =  120*10^(-9);    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa 
LpCalc3 = fminsearch(@(x) calcProtLp( Forcetest*1E-12, xprot, kB, T, P, x(1)), 120E-9, fMinOptions);

fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(dnaCalcforce2)-1):1.5*Lp*1E6].*(1E-6);
for ii = 1:length(dnaCalcforce2)
    xprot(ii) = fminsearch(@(x) fitEWLC( dnaCalcforce2(ii), x(1), kB, T, P, LpCalc3, K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end

plot((xDNA2+xprot')*1E6,dnaCalcforce2.*1E12,'linewidth',3)