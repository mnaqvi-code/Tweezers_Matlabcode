clear all;
close all;
clc;

% load wormlikechains % wormlike chains
% 
% KB = 1.3806488*10^(-23); %J/K
% kB = 1.3806488*10^(-23); %J/K
% T = 293.15; %K
% Ld = 920E-9; %nm
% Lp = 120E-9; %nm
% pP = 1.5E-9; %nm
% pD = 45E-9; %nm
% K = 1.2E-9;
% KB = 1.3806503*10^(-23);
% P = pD;
% L = Ld;
% xDNA = ext*1E-6;%[0.01:0.01:1.5*Ld*1E6].*(1E-6);
% 
% dnaCalcforce = calcEWLC( xDNA, KB, T, P, L, K  );
% K = 1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
% P = 1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
% L =  120*10^(-9);    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa 
%  
% 
% fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
% xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(dnaCalcforce)-1):1.5*Lp*1E6].*(1E-6);
% for ii = 1:length(dnaCalcforce)
%     xprot(ii) = fminsearch(@(x) fitTest2( dnaCalcforce(ii), x(1), KB, T, P, L, K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
% end
% 
% LpCalc = fminsearch(@(x) fitTest2( dnaCalcforce, xprot, KB, T, P, x(1), K ), 120E-9, fMinOptions); 
% LpCalc2 = fminsearch(@(x) calcProtLp( singleMBPForce, xprot, KB, T, P, x(1)), 120E-9, fMinOptions); 
% 
% 
% figure();
% plot(ext,DNAforce,'linewidth',5)
% hold on;
% plot((xDNA)*1E6,dnaCalcforce.*1E12,'linewidth',3)
% 
% plot(ext,singleMBPForce,'linewidth',5)
% 
% plot((xDNA+xprot)*1E6,dnaCalcforce.*1E12,'linewidth',3)
% hold off;
% axis([0 1.5 0 85]);

load WLCfittingtest2.mat
fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
KB = 1.3806488*10^(-23); %J/K
kB = 1.3806488*10^(-23); %J/K
T = 293.15; %K
Ld = 920E-9; %nm
Lp = 120E-9; %nm
pP = 1.5E-9; %nm
pD = 45E-9; %nm
K = 1.2E-9;
KB = 1.3806503*10^(-23);
P = pD;
L = Ld;

figure;
hold on;
plot(data.scombined{1}+data.shiftX,data.force{1}+data.shiftY)
plot(data.scombined{1}(75:140)+data.shiftX,data.force{1}(75:140)+data.shiftY,'r')
%plot(ext,DNAforce,'linewidth',2)
axis([0.5 1.1 0 35]);
Forcetest = data.force{1}(75:140)+data.shiftY;

for ii = 1:length(Forcetest)
    xDNA2(ii) = fminsearch(@(x) fitTest2( Forcetest(ii)*1E-12, x(1), KB, T, P, L, K ), 1E-6, fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end

plot(xDNA2*1E6,Forcetest,'linewidth',2)
xprot = (data.scombined{1}(75:140)+data.shiftX)*1E-6-xDNA2';
xprot = xprot;
%xDNA2 = data.scombined{1}(75:140)+data.shiftX;
dnaCalcforce2 = calcEWLC( xDNA2, KB, T, P, L, K  );
K = 1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
P = 1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
L =  120*10^(-9);    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa 
LpCalc3 = fminsearch(@(x) calcProtLp( Forcetest*1E-12, xprot, KB, T, P, x(1)), 120E-9, fMinOptions);
%plot((xprot'+xDNA2)*1E6,dnaCalcforce2*1E12,'linewidth',3)

fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(dnaCalcforce2)-1):1.5*Lp*1E6].*(1E-6);
for ii = 1:length(dnaCalcforce2)
    xprot(ii) = fminsearch(@(x) fitTest2( dnaCalcforce2(ii), x(1), KB, T, P, LpCalc3, K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end

plot((xDNA2+xprot')*1E6,dnaCalcforce2.*1E12,'linewidth',3)
