clear all;
close all;
clc;

kB = 1.38064852*10^(-23); %J/K
T = 293.15; %K
Ld = 2666*0.34*1E-9;%700.1E-9;%920E-9; %nm %1333 basepairs for both dual DNA handles
Lp = 495E-9;%10E-9;%120E-9;%370*0.%120E-9; %nm
Lp = 1500*0.33E-9; %Newly Adjusted
pP = 1.5E-9; %nm
pD = 45E-9; %nm
K = 1.2E-9;
P = pD;
L = Ld;


%fourMBPL = ([112, 204, 296, 388, 480]+15).*10^(-9); 

x = 0.4:0.001:1.1;
x = x.*1E-6

[ force ] = calcEWLC( x, kB, T, P, L, K  )
K = 2.0E-9;%1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
P = 0.8E-9;%1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
L =  Lp;%140*0.36E-9;%120*10^(-9);    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa 
LpCalc3 = ([368, 276, 184, 92, 0]).*10^(-9); 

%Cterminal length is 94 AA, 5*4 AA linkers between. Nterminal length is 276
%AA. Total protein length is 1500 AA (1501 according to Eline), MBP protein
%length is 370 AA.

LpCalc3 = ([396 672 948 1224 1500]).*0.33E-9;
xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(force)-1):1.5*Lp*1E6].*(1E-6);
fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
for ii = 1:length(force)
    xprot1(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, LpCalc3(1), K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end
for ii = 1:length(force)
    xprot2(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, LpCalc3(2), K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end
for ii = 1:length(force)
    xprot3(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, LpCalc3(3), K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end
for ii = 1:length(force)
    xprot4(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, LpCalc3(4), K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end
for ii = 1:length(force)
    xprot5(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, LpCalc3(5), K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end



plot(x*1E6,force*1E12,(x+xprot1)*1E6,force*1E12,(x+xprot2)*1E6,force*1E12,(x+xprot3)*1E6,force*1E12,(x+xprot4)*1E6,force*1E12,(x+xprot5)*1E6,force*1E12);