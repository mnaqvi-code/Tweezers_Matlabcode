clear all;
close all;
clc;

kB = 1.38064852*10^(-23); %J/K
T = 293.15; %K
Ld = 2666*0.338*1E-9;%2553*0.338*1E-9;%Philip-Ruud thesis%2687*0.34*1E-9;%Peter Thesis basepair length; %2666*0.34*1E-9;% Mario/Eline dual ruler %700.1E-9;%920E-9; %nm %1333 basepairs for both dual DNA handles
Lp = 376*0.365*1E-9;%376*0.33E-9;%120E-9;%10E-9;%120E-9;%370*0.%120E-9; %nm
pP = 1.5E-9; %nm
pD = 45E-9; %nm
K = 1.2E-9;
P = pD;
L = Ld;

x = 0.4:0.001:1.1;
x = x.*1E-6

[ force ] = calcEWLC( x, kB, T, P, L, K  )
K = 2.0E-9;%1.5*10^(-9);  % ref; PRL 94, 048301 (2005)
P = 0.8E-9;%1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
L =  Lp;%120E-9;%140*0.36E-9;%120*10^(-9);    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa 
LpCalc3 = 92E-9;
xprotGuess = [0.01:(1.5*Lp*1E6-0.01)/(length(force)-1):1.5*Lp*1E6].*(1E-6);
fMinOptions = optimset('MaxFunEvals', 2000, 'Display', 'Iter', 'TolFun', 1E-13, 'TolX', 1E-16); % fminsearch options.
for ii = 1:length(force)
    xprot1(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, L-LpCalc3, K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end
for ii = 1:length(force)
    xprot2(ii) = fminsearch(@(x) fitEWLC( force(ii), x(1), kB, T, P, L, K ), xprotGuess(ii), fMinOptions); %linspace(0, 0.9*Lp, length(extension))
end


plot(x*1E6,force*1E12,(x+xprot1)*1E6,force*1E12,(x+xprot2)*1E6,force*1E12);