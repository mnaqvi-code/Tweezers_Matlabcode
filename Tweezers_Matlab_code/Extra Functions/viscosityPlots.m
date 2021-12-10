clear all;
close all;
clc;

viscC = [0.1256735; 1.265347; -1.105369; 0.2044679; 1.308779];
alphaC = [-0.04296718; 0.3710073; 0.4230889; -0.3259828];

T = 0:0.1:50;
m = [0; 100E-3; 250E-3; 1];

for ii = 1:length(m)
    eta{ii} = viscC(1) +  viscC(2)*exp(alphaC(1).*T) + viscC(3)*exp(alphaC(2).*m(ii)) + viscC(4)*exp(alphaC(3).*(0.01.*T + m(ii))) + viscC(5)*exp(alphaC(4).*(0.01.*T - m(ii)));
end

figure;
hold on;
for ii = 1:length(eta)
    plot(T,eta{ii});
end
hold off;
axis([0 50 0.5 2]);
legend('Water', '100mM NaCl', '250mM NaCl', '1M NaCl')