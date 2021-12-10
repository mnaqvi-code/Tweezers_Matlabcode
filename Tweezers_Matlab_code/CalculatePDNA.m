clear all;
close all;
clc;

%From Manning et al. (2006) Eq 25 implemented.

bDNA = 0.17E-9; %m
%pStarDNA = 3.07*L^(0.083);
RDNA = 1E-9; %m
pStarDNA = 7.4E-9;%6.6E-9;%7.4E-9; %m
cSalt = -6:0.01:6; %M
cSalt = 10.^(cSalt); %M6
zVal = 1;
dzetaDNA = 4.2;%10.4;%*zVal;
lB = dzetaDNA*bDNA;%
%kap = sqrt(cSalt)./0.304;
kap = sqrt(cSalt.*8.*pi.*lB*1E9);
kap = kap.*1E9;
%kap = 1/0.96E-9;

PDNA1 = (pi./2).^(2./3).*RDNA.^(4./3).*pStarDNA.^(2./3).*zVal.^(-2).*lB.^(-1).*((2.*zVal.*dzetaDNA-1).*(kap.*bDNA.*exp(-kap.*bDNA))./(1-exp(-kap.*bDNA))-1-log(1-exp(-kap.*bDNA)));
PDNA1 = PDNA1.*1E9;

zVal = 2;
pStarDNA = 7.4E-9;%14E-9;%7.4E-9; %m
dzetaDNA = 4.2;%*zVal;
lB = dzetaDNA*bDNA;
%kap = sqrt(cSalt)./0.304;
kap = sqrt(cSalt.*8.*pi.*lB*1E9);
kap = kap.*1E9;
PDNA2 = (pi./2).^(2./3).*RDNA.^(4./3).*pStarDNA.^(2./3).*zVal.^(-2).*lB.^(-1).*((2.*zVal.*dzetaDNA-1).*(kap.*bDNA.*exp(-kap.*bDNA))./(1-exp(-kap.*bDNA))-1-log(1-exp(-kap.*bDNA)));
PDNA2 = PDNA2.*1E9;
figure;
%plot(1./cSalt,PDNA1)
plot(cSalt,PDNA1)
hold on;
%plot(1./cSalt,PDNA2)
plot(cSalt,30.1+(53.9-30.1)./(1+(cSalt./0.174).^0.994));
plot(cSalt,33.8+(51.7-33.8)./(1+(cSalt./0.045).^0.546));
plot(cSalt,PDNA2)
legend('ThZ1','ExZ1','ExZ2','ThZ2');
axis([0 1 0 60]);
%axis([10^(-3) 10^(3) 0 110])