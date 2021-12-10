function [ PDNA ] = calcPDNA( cSalt, zVal )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%From Manning et al. (2006) Eq 25 implemented.

bDNA = 0.17E-9;
RDNA = 1E-9;
pStarDNA = 7.4E-9;
dzetaDNA = 4.2;
lB = dzetaDNA*bDNA;
kap = sqrt(cSalt.*8.*pi.*lB*1E9);
kap = kap.*1E9;
PDNA = (pi./2).^(2./3).*RDNA.^(4./3).*pStarDNA.^(2./3).*zVal.^(-2).*lB.^(-1).*((2.*zVal.*dzetaDNA-1).*(kap.*bDNA.*exp(-kap.*bDNA))./(1-exp(-kap.*bDNA))-1-log(1-exp(-kap.*bDNA)));
PDNA = PDNA.*1E9;

end

