function [ error ] = fitWLCMBP( extension, pExtension, DNAforce, force, kB, T, Ld, Lp, pD, pP, K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

error = sum(abs((kB.*T)./(4.*pD.*(1-extension./Ld+DNAforce./K).^2)-(kB.*T)./(4.*pD)+(kB.*T.*extension)./(Ld.*pD)-(kB.*T.*DNAforce)./(K.*pD)+(kB.*T)./(4.*pP.*(1-(pExtension)./Lp).^2)-(kB.*T)./(4.*pP)+(kB.*T.*(pExtension))./(Lp.*pP)-(force)));
