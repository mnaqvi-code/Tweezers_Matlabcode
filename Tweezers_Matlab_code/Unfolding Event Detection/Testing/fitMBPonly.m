function [ error ] = fitMBPonly( extension, pExtension, force, kB, T, Ld, Lp, pD, pP, K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%@(kB,T,Ld,pD,x,F) (kB.*T)./(4.*(1-extension./Ld).^2)-(kB.*T)./(4.*pD)+(kB.*T.*extension)./(Ld.*pD);

%keyboard;
error = sum(abs((kB.*T)./(4.*pP.*(1-pExtension./Lp).^2)-(kB.*T)./(4.*pP)+(kB.*T.*pExtension)./(Lp.*pP)-(force)));
