function [ error, bDs ] = fitWLC( extension, force, kB, T, Ld, pD, K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%@(kB,T,Ld,pD,x,F) (kB.*T)./(4.*(1-extension./Ld).^2)-(kB.*T)./(4.*pD)+(kB.*T.*extension)./(Ld.*pD);

%keyboard;
error = sum(abs(((kB.*T)./(4.*pD.*(1-extension./Ld+force./K).^2)-(kB.*T)./(4.*pD)+(kB.*T.*extension)./(Ld.*pD)-(kB.*T.*force)./(K.*pD)-(force))));
