function [ error ] = calcDNAF( force, x, KB, T, P, L, K)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    error = sum(abs(force - ((KB.*T./P).*(x./L+1./(4.*(1-x./L).^2)-1/4-force./K))));


end

