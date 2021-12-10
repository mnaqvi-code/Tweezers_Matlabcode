function [ error ] = fitEWLC( dnaCalcforce, x, KB, T, P, L, K)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
error = sum(abs(dnaCalcforce - calcEWLC( x, KB, T, P, L, K)));

end

