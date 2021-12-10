function [ Error ] = LorentzianError( params, X, Y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        Error = sum((Y - (params(1)./(X.^2 + params(2)) + params(3))).^2);
end

