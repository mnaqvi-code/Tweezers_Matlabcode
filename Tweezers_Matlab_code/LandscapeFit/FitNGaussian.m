clear all;
close all;
clc;

ext = load('ExtensionDistrubution.mat');
ext = ext.AAA*1E9;
kB = 1.3806488*10^(-23); %J/K
T = 294.15; %K (26 C )
%ext = ext(250:end-200);
figure;
plot(ext);
[n,X] = hist(ext,150);
FXU = 10E-12; %pN;
FXI = 9E-12; %pN;
FXN = 8E-12; %pN;
f = fit(X',n','gauss3')
figure;
%hist(ext,150)
hold on;
plot(f,X',n')