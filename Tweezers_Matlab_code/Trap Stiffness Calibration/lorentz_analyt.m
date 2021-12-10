function [fc,D,sfc,sD,Pfit] = lorentz_analyt(f,P,nblock);
%Finds the cornerfrequency and Diffusion constant from closed formulas
%The resulting values are used as initial parameters for the fit
global T datVr

for p = 0 : 2
    for q = 0 : 2
        eval(['s' num2str(p) num2str(q) ' = sum ((f .^ (2*' num2str(p) ')) .* (P .^ (' num2str(q) ')));']);
    end;
end;

a   = (s01 * s22 - s11 * s12) / (s02 * s22 - s12.^2);
b   = (s11 * s02 - s01 * s12) / (s02 * s22 - s12.^2);

fc  = sqrt(a/b);
D   = (1/b) * 2 * (pi.^2) * datVr.Nsplit/(datVr.Nsplit+1);

Pfit = 1 ./ (a + b .* f.^2);
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% S.D.

x   = min(f) / fc;
y   = max(f) / fc;
s   = sqrt(pi) * ( (2*y) / (1 + y^2) - (2*x) / (1 + x^2) + 2 * atan((y - x) / (1 + x*y)) - ...
        (4/(y - x)) * (atan((y - x) / (1 + x*y)))^2) ^ (-1/2); 

sfc = s * fc / sqrt(pi * fc * T);

g   = sqrt( ((2*y)/(1 + y^2)-(2*x)/(1 + x^2) + 2*atan((y - x) / (1 + x*y)) )/((1 + pi/2)*(y - x)) );

sD  = D * sqrt( (1 + pi/2) / (pi * fc * T) )*g*s;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
