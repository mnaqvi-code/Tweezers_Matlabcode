function force = alirezawlcdiffpers(x, P, s)
%Modified Marko-Siggia WLC 
% SI units except for extension (um) and force (pN) 
% P = persistence length
% L = contour length
% KB = Boltzmann constant
% T = temp
% K = elastic modulus
% Refrenece: Steven Bock, Stretching DNA with optical tweezer. Biophys J 72, 1335-1346, 1997.

x = x*(10^-6);
KB = 1.3806503*10^(-23);

if s ==1,
    % dsDNA
    K = 1.2*10^(-9);  % Ruud thesis and others
    P = 45*10^(-9);  
   % P = 54*10^(-9);  % Ruud thesis  

    L =  920*10^(-9);   % 0.34*10^(-9)*2500; where 2500 is the number of bases
    
elseif s == 3,
    % 2*ssDNA : typical values for persistence length are either 0.75 or
    % 0.85
    
K = 0.8*10^(-9);

P = 0.75*10^(-9);

L =  2*0.7*10^(-9)*2500;
   
elseif s == 2, 
    % S-DNA
%    Ref:
%    Europhys. Lett., 62 (5), p. 760 (2003)
%    The bend stiffness of S-DNA
 K = 0.8*10^(-9);  % taken the value for ssDNA
 P = 12.32*10^(-9);
 L = 1.7*0.34*10^(-9)*2500;  
else 
% Protein:   
% For elongated peptide backbone, the young modulus is 50 GPa (Journal of Applied Physics, Vol. 90, No. 6, pp. 3095–3099, 15 September 2001)
%      MBP 360  
%      Myc 10 
%      4 Myc 40 

 K = 27*10^(-9);  % ref; PRL 94, 048301 (2005)
 % K = 250*1.2*10^(-9); % my very rough guess based on the fact that cross sectional area of the peptide is only a few angstrom square and: An AFM study of the elasticity of DNA molecules , Thin Solid films 2004.
 %K = 10^(-9); % Engineering Analysis with Boundary Elements. Volume 31, Issue 5, May 2007, Pages 402-409 
 %K = 1.2*10^(-9); % Ruud's thesis  
 P = 1.5*10^(-9);   % 0.5-2 nm Ruud thesis, Butsamante: 2 pN
% L =  120*10^(-9);    % Ruud thesis: 120 nm  i.e. 0.33 per aa and 370 aa  / other references I have seen that people use 0.36 per aa 

  L = (181.5)*10^(-9); %Luciferase
  
%L = 445 *10^(-9); % [160, 225, 250, 445]
%L =  28*10^(-9); for one c-terminal 
%   fourMBPL = ([112, 204, 296, 388, 480]+15).*10^(-9);  
%   L = fourMBPL(1);  
  

end
T = 293;

a = (1+ P*K./(KB*T))./K.^3;
b = ((0.25 - x./L)./K.^2) + 2*(1-x./L).*(P./(KB*T) + 1./K)./K;
c = (P./(KB*T) + 1./K).*((1-x./L).^2) - 2*(1-x./L).*((-0.25 + x./L)./K);
d = -1.*(-0.25 + x./L).*((1-x./L).^2) - 0.25;

for ix = 1:length(b),
    
[root,nroot] = alirezacubicfcn(a, b(ix), c(ix), d(ix));
if imag(root(1)) == 0,
force(ix) = (10^12).*root(1);
else 
    force(ix) = (10^12).*root(3);
end
    
root = [];
end
keyboard;


