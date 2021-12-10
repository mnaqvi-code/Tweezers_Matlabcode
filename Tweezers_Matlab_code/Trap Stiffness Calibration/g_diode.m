function g = g_diode(f,fdiode,alpha);

if nargin == 1 | isempty(fdiode),
        g   =   1;
elseif nargin == 2 | isempty(alpha),
    g   =   1 ./ (1 + (f/fdiode).^2);
else,
    g   =   alpha^2 + (1-alpha^2) ./ (1 + (f/fdiode).^2);
end;