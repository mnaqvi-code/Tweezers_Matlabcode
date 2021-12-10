function Ph = Phydro(f,R,l,nu,rho,rhol,fc,D);

global T


gamma   =   6 * pi * rho * nu * R;
m_star  =   4/3 * R^3 * rho + 2/3 * pi * R^3 * rhol;
fm      =   gamma / (2 * pi * m_star);
fnu     =   nu / (pi * R^2);

sfnu    =   sqrt(abs(f)/fnu);


Re1 = 1 + sfnu;

Im1 = -sfnu - 2/9*abs(f)/fnu;

Re2 = 1 - 9/16 * R/l*(1 - sfnu/3 - 4/3*(1 - exp(-(2*l-R)*sfnu/R).*cos((2*l-R)*sfnu/R )));

Im2 = - 9/16 * R/l*(sfnu/3 + 2/9*abs(f)/fnu + 4/3*exp(-(2*l-R)*sfnu/R).*sin((2*l-R)*sfnu/R ));

Norm2 = Re2.*Re2 + Im2.*Im2;

Reg = (Re1.*Re2 + Im1.*Im2)./Norm2;

Img = (Re1.*Im2 - Im1.*Re2)./Norm2;


Ph  =   (D /(2*pi^2) * Reg) ./ ((fc + f.*Img - f.^2/fm).^2 + (f .* Reg).^2);