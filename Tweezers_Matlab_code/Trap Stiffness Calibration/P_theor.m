function [P_theor, alias_check] = P_theor(scal_fit,parameters,x,check_flag);
global n fNyq T want_alias want_hydro R l nu rho rhol want_filter elec_filters nblock sampling_f delta_t
check_flag; % is set to 1 when we want to check if more aliasing terms should be included, otherwise check_flag = 0
alias_check = 0;
fc      =   parameters(1)*scal_fit(1);
D       =   parameters(2)*scal_fit(2);
if length(parameters) > 2, fdiode = parameters(3)*scal_fit(3); else fdiode = []; end;
if length(parameters) > 3, a0 = parameters(4)*scal_fit(4); else alpha = []; end;
if length(parameters) > 3, alpha = 1 / sqrt(1 + a0^2); else alpha = []; end;

if want_filter == 1,
    if isempty(findstr('f',elec_filters)), errordlg('Missing ''f'' in the filter function!'); end;
    g_elec  =   inline(elec_filters,'f');
end; %if want_filter == 1

gamma   =   6 * pi * rho * nu * R;
m_star  =   4/3 * R^3 * rho + 2/3 * pi * R^3 * rhol;
fm      =   gamma / (2 * pi * m_star);
fnu     =   nu / (pi * R^2);

if (check_flag == 0 & n == 0)
    n = 10;         %Default value for number of aliasing terms
end

ll = length([-n:n]);
xx = zeros(length(x),ll);
for i= 1: length(x)
    xx(i,:) = x(i) +2*[-n:n]*fNyq;
    if (check_flag == 1)
        xx_check(i,:) = x(i) +2*[-(n+1) n+1]*fNyq;
    end
end
switch want_hydro;
    case 0;
        switch want_filter;
            case 0;
                
                if want_alias == 1,
                    P_theor = sum((D/(2*pi^2)) ./ (xx.^2 + fc^2) .* g_diode(xx, fdiode, alpha),2);
                    if (check_flag == 1)
                        %P_theor_check = sum((1-(exp(-pi*fc/fNyq))^2)*D./(2*pi*fc).*delta_t./(1+(exp(-pi*fc/fNyq)).^2-2.*(exp(-pi*fc/fNyq)).*cos(2.*pi.*(x)./(sampling_f*1000))));%; .* g_diode(x, fdiode, alpha),2);
                        P_theor_check = sum((D/(2*pi^2)) ./ (xx_check.^2 + fc^2) .* g_diode(xx_check, fdiode, alpha),2);
                        alias_check = mean(P_theor_check(end-5:end)./P_theor(end-5:end));
                    end
                else
                    cP = exp(-pi*fc/fNyq);
                    %P_theor = (1-(exp(-pi*fc/fNyq))^2)*D./(2*pi*fc).*delta_t./(1+(exp(-pi*fc/fNyq)).^2-2.*(exp(-pi*fc/fNyq)).*cos(2.*pi.*(x)./(sampling_f*1000)));%; .* g_diode(x, fdiode, alpha);
                    P_theor =  (D/(2*pi^2)) ./ (x.^2 + fc^2) .* g_diode(x, fdiode, alpha); %.*(1 +(1 ./ (1 + (x/17000).^16)) .* ((((1-(exp(-pi*fc/fNyq))^2)*D./(2*pi*fc).*delta_t./(1+(exp(-pi*fc/fNyq)).^2-2.*(exp(-pi*fc/fNyq)).*cos(2.*pi.*(x)./(sampling_f*1000))))./((D/(2*pi^2)) ./ (x.^2 + fc^2))) - 1));
                end; %want_alias == 1
                
            case 1;
                if want_alias == 1,
                    P_theor = sum((D/(2*pi^2)) ./ (xx.^2 + fc^2) .* g_diode(xx, fdiode, alpha) .*...
                        g_elec(xx),2);
                    if (check_flag == 1)
                        P_theor_check = sum((D/(2*pi^2)) ./ (xx_check.^2 + fc^2) .* g_diode(xx_check, fdiode, alpha) .*...
                            g_elec(xx_check),2);
                        alias_check = mean(P_theor_check(end-5:end)./P_theor(end-5:end));
                    end
                else                    
                    cP = (exp(-pi*fc/fNyq));
                    %P_theor = (1-(exp(-pi*fc/fNyq))^2)*D./(2*pi*fc).*delta_t./(1+(exp(-pi*fc/fNyq)).^2-2.*(exp(-pi*fc/fNyq)).*cos(2.*pi.*(x)./(sampling_f*1000)));% .* g_diode(x, fdiode, alpha);
                    P_theor = (D/(2*pi^2)) ./ (x.^2 + fc^2) .* g_diode(x, fdiode, alpha) .* g_elec(x);
                end; %want_alias == 1
        end; % switch want_filter
        
        
    case 1;
        switch want_filter;
            case 0;
                sfnu    =   sqrt(abs(xx)/fnu);
                Reg =   1 + sfnu - 3/16 * R/l + 3/4 * R/l * exp(-2*l/R*sfnu) .* cos(2*l/R*sfnu);
                Img =   -sfnu + 3/4 * R/l * exp(-2*l/R*sfnu) .* sin(2*l/R*sfnu);
                Ph  =   (D /(2*pi^2) * Reg) ./ ((fc + xx.*Img - xx.^2/fm).^2 + (xx .* Reg).^2);
                if want_alias == 1,
                    P_theor = sum(Ph .* g_diode(xx,fdiode,alpha),2);
                    if (check_flag == 1)
                        sfnu    =   sqrt(abs(xx_check)/fnu);
                        Reg =   1 + sfnu - 3/16 * R/l + 3/4 * R/l * exp(-2*l/R*sfnu) .* cos(2*l/R*sfnu);
                        Img =   -sfnu + 3/4 * R/l * exp(-2*l/R*sfnu) .* sin(2*l/R*sfnu);
                        Ph  =   (D /(2*pi^2) * Reg) ./ ((fc + xx_check.*Img - xx_check.^2/fm).^2 + (xx_check .* Reg).^2);
                        P_theor_check = sum(Ph .* g_diode(xx_check,fdiode,alpha),2);
                        alias_check = mean(P_theor_check(end-5:end)./P_theor(end-5:end));
                    end
                else
                    P_theor = P_hydro(x, R, l, nu, rho, rhol, fc, D) .* g_diode(x, fdiode, alpha);
                end; %want_alias == 1
                
                
            case 1;
                sfnu    =   sqrt(abs(xx)/fnu);
                Reg =   1 + sfnu - 3/16 * R/l + 3/4 * R/l * exp(-2*l/R*sfnu) .* cos(2*l/R*sfnu);
                Img =   -sfnu + 3/4 * R/l * exp(-2*l/R*sfnu) .* sin(2*l/R*sfnu);
                Ph  =   (D /(2*pi^2) * Reg) ./ ((fc + xx.*Img - xx.^2/fm).^2 + (xx .* Reg).^2);
                if want_alias == 1,
                    P_theor = sum(Ph .* g_diode(xx,fdiode,alpha) .* g_elec(xx),2);
                    if (check_flag == 1)
                        sfnu    =   sqrt(abs(xx_check)/fnu);
                        Reg =   1 + sfnu - 3/16 * R/l + 3/4 * R/l * exp(-2*l/R*sfnu) .* cos(2*l/R*sfnu);
                        Img =   -sfnu + 3/4 * R/l * exp(-2*l/R*sfnu) .* sin(2*l/R*sfnu);
                        Ph  =   (D /(2*pi^2) * Reg) ./ ((fc + xx_check.*Img - xx_check.^2/fm).^2 + (xx_check .* Reg).^2);
                        P_theor_check = sum(Ph .* g_diode(xx_check,fdiode,alpha) .* g_elec(xx_check),2);
                        alias_check = mean(P_theor_check(end-5:end)./P_theor(end-5:end));
                    end
                else
                    P_theor = P_hydro(x, R, l, nu, rho, rhol, fc, D) .* g_diode(x, fdiode, alpha) .* g_elec(x);
                end; %want_alias == 1
                
        end; % switch want_filter
        
end; % switch want_hydro
P_theor =   P_theor(:);
