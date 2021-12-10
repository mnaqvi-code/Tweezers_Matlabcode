function hh = funn(scal_fit,parameters,xfin,yfin,sfin,check_flag);
%The function minus the data divided by the error sfin
hh = (1./P_theor(scal_fit,parameters,xfin,check_flag) - 1./yfin) ./ sfin;

