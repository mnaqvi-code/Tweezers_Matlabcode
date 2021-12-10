function [f,P,T] = calc_powersp(X,sampling_f)

global T datVr fNyq delta_t

fNyq    =   sampling_f / 2;
delta_t =   1 / sampling_f;
X(end-mod(length(X),datVr.Nsplit)+1:end) = [];
if mod(length(X)/datVr.Nsplit,2) ~= 0
    X(end-datVr.Nsplit+1:end) = [];
end
for ii = 1:datVr.Nsplit

time{ii}    =   [ceil(length(X)*(ii-1)/datVr.Nsplit)*delta_t : delta_t : (ceil(length(X)*ii/datVr.Nsplit)-1)*delta_t]';
Tmax{ii}       =   max(time{ii})-min(time{ii});
f{ii}       =   (([ceil(length(X)*(ii-1)/datVr.Nsplit) : floor(length(X)*ii/datVr.Nsplit)+1] - ceil(length(X)*(ii-1)/datVr.Nsplit)) / (Tmax{ii}))';

FT{ii}      =   delta_t*fft(X(ceil(length(X)*(ii-1)/datVr.Nsplit)+1 : floor(length(X)*ii/datVr.Nsplit)));
P{ii}       =   FT{ii} .* conj(FT{ii}) / Tmax{ii};

ind     =   find(f{ii} <= fNyq); % only to the Nyquist f 
f{ii}       =   f{ii}(ind);
P{ii}       =   P{ii}(ind);

end
f = mean([f{:}],2);
P = mean([P{:}],2);
T = mean([Tmax{:}],2);

% global T
% 
% fNyq    =   sampling_f / 2;
% delta_t =   1 / sampling_f;
% 
% time    =   [0 : delta_t : (length(X)-1)*delta_t]';
% T       =   max(time);
% f       =   ([0 : length(X)-1] / (T+delta_t))';
% 
% FT      =   delta_t*fft(X);
% P       =   FT .* conj(FT) / T;
% 
% ind     =   find(f <= fNyq); % only to the Nyquist f 
% f       =   f(ind);
% P       =   P(ind);
%keyboard
   
