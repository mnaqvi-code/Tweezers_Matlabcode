function [ force ] = calcEWLC( x, KB, T, P, L, K  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% if length(P) == length(K) && length(P) > 1
%     for ii = 1:length(P)
%         a(ii,:) = (1+ P(ii)*K(ii)./(KB*T))./K(ii).^3;
%         b(ii,:)  = ((0.25 - x./L)./K(ii).^2) + 2*(1-x./L).*(P(ii)./(KB*T) + 1./K(ii))./K(ii);
%         c(ii,:)  = (P(ii)./(KB*T) + 1./K(ii)).*((1-x./L).^2) - 2*(1-x./L).*((-0.25 + x./L)./K(ii));
%         d(ii,:)  = -1.*(-0.25 + x./L).*((1-x./L).^2) - 0.25;
%     end
%     a = sum(a);
%     b = sum(b);
%     c = sum(c);
%     d = sum(d);
% else
    a = (1+ P*K./(KB*T))./K.^3;
    b = ((0.25 - x./L)./K.^2) + 2*(1-x./L).*(P./(KB*T) + 1./K)./K;
    c = (P./(KB*T) + 1./K).*((1-x./L).^2) - 2*(1-x./L).*((-0.25 + x./L)./K);
    d = -1.*(-0.25 + x./L).*((1-x./L).^2) - 0.25;
% end

for ii = 1:length(x)
    yTemp = roots([a b(ii) c(ii) d(ii)]);
    if imag(yTemp(1,1)) == 0,
        force(ii) = yTemp(1,1);
    else
        force(ii) = yTemp(3,1);
    end
end

end

