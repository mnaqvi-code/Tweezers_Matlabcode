function y = free_gauss(beta,x);

mu    = beta(1);
sigma = beta(2);
A     = beta(3);

y = (A/(sqrt(2*pi)*sigma)) * exp(- ((x-mu)/sigma).^2 / 2);

