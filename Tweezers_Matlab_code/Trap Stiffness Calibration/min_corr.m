function M = min_corr(parameters,Px,Py,Pxy);

b = parameters(1);
c = parameters(2);

M =  ((1 + b*c) .* Pxy + c*Px + b*Py) ./ sqrt((Px + 2*b*Pxy + b^2*Py) .* (Py + 2*c*Pxy + c^2*Px));
