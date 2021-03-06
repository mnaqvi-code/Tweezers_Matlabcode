function y = std_nan(x)



%%Protect against bad input
%**********************************************
if (isnan(x) == 1)
    y = NaN;
    return
end

if isempty(x) % if empty return NaN
    y = NaN;
    return
end

if (size(x,1)> 1 & size(x,2) > 1) % check for matrix input
    disp('mean_nan requires a vector input and cannot handle matrix input')
    return
end

%**********************************************




not_nums = isnan(x);
index = find(not_nums);

% Find mean
mean_value = mean_nan(x);

   length_x = length(x)-sum(not_nums);
   x = x - mean_value;
   
% Replace NaNs with zeros.
x(index) = zeros(size(index));

% Protect against a column of all NaNs
y = sqrt(sum(x.*x)./max(length_x-1,1));

