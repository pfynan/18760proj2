function [f_prime_out] = dfunc(x,state)

delta = 0.01;

f_prime_out = zeros(size(x,1),1);

for i = 1:size(x,1)
    xtmp = x;
    x(i) = x(i) + delta;
    f_prime_out(i) = (func(xtmp,state) - func(x,state))/delta;
end


%here I tried to maximise the function
%f_prime_out = -1.*f_prime_out';
