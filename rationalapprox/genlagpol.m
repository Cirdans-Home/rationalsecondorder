function numpol = genlagpol(n)
%GENLAGPOL Generates the poles for sinc(z) via the generalized Laguerre
%polynomials.

p = @(n,x) laguerreL(n,-2*n-2,2*1i*x);
syms x
pol = root(p(n,x),x);
numpol = arrayfun(@(x) sym2poly(x),vpa(pol)).';

end