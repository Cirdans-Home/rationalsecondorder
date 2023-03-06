function numpol = genpadeexppol(n)
%GENLAGEXPPOL Generates the poles for exp(-z) via the generalized Laguerre
%polynomials in for the exponential function

p = @(n,x) laguerreL(n,-2*n-1,x);
syms x
pol = root(p(n,x),x);
numpol = arrayfun(@(x) sym2poly(x),vpa(pol)).';

end