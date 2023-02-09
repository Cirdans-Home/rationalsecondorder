function numpol = genlagexppol(n)
%GENLAGEXPPOL Generates the poles for sinc(z) via the generalized Laguerre
%polynomials in for the exponential function

p = @(n,x) laguerreL(n,-2*n-1,1i*x);
syms x
pol = root(p(n,x),x);
numpol = arrayfun(@(x) sym2poly(x),vpa(pol)).';
numpol = [numpol,0,conj(numpol)];

end