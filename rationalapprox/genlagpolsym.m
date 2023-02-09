function numpol = genlagpolsym(n)
%GENLAGPOLSYM Generates the poles for sinc(z) via the generalized Laguerre
%polynomials and doubles them to enforce symmetry.

p = @(n,x) laguerreL(n,-2*n-2,2*1i*x);
syms x
pol = root(p(n,x),x);
numpol = arrayfun(@(x) sym2poly(x),vpa(pol)).';
numpol = [numpol,conj(numpol)];

end