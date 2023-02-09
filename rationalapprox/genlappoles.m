function poles = genlappoles(N,t0,t1)
%%GENLAPPOLES Poles obtained from the inversion of Laplace transform.
Lambda = t1/t0;
h = sqrt(8*Lambda+1)/N;
mu = pi*N/(4*t1*sqrt(8*Lambda+1));
s = @(u) -mu*(1i*u + 1).^2;
uk = (-N:N)*h;
poles = [1i*s(uk),-1i*s(uk)];

end