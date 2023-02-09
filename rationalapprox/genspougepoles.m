function poles = genspougepoles(N)
%%GENSPOUGEPOLES Generates the poles for the approximation of the sinc(A)v
% based on the Spouge approximation of the Gamma(1+x) function.

poles = [pi*(1:N),-pi*(1:N)];


end