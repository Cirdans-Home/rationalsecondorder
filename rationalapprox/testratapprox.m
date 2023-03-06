function results = testratapprox(A,v,ytrue,poles)
%%TESTRATAPPROX Perform tests on the rational approximation to the sinc
%%function for computing y = sinc(A)v within a Rational-Krylov method built
%%on the given poles.

numpoles = length(poles.poles);
results.times = zeros(numpoles,1);
results.abserrors = zeros(numpoles,1);
results.relerrors = zeros(numpoles,1);
results.numpoles = zeros(numpoles,1);
normy = norm(ytrue,2);

for i=1:numpoles
    tic;
    y = sincm(A,v,poles.poles{i});
    results.times(i) = toc;
    results.abserrors(i) = norm(y-ytrue,2);
    results.relerrors(i) = results.abserrors(i)./normy;
    results.numpoles(i) = length(poles.poles{i});
end


end

function P = msincf(A)
%%MSINCF dense evaluation of matrix sinc function
P = funm(A,@sin);
P = A\P;
end