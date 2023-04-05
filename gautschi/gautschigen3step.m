function [T,Y] = gautschigen3step(odefun,t,yzero,yprime,A,h)
%GAUTSCHIGEN3STEP General integrator of Gautschi type for a (system) of second
%order ODE of the form:
%       y''(t) = -Ay(t) + g(t)      t \in [t(1),t(2)]
%       y(t(1)) = yzero
%       y(t(1)) = yprime

%% Check the inputs
if ~iscolumn(yzero)
    yzero = yzero.';
end
if ~iscolumn(yprime)
    yprime = yprime.';
end

%% Generate the time-grid
if length(t) > 2
    T = t;
    h = T(2)-T(1);
    nt = length(T);
elseif length(t) == 2
    T = t(1):h:t(2);
    nt = length(T);
else
    error("t has to be length 2 or more")
end

%% Initial conditions
ndof = length(yzero);
Y = zeros(ndof,nt);
V = zeros(ndof,nt);

Y(:,1) = yzero;
V(:,1) = yprime;
V(:,1) = sigma(h^2*A)*V(:,1);

psiA = psi(h^2*A);
%% Time loop
for i=1:nt-1
    % Matrix function computation
    if abs((T(i+1)-T(i))-h) > 10*eps
        h = T(i+1)-T(i);
        psiA = psi(h^2*A);
    end
    % Time marching
    vtemp    = V(:,i) + 0.5*h*psiA*odefun(T(i),Y(:,i));
    Y(:,i+1) = Y(:,i) + h*vtemp;
    V(:,i+1) = vtemp + 0.5*h*psiA*odefun(T(i+1),Y(:,i+1));
end

end

%% Auxiliary routines for the computation of the matrix function
function S = sigma(A)
%%SIGMA computes S = sqrt(A)^{-1} sin(sqrt(A))
sqA = sqrtm(A);
S = sqA\funm(sqA,@sin);
end

function P = psi(A)
%%PSI computes P = ( (sqrt(A)/2)^{-1} sin(sqrt(A)/2) )^2 
sqA2 = sqrtm(A)/2;
P = sqA2\funm(sqA2,@sin);
P = P*P;
end