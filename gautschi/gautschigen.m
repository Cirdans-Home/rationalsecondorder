function [T,Y,mfunoptions] = gautschigen(odefun,t,yzero,yprime,A,h,mfunoptions)
%GAUTSCHIGEN General integrator of Gautschi type for a (system) of second
%order ODE of the form:
%       y''(t) = -Ay(t) + g(t)      t \in [t(1),t(2)]
%       y(t(1)) = yzero
%       y(t(1)) = yprime
%   INPUT:  odefun = @(t,y) Dynamics of the system
%           t either the time interval [t(1),t(2)] or the sequence of time
%           steps
%           yzero,yprime Initial data for state and velocity
%           A matrix
%           h time-step if t is a sequence is overwritten
%           mfunoptions options for the matrix function routine, if the
%           routine gautschigen is called without inputs in generates a
%           suitable mfunoptions with all the relevant fields
%   OUTPUT: [T,Y] time vector for the evaluated solution and solution
%               vector
%           mfunoptions matrix-function options used

%% No input arguments
if nargin == 0
    T = [];
    Y = [];
    mfunoptions.type = "direct";
    mfunoptions.poltype = "";       % Used only if type is rational
    mfunoptions.numpoles = 15;
    mfunoptions.expsumterms = 15;
    return
end


%% Check the inputs
if ~iscolumn(yzero)
    yzero = yzero.';
end
if ~iscolumn(yprime)
    yprime = yprime.';
end

hbar = waitbar(0,"Integration underway");

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

if strcmpi(mfunoptions.type,"direct")
    psiA = psi(h^2*A);
    % Initial time step
    V(:,1) = sigma(h^2*A)*V(:,1) + 0.5*h*( psiA*odefun(T(1),Y(:,1)) );
    for i=1:nt-1 %% Time loop
        waitbar(i/(nt-1),hbar,sprintf("Integration underway: %d of %d",i,nt-1));
        % Matrix function computation
        if abs((T(i+1)-T(i))-h) > 10*eps
            h = T(i+1)-T(i);
            psiA = psi(h^2*A);
        end
        % Time marching done with just two step and a single matrix-function
        % evaluation:
        Y(:,i+1) = Y(:,i) + h*V(:,i);
        V(:,i+1) = V(:,i) + h*psiA*odefun(T(i+1),Y(:,i+1));
    end
elseif strcmpi(mfunoptions.type,"rational")
    % We use here the Rational Krylov method. First of all we have to
    % decide what poles we want to use:
    switch upper(mfunoptions.poltype)
        case "PADEEXP"
            poles = genlagexppol(mfunoptions.numpoles);
        case "PADEHYPERGEOM"
            poles = genlagpolsym(ceil(mfunoptions.numpoles/2));
    end
    if iscolumn(poles)
        % rat_krylov requires poles to be a row vector!
        poles = poles.';
    end
    % Computed the poles we can proceed with the time stepping
    havev1 = false;
    havey1 = false;
    if norm(V(:,1)) > 10*eps
        Vk = rat_krylov(h^2*A, V(:,1), poles); % Build rational krylov subspace
        havev1 = true;
    end
    ytemp = odefun(T(1),Y(:,1));
    if norm(ytemp) > 10*eps
        Wk = rat_krylov(h^2*A, ytemp, poles);
        havey1 = true;
    end

    % Initial time step
    if havev1 && havey1
        V(:,1) = Vk*(sigma(Vk'*(h^2*A)*Vk)*(Vk'*V(:,1))) + 0.5*h*...
            Wk*(psi(h^2*(Wk'*A*Wk))*(Wk'*ytemp));
    elseif havey1 && ~havev1 
        V(:,1) = 0.5*h*Wk*(psi(h^2*(Wk'*A*Wk))*(Wk'*ytemp));
    elseif ~havey1 && havev1
        V(:,1) = Vk*(sigma(Vk'*(h^2*A)*Vk)*(Vk'*V(:,1)));
    end
    for i=1:nt-1 %% Time loop
        waitbar(i/(nt-1),hbar,sprintf("Integration underway: %d of %d",i,nt-1));
        % Matrix function computation
        h = T(i+1)-T(i);
        % Time marching done with just two step and a single matrix-function
        % evaluation:
        Y(:,i+1) = Y(:,i) + h*V(:,i);
        % New Rational Krylov Computation
        ytemp = odefun(T(i+1),Y(:,i+1));
        Wk = rat_krylov(h^2*A, ytemp, poles);
        % Update
        V(:,i+1) = V(:,i) + h*...
            Wk*(psi(h^2*(Wk'*A*Wk))*(Wk'*ytemp));
    end

    mfunoptions.numpoles = length(poles);
elseif strcmpi(mfunoptions.type,"expsum")
    % This branch uses the exponential-sum approximation of the sinc
    % function for the computation.

    % First of all we compute the poles to be used in the computation of
    % the various projections in the algorithm
    [x,w] = legpts(mfunoptions.expsumterms);
    [xpsi,wpsi] = legpts(mfunoptions.expsumterms,[-2,0]);
    poles = genpadeexppol(mfunoptions.numpoles).';

    % Initial time step
    V(:,1) = expsumsigma(h^2*A,V(:,1),poles,x,w) ....
        + 0.5*h*( expsumpsi(h^2*A,odefun(T(1),Y(:,1)),poles,xpsi,wpsi) );
    for i=1:nt-1 %% Time loop
        waitbar(i/(nt-1),hbar,sprintf("Integration underway: %d of %d",i,nt-1));
        h = T(i+1)-T(i);
        % Time marching done with just two step and a single matrix-function
        % evaluation:
        Y(:,i+1) = Y(:,i) + h*V(:,i);
        V(:,i+1) = V(:,i) + ...
            h*expsumpsi(h^2*A,odefun(T(i+1),Y(:,i+1)),poles,xpsi,wpsi);
    end

else
    close(hbar)
    error("Unkwnown matrix-function type: %d it has to be: 'direct'," + ...
        "'rational' or 'expsum'",mfunoptions.type);
end

close(hbar)

end

%% Auxiliary routines for the computation of the matrix function
function S = sigma(A)
%SIGMA computes S = sqrt(A)^{-1} sin(sqrt(A))
if issparse(A)
    A = full(A);
end
sqA = sqrtm(A);
S = sqA\funm(sqA,@sin);
end

function P = psi(A)
%PSI computes P = ( (sqrt(A)/2)^{-1} sin(sqrt(A)/2) )^2
if issparse(A)
    A = full(A);
end
sqA2 = sqrtm(A)/2;
P = sqA2\funm(sqA2,@sin);
P = P*P;
end
function y = expsumsigma(A,v,poles,x,w)
%EXPSUMSIGMA computes S = sqrt(A)^{-1} sin(sqrt(A)) usig the exponential
%sum approach.

y = zeros(size(v));
if norm(v,2) < 10*eps
    return
end

N = length(x);
V = rat_krylov(A, v, poles.'); % Build rational krylov subspace
vk = V'*v;
Ap = sqrtm(V'*A*V);
for i=1:N
    Ak = -1i*x(i)*Ap;
    ytemp = V*(expm(Ak)*vk);
    y = y + w(i)*ytemp;
end
y = real(0.5*y);


end

function y = expsumpsi(A,v,poles,x,w)
%EXPSUMPSI computes P = ( (sqrt(A)/2)^{-1} sin(sqrt(A)/2) )^2 using the
%exponential sum approach.

y = zeros(size(v));
if norm(v,2) < 10*eps
    return
end

N = length(x);
V = rat_krylov(A, v, poles.'); % Build rational krylov subspace
vk = V'*v;
Ap = sqrtm(V'*A*V)/2;
I = eye(size(Ap));
for i=1:N
    Ak1 = -1i*x(i)*Ap;
    Ak2 = +1i*x(N-i+1)*Ap;
    ytemp = V*((w(i)*(4 + 2*x(i))*expm(Ak1)+w(i)*(4 + 2*x(N-i+1))*expm(Ak2))*vk);
    y = y + ytemp/8;
end
y = real(y);


end
