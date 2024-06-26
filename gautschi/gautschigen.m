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
    mfunoptions.verbose = false;
    mfunoptions.solvetype = "backslash";
    return
end


%% Check the inputs
if ~iscolumn(yzero)
    yzero = yzero.';
end
if ~iscolumn(yprime)
    yprime = yprime.';
end

if (mfunoptions.verbose)
    hbar = waitbar(0,"Integration underway");
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
isrealyes = false;
if isreal(yzero) && isreal(yprime)
    isrealyes = true;
end

if strcmpi(mfunoptions.type,"direct")
    psiA = psi(h^2*A);
    % Initial time step
    V(:,1) = sigma(h^2*A)*V(:,1) + 0.5*h*( psiA*odefun(T(1),Y(:,1)) );
    for i=1:nt-1 %% Time loop
        if (mfunoptions.verbose)
            waitbar(i/(nt-1),hbar,sprintf("Integration underway: %d of %d",i,nt-1));
        end
        % Matrix function computation
        if abs((T(i+1)-T(i))-h) > 10*eps
            h = T(i+1)-T(i);
            psiA = psi(h^2*A);
        end
        % Time marching done with just two step and a single matrix-function
        % evaluation:
        if isrealyes
            Y(:,i+1) = real(Y(:,i) + h*V(:,i));
            V(:,i+1) = real(V(:,i) + h*psiA*odefun(T(i+1),Y(:,i+1)));
        else
            Y(:,i+1) = Y(:,i) + h*V(:,i);
            V(:,i+1) = V(:,i) + h*psiA*odefun(T(i+1),Y(:,i+1));
        end
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
    if isfield(mfunoptions,"solvetype")
        switch upper(mfunoptions.solvetype)
            case "BACKSLASH"
                AB = h^2*A;
            case "PARDISO"
                AB = struct;
                AB.isreal = isreal(A);
                AB.A = h^2*A;
                AB.multiply = @(rho, eta, x) h^2*rho*(A*x) - eta*x;
                factorization = struct;
                xi_unique = unique(poles);
                factorization.xi_unique = xi_unique;
                info = cell(3*length(xi_unique),1);
                for i=1:length(xi_unique)
                    Aop_shifted = AB.A - xi_unique(i)*speye(size(AB.A));
                    if ~isreal(xi_unique)
                        if issymmetric(Aop_shifted)
                            info{3*i-2} = pardisoinit(6,0);
                            info{3*i-2} = pardisoreorder(tril(Aop_shifted),info{3*i-2},false);
                            info{3*i-2} = pardisofactor(tril(Aop_shifted),info{3*i-2},false);
                        else
                            info{3*i-2} = pardisoinit(13,0);
                            info{3*i-2} = pardisoreorder(Aop_shifted,info{3*i-2},false);
                            info{3*i-2} = pardisofactor(Aop_shifted,info{3*i-2},false);
                        end
                    else
                        if issymmetric(Aop_shifted)
                            info{3*i-2} = pardisoinit(2,0);
                            info{3*i-2} = pardisoreorder(tril(Aop_shifted),info{3*i-2},false);
                            info{3*i-2} = pardisofactor(tril(Aop_shifted),info{3*i-2},false);
                        else
                            info{3*i-2} = pardisoinit(11,0);
                            info{3*i-2} = pardisoreorder(Aop_shifted,info{3*i-2},false);
                            info{3*i-2} = pardisofactor(Aop_shifted,info{3*i-2},false);
                        end
                    end
                end
                factorization.info = info;
                AB.solve = @(nu,mu,b) solve_pardiso(AB.A,nu,mu,b,factorization);
            case "GMRES"
                AB = struct;
                AB.isreal = isreal(A);
                AB.A = h^2*A;
                AB.multiply = @(rho, eta, x) h^2*rho*(A*x) - eta*x;
                AB.solve = @(nu,mu,b) solve_gmres(AB.A,nu,mu,b);
            case "BICGSTAB"
                AB = struct;
                AB.isreal = isreal(A);
                AB.A = h^2*A;
                AB.multiply = @(rho, eta, x) h^2*rho*(A*x) - eta*x;
                AB.solve = @(nu,mu,b) solve_bicgstab(AB.A,nu,mu,b);
            case "PGMRES"
                AB = struct;
                AB.isreal = isreal(A);
                AB.A = h^2*A;
                AB.multiply = @(rho, eta, x) h^2*rho*(A*x) - eta*x;
                xi_unique = unique(poles);
                LP = cell(length(xi_unique));
                UP = cell(length(xi_unique));
                for i=1:length(xi_unique)
                    [LP{i},UP{i}] = ilu(AB.A - xi_unique(i)*speye(size(AB.A)));
                end
                AB.solve = @(nu,mu,b) solve_pgmres(AB.A,nu,mu,b,LP,UP,xi_unique);
            case "PBICGSTAB"
                AB = struct;
                AB.isreal = isreal(A);
                AB.A = h^2*A;
                AB.multiply = @(rho, eta, x) h^2*rho*(A*x) - eta*x;
                xi_unique = unique(poles);
                LP = cell(length(xi_unique));
                UP = cell(length(xi_unique));
                for i=1:length(xi_unique)
                    [LP{i},UP{i}] = ilu(AB.A - xi_unique(i)*speye(size(AB.A)));
                end
                AB.solve = @(nu,mu,b) solve_pbicgstab(AB.A,nu,mu,b,LP,UP,xi_unique);
            otherwise
                AB = h^2*A;
        end
    else
        AB = h^2*A;
    end
    if iscolumn(poles)
        % rat_krylov requires poles to be a row vector!
        poles = poles.';
    end
    % Computed the poles we can proceed with the time stepping
    havev1 = false;
    havey1 = false;
    if norm(V(:,1)) > 10*eps
        Vk = rat_krylov(AB, V(:,1), poles); % Build rational krylov subspace
        havev1 = true;
    end
    ytemp = odefun(T(1),Y(:,1));
    if norm(ytemp) > 10*eps
        Wk = rat_krylov(AB, ytemp, poles);
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
        if (mfunoptions.verbose)
            waitbar(i/(nt-1),hbar,sprintf("Integration underway: %d of %d",i,nt-1));
        end
        % Matrix function computation
        h = T(i+1)-T(i);
        % Time marching done with just two step and a single matrix-function
        % evaluation:
        Y(:,i+1) = Y(:,i) + h*V(:,i);
        % New Rational Krylov Computation
        ytemp = odefun(T(i+1),Y(:,i+1));
        Wk = rat_krylov(AB, ytemp, poles);
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
        if (mfunoptions.verbose)
            waitbar(i/(nt-1),hbar,sprintf("Integration underway: %d of %d",i,nt-1));
        end
        h = T(i+1)-T(i);
        % Time marching done with just two step and a single matrix-function
        % evaluation:
        Y(:,i+1) = Y(:,i) + h*V(:,i);
        V(:,i+1) = V(:,i) + ...
            h*expsumpsi(h^2*A,odefun(T(i+1),Y(:,i+1)),poles,xpsi,wpsi);
    end

else
    if (mfunoptions.verbose)
        close(hbar)
    end
    error("Unkwnown matrix-function type: %d it has to be: 'direct'," + ...
        "'rational' or 'expsum'",mfunoptions.type);
end
if (mfunoptions.verbose)
    close(hbar)
end

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
T = schur(A);
if any(abs(diag(T)) < 10*eps ) && issymmetric(A)
    [V,D,W] = eig(A);
    P = V*( diag(sincsq(diag(D))) )*W';
else
    sqA2 = sqrtm(T)/2;
    P = sqA2\funm(sqA2,@sin);
    P = P*P;
end

end

function y = sincsq(d)
y = d;
for i=1:length(d)
    if abs(d(i)) < 10*eps
        y(i) = 1;
    else
        y(i) = (sin(sqrt(d(i))/2)/sqrt(d(i))/2)^2;
    end
end
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
%I = eye(size(Ap));
for i=1:N
    Ak1 = -1i*x(i)*Ap;
    Ak2 = +1i*x(N-i+1)*Ap;
    ytemp = V*((w(i)*(4 + 2*x(i))*expm(Ak1)+w(i)*(4 + 2*x(N-i+1))*expm(Ak2))*vk);
    y = y + ytemp/8;
end
y = real(y);


end

function y = solve_pardiso(A,nu,mu,b,factorization)
%SOLVE_PARDISO Solve shifted linear system with PARDISO

if nu ==0 && mu==1
    % Linear system solve with identity
    y = b;
else
    % Solve with PARDISO
    if issymmetric(A)
            [y, ~] = pardisosolve(tril(nu*A - mu*speye(size(A))), b, ...
            factorization.info{3*find(factorization.xi_unique==mu)-2}, ...
            false);
    else
        [y, ~] = pardisosolve(nu*A - mu*speye(size(A)), b, ...
            factorization.info{3*find(factorization.xi_unique==mu)-2}, ...
            false);
    end
end

end

function y = solve_gmres(A,nu,mu,b)
%SOLVE_GMRES Solve shifted linear system with GMRES

if nu ==0 && mu==1
    % Linear system solve with identity
    y = b;
else
    % Solve with GMRES
    matvec = @(x) nu*(A*x) - mu*x;
    [y, ~] = gmres(@(x) matvec(x), b,[],1e-6,100);
end

end

function y = solve_bicgstab(A,nu,mu,b)
%SOLVE_BICGSTAB Solve shifted linear system with BiCGStab

if nu ==0 && mu==1
    % Linear system solve with identity
    y = b;
else
    % Solve with BICGSTAB
    matvec = @(x) nu*(A*x) - mu*x;
    [y, ~] = bicgstab(@(x) matvec(x), b,1e-6,100);
end

end

function y = solve_pgmres(A,nu,mu,b,LP,UP,xi_unique)
%SOLVE_GMRES Solve shifted linear system with preconditioned GMRES

if nu ==0 && mu==1
    % Linear system solve with identity
    y = b;
else
    % Solve with GMRES
    matvec = @(x) nu*(A*x) - mu*x;
    [y, ~] = gmres(@(x) matvec(x), b,[],1e-6,100,LP{xi_unique==mu},...
        UP{xi_unique==mu});
end

end

function y = solve_pbicgstab(A,nu,mu,b,LP,UP,xi_unique)
%SOLVE_BICGSTAB Solve shifted linear system with preconditioned BiCGStab

if nu ==0 && mu==1
    % Linear system solve with identity
    y = b;
else
    % Solve with BICGSTAB
    matvec = @(x) nu*(A*x) - mu*x;
    [y, ~] = bicgstab(@(x) matvec(x), b,1e-6,100,LP{xi_unique==mu},...
        UP{xi_unique==mu});
end

end

