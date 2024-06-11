function [T,Y,krylsize] = polgautschi(odefun,t,yzero,yprime,A,h)
%POLGAUTSCHI Gautschi scheme with the Krylov subspace matrix function
% evaluation and adaptive choice of the Krylov dimension as in
% Botchev, M. A.; Harutyunyan, D.; van der Vegt, J. J. W.
% The Gautschi time stepping scheme for edge finite element
% discretizations of the Maxwell equations. J. Comput. Phys.216(2006),
% no.2, 654–686.

% Allocate solution space
n = length(yzero);
T = t(1):h:t(2);
k = length(T);
Y = zeros(n,k);
krylsize = zeros(k-2,1);

Y(:,1) = yzero;
% Need to compute the second step outside the loop, for this we use one
% step of the standard leap-frog integrator:
vtemp = sigma_matvec(h^2*A,yprime,h^2) + ...
    0.5*h*psi_matvec(h^2*A,odefun(T(1),Y(:,1)),h^2);
Y(:,2) = Y(:,1) + h*vtemp;
% Loop using the Polynomial Krylov routine
if issymmetric(A)
    for i=2:k-1
        [Y(:,i+1),krylsize(i)] = gautschistep_lanczos(A,Y(:,i),Y(:,i-1),odefun,T(i),h);
    end
else
    for i=2:k-1
        [Y(:,i+1),krylsize(i)] = gautschistep_arnoldi(A,Y(:,i),Y(:,i-1),odefun,T(i),h);
    end
end

end

function [ynp1,krylsize] = gautschistep_arnoldi(A,yn,ynm1,odefun,t,h)
%GAUTSCHISTEP performs a step of the Gautschi integrator given y^{(n)} and
%y^{(n-1)}.

mmax = min(size(A,2)-1,100);

V = zeros(size(A,1),mmax);
H = zeros(mmax+1,mmax);

v = -odefun(t,yn);              % v = A y^{(n)} - f^{(n)}
beta = norm(v,2);           % \beta = ||v||_2
ynp1_zero = 2*yn -ynm1 -h^2*v;   % y^{(n+1)}_0 = 2 y^{(n)}
for m = 1:mmax-1
    % One-step of Arnoldi
    if m == 1
        V(:,1) = v/beta;
        ynp1_old = ynp1_zero;
    else
        ynp1_old = ynp1_new;
    end
    w = h^2*A*V(:,m);
    % Arnoldi cycle
    for i = 1:m
        H(i,m) = w'*V(:,m);
        w = w - H(i,m)*V(:,i);
    end
    H(m+1,m) = norm(w,2);
    V(:,m+1) = w/H(m+1,m);
    % Matrix function evaluation
    psival = beta*psi(H(1:m,1:m));
    u = V(:,1:m)*psival(:,1);
    ynp1_new = 2*yn - ynm1 - h^2*u;
    % Compute the error estimator
    E = norm(ynp1_new - ynp1_old,"inf")/norm(ynp1_new - ynp1_zero,"inf");
    if E < min(max(h^2,10*eps),1e-6)
        break
    end
end
krylsize = i;
ynp1 = ynp1_new;

end

function [ynp1,krylsize] = gautschistep_lanczos(A,yn,ynm1,odefun,t,h)
%GAUTSCHISTEP performs a step of the Gautschi integrator given y^{(n)} and
%y^{(n-1)}.

mmax = min(size(A,2)-1,100);

V = zeros(size(A,1),mmax);
T = zeros(mmax,mmax);

v = -odefun(t,yn);              % v = A y^{(n)} - f^{(n)}
beta = norm(v,2);               % \beta = ||v||_2
beta0 = beta;
ynp1_zero = 2*yn -ynm1 -h^2*v;  % y^{(n+1)}_0 = 2 y^{(n)}

V(:,1) = v/beta;
for m = 1:mmax
    if m == 1
        ynp1_old = ynp1_zero;
    else
        ynp1_old = ynp1_new;
    end
    w = h^2*A*V(:,m);
    if m > 1
        w = w - beta*V(:,m-1);
    end
    alpha = w' * V(:,m);
    w = w - alpha*V(:,m);
    beta = norm(w,2);
    V(:,m+1) = w / beta;
    % Update tridiagonal matrix
    T(m,m) = alpha;
    T(m,m+1) = beta;
    T(m+1,m) = beta;
    % Matrix function evaluation
    psival = beta0*psi(T(1:m,1:m));
    u = V(:,1:m)*psival(:,1);
    ynp1_new = 2*yn - ynm1 - h^2*u;

    % Compute the error estimator
    E = norm(ynp1_new - ynp1_old,"inf") / norm(ynp1_new - ynp1_zero,"inf");
    if E < min(max(h^2,10*eps),1e-6)
        break
    end
end

krylsize = m;
if (krylsize == mmax)
    warning("Reached maximum number of iterations: err = %e\n",E)
end
ynp1 = ynp1_new;
end

function y = psi_matvec(A,v,tol)
%PSI_MATVEC routine that uses a polynomial Krylov method to perform the
%matrix vector product psi(A)v
if norm(v,"inf") < 10*eps
    y = zeros(size(v,1),size(v,2));
    return
end

if issymmetric(A)
    % Lanczos
    mmax = 200;
    n = length(v);
    V = zeros(n, mmax);
    T = zeros(mmax, mmax);
    beta = norm(v);
    beta0 = beta;
    V(:,1) = v / beta;
    for m = 1:mmax
        w = A * V(:,m);
        if m > 1
            w = w - beta * V(:,m-1);
        end
        alpha = V(:,m)' * w;
        w = w - alpha * V(:,m);
        beta = norm(w);
        if beta < eps
            break;
        end
        if m < mmax
            V(:,m+1) = w / beta;
        end
        T(m,m) = alpha;
        T(m,m+1) = beta;
        T(m+1,m) = beta;
        if abs(T(m,m)) < tol
            break
        end
    end
    % Adjust size of T if we stopped early
    psiT =  beta0*psi(T(1:m,1:m));
    y  = V(:,1:m)*psiT(:,1);
else
    % Arnoldi
    mmax = 200;
    n = length(b);
    V = zeros(n, mmax+1);
    H = zeros(mmax+1, mmax);

    beta0 = norm(b);
    V(:,1) = b / beta0;

    for k = 1:mmax
        w = A * V(:,k);

        for j = 1:k
            H(j,k) = V(:,j)' * w;
            w = w - H(j,k) * V(:,j);
        end

        H(k+1,k) = norm(w);

        if abs(H(k+1,k)) > tol
            V(:,k+1) = w / H(k+1,k);
        else
            break;
        end
    end
    % Adjust size of T if we stopped early
    psiH =  beta0*psi(H(1:k, 1:k));
    y  = V(:,1:k)*psiH(:,1);
end
end

function y = sigma_matvec(A,v,tol)
%SIGMA_MATVEC routine that uses a polynomial Krylov method to perform the
%matrix vector product sigma(A)v
if norm(v,"inf") < 10*eps
    y = zeros(size(v,1),size(v,2));
    return
end

if issymmetric(A)
    mmax = 200;
    n = length(v);
    V = zeros(n, mmax);
    T = zeros(mmax, mmax);
    beta = norm(v);
    beta0 = beta;
    V(:,1) = v / beta;
    for m = 1:mmax
        w = A * V(:,m);
        if m > 1
            w = w - beta * V(:,m-1);
        end
        alpha = V(:,m)' * w;
        w = w - alpha * V(:,m);
        beta = norm(w);
        if beta < eps
            break;
        end
        if m < mmax
            V(:,m+1) = w / beta;
        end
        T(m,m) = alpha;
        T(m,m+1) = beta;
        T(m+1,m) = beta;
        if abs(T(m,m)) < tol
            break
        end
    end
    % Adjust size of T if we stopped early
    psiT = beta0*sigma(T(1:m,1:m));
    y  = V(:,1:m)*psiT(:,1);
else
    % Arnoldi
    mmax = 200;
    n = length(b);
    V = zeros(n, mmax+1);
    H = zeros(mmax+1, mmax);

    beta0 = norm(b);
    V(:,1) = b / beta0;

    for k = 1:mmax
        w = A * V(:,k);

        for j = 1:k
            H(j,k) = V(:,j)' * w;
            w = w - H(j,k) * V(:,j);
        end

        H(k+1,k) = norm(w);

        if abs(H(k+1,k)) > tol
            V(:,k+1) = w / H(k+1,k);
        else
            break;
        end
    end
    % Adjust size of T if we stopped early
    psiH =  beta0*sigma(H(1:k, 1:k));
    y  = V(:,1:k)*psiH(:,1);
end
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

function S = sigma(A)
%SIGMA computes S = sqrt(A)^{-1} sin(sqrt(A))
if issparse(A)
    A = full(A);
end
sqA = sqrtm(A);
S = sqA\funm(sqA,@sin);
end