function [y,varargout] = sincfourier(A,v,N,type,exptype,cfpoles,a,b)
%SINCFOURIER computes the sinc(A)v product using the inverse Fourier
%transform i.e., y = sinc(A)v
% INPUT: A square matrix
%        v vector
%        N Number of quadrature nodes
% OPTIONAL: these have defaults or can be put to []
%        type Type of quadrature 'Gauss' (Default) is a Gauss-Legendre
%           quadrature or 'Clenshaw-Curtis', for Clenshaw-Curtis the number
%           of nodes is selected to be the nearest power of two to the
%           supplied N
%        exptype 'direct' (default) uses expm (cubic cost), 'polynomial'
%           uses a Polynomial Krylov method or 'rational' using a
%           Rational Krylov method for the computation of the exponential
%           in the exponential sums.
%       cfpoles 15 (default) number of Carath́eodory–Fejer poles for the
%           rational approximation of the exponential function.
%       [a,b] spectral interval containing the spectrum of A
% OUTPUT: y the approximation to sinc(A)v
%         info (optional) structure containing information on the execution

%% Check of the inputs
if ~exist("type","var")
    type = 'gauss';
else
    if isempty(type)
        type = 'gauss';
    end
end

if ~exist("exptype","var")
    exptype = 'direct';
else
    if isempty(exptype)
        exptype = 'direct';
    end
end

if ~exist("cfpoles","var")
    cfpoles = 15;
else
    if isempty(cfpoles)
        cfpoles = 15;
    end
end

if strcmpi(exptype,"SAI") && (~exist("a","var") || ~exist("b","var") )
    a = eigs(A,1,"smallestabs");
    b = eigs(A,1,"largestabs");
end

%% Selecting type of exponential sum
tic;
switch upper(type)
    case "GAUSS"
        [x,w] = legpts(N);
    case "CLENSHAW-CURTIS"
        % This requires having a number of nodes that is a power of two
        t = round(log2(N));
        n = 2^t;
        epsn = arrayfun(@(k) epsilon(k,n),0:n);
        a = zeros(1,n+1);
        for j = 0:n/2
            a(2*j+1) = epsilon(2*j,n)*2/(1-4*j^2);
        end
        ahat = dct(a,n+1,'Type',1);
        w = 2/(sqrt(2*n)).*epsn.*ahat;
        x = cos( (0:n)*pi/n);
        N = length(x);
    otherwise
        warning("Unknown quadrature: defaulting to Gauss!");
        type = 'gauss';
        [x,w] = legpts(N);
end
quadpoints = toc;

%% Check for the outputs
if nargout == 2
    info.type = type;
    info.exptype = exptype;
    info.quadpoints = quadpoints;
end

%% Matrix functions computations
y = zeros(size(v));
switch upper(exptype)
    case "DIRECT"
        tic;
        % Feasible only for small matrices it has cubic cost!
        for i=1:N
            ytemp = expm(-1i*A*x(i))*v; % Cubic cost here!
            y = y + w(i)*ytemp;
        end
        computetime = toc;
    case "POLYNOMIAL"
        % Uses the expmv code from:
        % A. H. Al-Mohy and N. J. Higham, "[Computing the action of the matrix
        % exponential, with an application to exponential
        % integrators](https://doi.org/10.1137/100788860)" SIAM
        % J. Sci. Comput., 33(2):488--511, 2011.
        tic;
        M = select_taylor_degree(-A,v,[],[],'double',false);
        for i=1:N
            ytemp = expmv(1i*x(i),-A,v,M);
            y = y + w(i)*ytemp;
        end
        computetime = toc;
    case "CF"
        tic;
        [poles,~] = cf(cfpoles); % Carath́eodory–Fejer method for the type
        % (n,n) best approximation r to exp(z) on R^{-}.
        V = rat_krylov(-A, v, -poles.'); % Build rational krylov subspace
        vk = V'*v;
        Ap = -V'*A*V;
        for i=1:N
            Ak = 1i*x(i)*Ap;
            ytemp = V*(expm(Ak)*vk);
            y = y + w(i)*ytemp;
        end
        computetime = toc;
    case "RATIONAL"
        poles = genpadeexppol(cfpoles).';
        tic;
        V = rat_krylov(-A, v, poles.'); % Build rational krylov subspace
        vk = V'*v;
        Ap = -V'*A*V;
        for i=1:N
            Ak = 1i*x(i)*Ap;
            ytemp = V*(expm(Ak)*vk);
            y = y + w(i)*ytemp;
        end
        computetime = toc;
    case "SAI"
        % Shift-and-Invert method with the repetead pole from:
        % MARKO HUHTANEN Rational approximation of the unitary exponential
        % IMA Journal of Numerical Analysis (2010) 30, 512−524
        tic;
        saipoles = N + (2-(mod(N,2)));
        for i=1:N
            obj = @(ell,k,a,b) atan(sqrt( (2*k-ell)./ell )) + atan((b-a)./(2*ell)) ...
                - sqrt(ell.*(2*k-ell))./(2*k) - (b-a)./(4*k);
            ab = sort([x(i)*a,x(i)*b]);
            ell = fzero(@(ell) obj(ell,saipoles,ab(1),ab(2)),2*saipoles);
            if isnan(ell)
                ell = 2*saipoles;
            end
            xi = ones(1,saipoles)*(ell + 1i*x(i)*(a+b)/2);
            V = rat_krylov(A,v,xi);
            vk = V'*v;
            Ap = V'*A*V;
            Ak = -1i*x(i)*Ap;
            ytemp = V*(expm(Ak)*vk);
            y = y + w(i)*ytemp;
        end
        computetime = toc;
    otherwise
        tic;
        warning("Unknown exponential type defaulting to direct")
        for i=1:N
            ytemp = expm(-1i*A*x(i))*v;
            y = y + w(i)*ytemp;
        end
        computetime = toc;
end
tic;
y = real(0.5*y);
computetime = computetime+toc;
if nargout == 2
    info.computetime = computetime;
    varargout{1} = info;
end


end
%% AUXILIARY FUNCTIONS
% Auxiliary function for the computation of CLENSHAW-CURTIS quadrature
% weights and nodes
function y = epsilon(k,n)
if k == n || k == 0
    y = sqrt(2)/2;
else
    y = 1;
end
end