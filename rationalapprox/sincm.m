function [y,varargout] = sincm(A,v,varargin)
%SINCM Computes the matrix-vector product of the sinc function and a given
%vector v using Rational-Krylov methods
%   INPUT:  A sparse matrix
%           v vector
%           n (Optional: default 15) number of Laguerre Poles to use OR
%           poles (Optional: default Laguerre) vector of precomputed poles
%   OUTPUT: y approximation of sinc(A)v
%           approximation for growing number of poles

switch nargin
    case 2
        n = 20;
        poles = genlagpol(n);
    case 3
        if isscalar(varargin{1}) && round(varargin{1}) == varargin{1}
            n = varargin{1};
            poles = genlagpol(n);
        else
            poles = varargin{1};
            if iscolumn(poles)
                % rat_krylov requires poles to be a row vector!
                poles = poles.';
            end
        end
end

V = rat_krylov(A, v, poles); % Build rational krylov subspace
fAk = msincf(full(V'*A*V));
y = V*(fAk*(V'*v));

if nargout == 2
    Y = zeros(length(v),length(poles));
    for i=1:length(poles)-1
        Y(:,i) = V(:,i)*(fAk(1:i,1:i)*(V(:,i)'*v));
    end
    Y(:,end) = y;
    varargout{1} = Y;
end

end

function P = msincf(A)
%%MSINCF dense evaluation of matrix sinc function
P = funm(A,@sin);
P = A\P;
end