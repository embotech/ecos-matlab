function [L,Winv] = lino_factor_slow(varargin)
%Cholesky factorization of normal equations

if nargin == 4 %init
    A = varargin{1};
    G = varargin{2};
    W = varargin{3};
    EPS = varargin{4};
elseif nargin == 5
    A = varargin{1};
    G = varargin{2};
    scaling = varargin{3};
    dims = varargin{4};
    EPS = varargin{5};
    
    W = diag(scaling.l.wl);
    for k = 1:length(dims.q)
        W = blkdiag(W,scaling.q(k).W);
    end
else
    error('lino_factor_slow needs 4 or 5 input arguments')
end

W = W; %- EPS * eye(size(W));
% W(1,1) = W(1,1) - EPS
% W = W - [EPS zeros(1,size(W,2)-1);zeros(size(W,1)-1,1) eye(size(W)-1)*EPS];

I_n = eye(size(A,2));
I_p = eye(size(A,1));

Winv = W^(-1);
V = Winv*G;
Y = V'*V + I_n*EPS;
L_one = chol(Y,'lower');
%Z = L_one\A'; L_one*Z = A'
Z = conelp_forwardsub(L_one,A(1,:)');
for i=2:size(A,1)
    z = conelp_forwardsub(L_one,A(i,:)');
    Z = [Z z];
end

M = Z'*Z + I_p*EPS;
L = chol(M,'lower');