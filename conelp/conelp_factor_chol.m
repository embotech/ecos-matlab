function [L] = conelp_factor_chol(A,G,V,delta)
% Transform KKT to normal quations form and do cholesky factorization

if ( nargin == 3 ), delta = 0;
end

I_n = eye(size(A,2));
I_p = eye(size(A,1));

% transform to normal equation form
Y = A*(delta*I_n+G'*V^(-1)*G)^(-1)*A'+delta*I_p;

% cholesky factorization
L = chol(Y,'lower');

end

