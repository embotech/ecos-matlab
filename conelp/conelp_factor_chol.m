function [L] = conelp_factor_chol(A,G,V,delta)
% Transform KKT to normal equations form and do cholesky factorization

I_n = eye(size(A,2));
I_p = eye(size(A,1));
% spy(V);
% spy(G'*V^(-1)*G)
% transform to normal equations form
Y = A*(delta*I_n+G'*V^(-1)*G)^(-1)*A' + delta*I_p;

% cholesky factorization
L = chol(Y,'lower');

end

