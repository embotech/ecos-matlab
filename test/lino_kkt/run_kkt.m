%Run KKT-solvers

clear;
N = 10; % Number of problems to solver
NORM = Inf; % Norm of solution to analyze
EPS = 1e-7; % Regularization on normal equations solver

%Problem data dimensions: nStates>=nEqConstr and
%nEqConstr+nLP_SOC_Constry>=nStates
nEqConstr = 250; %Number of equality constraints (size(A,1))
nStates = 500; %Number of states (size(A,2), size(G,2))
nLP_SOC_Constr = 350; %Number of LP and SOC constraints (size(G,1))

assert(nStates>=nEqConstr,'nStates must be larger or equal to nEqConstraints!');
assert(nEqConstr+nLP_SOC_Constr>=nStates,'nEqConstr+nLP_SOC_Constr must be larger or equal to nStates, otherwise KKT-Matrix is singular!')

%Problem data conditioning
DELTA = 1e-12; % how much in the cone are s and z (or how regular is W)
condA = 100; % Condition number of A
condG = 100; % Condition number of G

res_ne_fastWinv_max = zeros(N,1);
res_ne_matlabWinv_max = zeros(N,1);
res_backslash_max = zeros(N,1);
res_ldlsparse_max = zeros(N,1);
for ntest=1:N
    [A,G,s,z,dims,bx,by,bz] = data_kkt(nEqConstr,nStates,nLP_SOC_Constr,condA,condG,DELTA);
    [res_ne_fastWinv,res_ne_matlabWinv,res_backslash,res_ldlsparse,t_ne,t_backslash] = lino_kkt(A,G,s,z,dims,bx,by,bz,EPS);
    res_ne_fastWinv_max(ntest) = norm(res_ne_fastWinv,NORM);
    res_ne_matlabWinv_max(ntest) = norm(res_ne_matlabWinv,NORM);
    res_backslash_max(ntest) = norm(res_backslash,NORM);
    res_ldlsparse_max(ntest) = norm(res_ldlsparse,NORM);
    ntest
end

str_plot = sprintf('N = %d, nEqConstr = %d, nStates = %d, nLPSOCConstr = %d, DELTA = %d, condA = %d, condG = %d',N,nEqConstr,nStates,nLP_SOC_Constr,DELTA,condA,condG);

fastWinv_mean_error = mean(res_ne_fastWinv_max);
matlabWinv_mean_error = mean(res_ne_matlabWinv_max);
backslash_mean_error = mean(res_backslash_max);
ldlsparse_mean_error = mean(res_ldlsparse_max);

str_mean1 = sprintf('Normal equations fast Winv, mean max error = %d',fastWinv_mean_error);
str_mean2 = sprintf('Normal equations matlab Winv, mean max error = %d',matlabWinv_mean_error);
str_mean3 = sprintf('Backslash, mean max error = %d',backslash_mean_error);
str_mean4 = sprintf('Ldlsparse, mean max error = %d',ldlsparse_mean_error);

figure
boxplot([res_ne_fastWinv_max,res_ne_matlabWinv_max,res_backslash_max,res_ldlsparse_max],'labels',{str_mean1,str_mean2,str_mean3,str_mean4})
title(str_plot)
xlabel('Solver')
ylabel('Residuals')

figure
boxplot([res_ne_fastWinv_max,res_backslash_max],'labels',{str_mean1,str_mean3})
title(str_plot)
xlabel('Solver')
ylabel('Residuals')