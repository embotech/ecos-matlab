%Run KKT-solvers
N = 100;
res_ne_fastWinv_max = zeros(N,1);
res_ne_matlabWinv_max = zeros(N,1);
res_backslash_max = zeros(N,1);
res_ldlsparse_max = zeros(N,1);
for ntest=1:N
    run('data.m')
    [res_ne_fastWinv,res_ne_matlabWinv,res_backslash,res_ldlsparse,t_ne,t_backslash] = lino_kkt(A,G,s,z,bx,by,bz);
    res_ne_fastWinv_max(ntest) = norm(res_ne_fastWinv,Inf);
    res_ne_matlabWinv_max(ntest) = norm(res_ne_matlabWinv,Inf);
    res_backslash_max(ntest) = norm(res_backslash,Inf);
    res_ldlsparse_max(ntest) = norm(res_ldlsparse,Inf);
end

boxplot([res_ne_fastWinv_max,res_ne_matlabWinv_max,res_backslash_max,res_ldlsparse_max],'labels',{'Normal equations fast Winv','Normal equations matlab Winv','Backslash','ldlsparse'})
title('Solver accuracy')
xlabel('Solver')
ylabel('Residuals')