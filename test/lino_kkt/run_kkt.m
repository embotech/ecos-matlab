%Run KKT-solvers
% clear
run('data.m')
[x,y,z,x_backslash,t_ne,t_backslash] = lino_kkt(A,G,s,z,bx,by,bz);
x_ne = [x;y;z];
%t_ne
%t_backslash
DELTA = 1e-5;
n_errors = 0;
for i=1:size(x_ne,1)
    if (x_ne(i)/x_backslash(i)<=1-DELTA) || (x_ne(i)/x_backslash(i)>=1+DELTA)
        n_errors = n_errors + 1;
    end
end
n_errors