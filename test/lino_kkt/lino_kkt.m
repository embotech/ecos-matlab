function [dx,dy,dz,x_backslash,t_ne,t_backslash] = lino_kkt(A,G,s,z,bx,by,bz)
%Solve KKT-System using normal equations

EPS = 5e-7; %Regularization on normal equations solver
EPS2 = 0; %Regularization on backslash-solver
nItref = 3; %number of iterative refinement steps
LINSYSACC = 1e-10;

I_k = eye(size(G,2));
I_p = eye(size(A,1));

t1 = tic;

if size(z,1)==1 %LP-case
    W = (s/z)^(1/2);
    WinvG = G/W;
else
    %scaling and fast inverse of W, SOCP-case
    zbar = (1/(z(1)^2-z(2:end)'*z(2:end))^(1/2))*z;
    sbar = (1/(s(1)^2-s(2:end)'*s(2:end))^(1/2))*s;
    gamma = ((1+zbar'*sbar)/2)^(1/2);
    wbar = 1/(2*gamma)*(sbar+[zbar(1); -zbar(2:end)]);
    eta = ((s(1)^2-s(2:end)'*s(2:end))/(z(1)^2-z(2:end)'*z(2:end)))^(1/4);
    zeta = wbar(2:end)'*G(2:end,1);
    WinvG = 1/eta*[wbar(1)*G(1,1)-zeta; G(2:end,1)+(-G(1,1)+zeta/(1+wbar(1)))*wbar(2:end)];
    for i=2:size(G,2)
        zeta = wbar(2:end)'*G(2:end,i);
        WinvG = [WinvG 1/eta*[wbar(1)*G(1,i)-zeta; G(2:end,i)+(-G(1,i)+zeta/(1+wbar(1)))*wbar(2:end)]];
    end
    W = eta*[wbar(1) wbar(2:end)'; wbar(2:end) eye(size(wbar(2:end),1))+wbar(2:end)*wbar(2:end)'/(1+wbar(1))];
end

%Factor LHS
Y = WinvG'*WinvG + EPS*I_k;
t2 = tic;
L_one = chol(Y,'lower');
t_chol1 = toc(t2)
%Z = L_one\A'; L_one*Z = A'
t3 = tic;
Z = zeros(size(L_one,1),size(A,1));
for i=1:size(A,1)
    z_temp = conelp_forwardsub(L_one,A(i,:)');
    Z(:,i) = z_temp;
end
t_Z = toc(t3)
M = Z'*Z + I_p*EPS;
L = chol(M,'lower');
%RHS
if size(z,1)==1 %LP-case
    Winvbz = bz/W;
    GtWinv2bz = G'*bz/W^2;
else % SOCP-case
    zeta = wbar(2:end)'*bz(2:end);
    Winvbz = 1/eta*[wbar(1)*bz(1)-zeta; bz(2:end)+(-bz(1)+zeta/(1+wbar(1)))*wbar(2:end)];
    GtWinv2bz = WinvG'*Winvbz;
end
A_new = bx+GtWinv2bz;
Z_new = conelp_forwardsub(L_one,A_new);
RHS_ne = Z'*Z_new-by;

%Solve
t4 = tic;
v = conelp_forwardsub(L,RHS_ne);
dy = conelp_backwardsub(L',v);
dx_temp = conelp_forwardsub(L_one,bx+GtWinv2bz-A'*dy);
dx = conelp_backwardsub(L_one',dx_temp);
%dz = W^(-2)*(G*dx-bz) = Winv*(Winv*G*dx-Winv*bz);
dz_temp = WinvG*dx-Winvbz;
if size(z,1)==1
    dz = (G*dx-bz)/W^2;
else
    zeta = wbar(2:end)'*dz_temp(2:end);
    dz = 1/eta*[wbar(1)*dz_temp(1)-zeta; dz_temp(2:end)+(-dz_temp(1)+zeta/(1+wbar(1)))*wbar(2:end)];
end
t_solve = toc(t4)
t_ne = toc(t1)

% iterative refinement
RHS = [bx;by;bz];
bnorm = 1+norm(RHS,inf);
for i = 1:nItref

    % errors
    ex = bx - EPS2*eye(size(A,2))*dx - A'*dy - G'*dz;
    ey = by - A*dx + EPS2*eye(size(A,1))*dy;
    ez = bz - G*dx + W^2*dz;
    
        if(norm(ex,inf)/bnorm < LINSYSACC && ...
           norm(ey,inf)/bnorm < LINSYSACC && ...
           norm(ez,inf)/bnorm < LINSYSACC ), break; end
    
    % solve for correction
    
    %RHS error ne
    if size(z,1)==1 %LP-case
        Winvez = ez/W;
        GtWinv2ez = G'*ez/W^2;
    else % SOCP-case
        zeta = wbar(2:end)'*ez(2:end);
        Winvez = 1/eta*[wbar(1)*ez(1)-zeta; ez(2:end)+(-ez(1)+zeta/(1+wbar(1)))*wbar(2:end)];
        GtWinv2ez = WinvG'*Winvez;
    end
    A_new = ex+GtWinv2ez;
    Z_new = conelp_forwardsub(L_one,A_new);
    RHS_error_ne = Z'*Z_new-ey;
    
    v = conelp_forwardsub(L,RHS_error_ne);
    ddy = conelp_backwardsub(L',v);
    ddx_temp = conelp_forwardsub(L_one,ex+GtWinv2ez-A'*ddy);
    ddx = conelp_backwardsub(L_one',ddx_temp);
    ddz_temp = WinvG*ddx-Winvez;
    if size(z,1)==1
        ddz = (G*ddx-ez)/W^2;
    else
        zeta = wbar(2:end)'*ddz_temp(2:end);
        ddz = 1/eta*[wbar(1)*ddz_temp(1)-zeta; ddz_temp(2:end)+(-ddz_temp(1)+zeta/(1+wbar(1)))*wbar(2:end)];
    end

    dy = dy + ddy;
    dx = dx + ddx;
    dz = dz + ddz;   
end

%Solve using backslash
tic;
if size(z,1)==1
    W = (s/z)^(1/2);
else
    W = eta*[wbar(1) wbar(2:end)'; wbar(2:end) eye(size(wbar(2:end),1))+wbar(2:end)*wbar(2:end)'/(1+wbar(1))];
end
Abig = [EPS2*eye(size(A,2)) A' G'; A -EPS2*eye(size(A,1)) zeros(size(A,1),size(G,1)); G zeros(size(G,1),size(A,1)) -W^2];
x_backslash = Abig\RHS;
t_backslash = toc;

end

