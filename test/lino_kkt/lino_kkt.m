function [res_ne_fastWinv,res_ne_matlabWinv,res_backslash,res_ldlsparse,t_ne,t_backslash] = lino_kkt(A,G,s,z,dims,bx,by,bz,EPS)
%Solve KKT-System using normal equations

EPS2 = 0; %Regularization on backslash-solver
nItref = 3; %number of iterative refinement steps
LINSYSACC = 1e-14;

I_k = eye(size(G,2));
I_p = eye(size(A,1));

t1 = tic;

    % LP-cone
    w_l = sqrt( s(1:dims.l) ./ z(1:dims.l) );
    V = spdiags(w_l.^2,0,dims.l,dims.l);
    scaling.l.wl = w_l;
    % W = diag(w_l);
    % Winv = diag( sqrt(z(1:dims.l)./s(1:dims.l)) );
    
    % Second-order cone [1, ?4.2] - do this construction for all cones
    for nk = 1:length(dims.q)
        
        % get variables for current cone
        coneidx = dims.l+sum(dims.q(1:nk-1))+1:dims.l+sum(dims.q(1:nk));
        zk = z(coneidx); sk = s(coneidx);
        conesize = length(coneidx);
        
        % get scaling matrix and beta for rank 1 update
        s0 = sk(1); s1 = sk(2:end);
        z0 = zk(1); z1 = zk(2:end);
        
        sres = s0^2 - s1'*s1;
        zres = z0^2 - z1'*z1;
        
        assert( sres > 0, 's not in second-order cone');
        assert( zres > 0, 'z not in second-order cone');        
        
        % scalings
        sbar = sk./sqrt(sres);
        zbar = zk./sqrt(zres);
        eta = (sres / zres)^(1/4);
        gamma = sqrt( (1+sbar'*zbar)/2 );
        wbar = 1/2/gamma*(sbar + [zbar(1); -zbar(2:end)]);
        q = wbar(2:end);
        w = q'*q;
        a = wbar(1);
        b = 1/(1+a);
        c = 1 + a + w / (1+a);
        d = 1 + 2/(1+a) + w/(1+a)^2;
        atilde = a^2 + w;
        D = eta^2*diag([atilde; ones(conesize-1,1)]);
        alpha = eta^2*d;
        beta = eta^2*c^2/d;
        
        % save stuff needed to apply scaling later in linear time
        scaling.q(nk).eta = eta;
        scaling.q(nk).a = a;
        scaling.q(nk).b = b;
        scaling.q(nk).q = q;
        scaling.q(nk).u = [c/d; q];
        scaling.q(nk).d = diag(D);
        scaling.q(nk).alpha = alpha; 
        scaling.q(nk).beta = beta;
        scaling.q(nk).W = eta*[a q'; q, eye(length(q)) + b*(q*q')];
        scaling.q(nk).V = scaling.q(nk).W*scaling.q(nk).W;
    end

% if size(z,1)==1 %LP-case
%     W = (s/z)^(1/2);
%     WinvG = G/W;
% else
%     %scaling and fast inverse of W, SOC-case
%     zbar = (1/(z(1)^2-z(2:end)'*z(2:end))^(1/2))*z;
%     sbar = (1/(s(1)^2-s(2:end)'*s(2:end))^(1/2))*s;
%     gamma = ((1+zbar'*sbar)/2)^(1/2);
%     wbar = 1/(2*gamma)*(sbar+[zbar(1); -zbar(2:end)]);
%     eta = ((s(1)^2-s(2:end)'*s(2:end))/(z(1)^2-z(2:end)'*z(2:end)))^(1/4);
%     zeta = wbar(2:end)'*G(2:end,1);
%     WinvG = 1/eta*[wbar(1)*G(1,1)-zeta; G(2:end,1)+(-G(1,1)+zeta/(1+wbar(1)))*wbar(2:end)];
%     for i=2:size(G,2)
%         zeta = wbar(2:end)'*G(2:end,i);
%         WinvG = [WinvG 1/eta*[wbar(1)*G(1,i)-zeta; G(2:end,i)+(-G(1,i)+zeta/(1+wbar(1)))*wbar(2:end)]];
%     end
%     W = eta*[wbar(1) wbar(2:end)'; wbar(2:end) eye(size(wbar(2:end),1))+wbar(2:end)*wbar(2:end)'/(1+wbar(1))];
% end
WinvG = fastWinv_kkt(G,scaling,dims);
W = diag(w_l);
for i=1:length(dims.q)
    W = blkdiag(W,scaling.q(i).W);
end

for p=1:2 %p=1: fast Winv, p==2: matlab Winv
    if p==1
        WinvG = WinvG;
    else
        WinvG = W\G;
    end
    %Factor LHS
    Y = WinvG'*WinvG + EPS*I_k;
    L_one = chol(Y,'lower');
    %Z = L_one\A'; L_one*Z = A'
    Z = L_one\A';
    M = Z'*Z + I_p*EPS;
    L = chol(M,'lower');
    if p==1
        Winvbz = fastWinv_kkt(bz,scaling,dims);%1/eta*[wbar(1)*bz(1)-zeta; bz(2:end)+(-bz(1)+zeta/(1+wbar(1)))*wbar(2:end)];
        GtWinv2bz = WinvG'*Winvbz;
    else
        Winvbz = W\bz;
        GtWinv2bz = WinvG'*Winvbz;
    end
    
    A_new = bx+GtWinv2bz;
    Z_new = conelp_forwardsub(L_one,A_new);
    RHS_ne = Z'*Z_new-by;

    %Solve
    v = conelp_forwardsub(L,RHS_ne);
    dy = conelp_backwardsub(L',v);
    dx_temp = conelp_forwardsub(L_one,bx+GtWinv2bz-A'*dy);
    dx = conelp_backwardsub(L_one',dx_temp);
    %dz = W^(-2)*(G*dx-bz) = Winv*(Winv*G*dx-Winv*bz);
    dz_temp = WinvG*dx-Winvbz;
    if p==1
        dz = fastWinv_kkt(dz_temp,scaling,dims);%1/eta*[wbar(1)*dz_temp(1)-zeta; dz_temp(2:end)+(-dz_temp(1)+zeta/(1+wbar(1)))*wbar(2:end)];
    else
        dz = W\dz_temp;
    end
    t_ne = toc(t1);
    
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

        %RHS_error_ne
        if p==1
            Winvez = fastWinv_kkt(ez,scaling,dims); %1/eta*[wbar(1)*ez(1)-zeta; ez(2:end)+(-ez(1)+zeta/(1+wbar(1)))*wbar(2:end)];
            GtWinv2ez = WinvG'*Winvez;
        else
            Winvez = W\ez;
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
        if p==1
            ddz = fastWinv_kkt(ddz_temp,scaling,dims); %1/eta*[wbar(1)*ddz_temp(1)-zeta; ddz_temp(2:end)+(-ddz_temp(1)+zeta/(1+wbar(1)))*wbar(2:end)];
        else
            ddz = W\ddz_temp;
        end
        
        dy = dy + ddy;
        dx = dx + ddx;
        dz = dz + ddz;   
    end
    if p==1
        x_ne_fastWinv= [dx;dy;dz];
    else
        x_ne_matlabWinv = [dx;dy;dz];
    end
end

%Solve using backslash
tic;
K = [EPS2*eye(size(A,2)) A' G'; A -EPS2*eye(size(A,1)) zeros(size(A,1),size(G,1)); G zeros(size(G,1),size(A,1)) -W^2];
x_backslash = K\RHS;
t_backslash = toc;

%Solve using ldlsparse
LINSOLVER = 'ldlsparse';
Gtilde = conelp_stretch(G,dims,2);
Vpattern = conelp_scaling(dims, LINSOLVER, 'pattern');
Kpattern = conelp_KKTmatrix(A,Gtilde,Vpattern,0);
P = conelp_getPerm(Kpattern~=0, 0);
[scaling,lambda,Vreg,Vtrue] = conelp_scaling(s,z,dims,LINSOLVER,EPS);
K_lifted = conelp_KKTmatrix(A, Gtilde, Vreg, EPS);
[L,D,PL,QL,P] = conelp_factor(K_lifted,P,LINSOLVER,size(G,2),size(A,1),dims,scaling);
[x,y,z,info.nitref1] = conelp_solve(L,D,P,PL,QL,bx,by,bz, A,G,Vtrue, dims, nItref,LINSOLVER,LINSYSACC);
x_ldlsparse = [x;y;z];

%Return residuals
res_ne_fastWinv = K*x_ne_fastWinv-RHS;
res_ne_matlabWinv = K*x_ne_matlabWinv-RHS;
res_backslash = K*x_backslash-RHS;
res_ldlsparse = K*x_ldlsparse-RHS;

end

