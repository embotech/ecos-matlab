
function [x_backslash, K, RHS] = test_factor_matlab(A,G,s,z,dims,EPS,bx,by,bz)

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
        s02z02stsztz = [s0^2 z0^2 s1'*s1 z1'*z1];
        
        assert( sres > 0, 's not in second-order cone');
        assert( zres > 0, 'z not in second-order cone');        
        % w0 new calculation
        zktsk = 0;
        for i=1:conesize
            zktsk = zktsk + zk(i)'*sk(i);
        end
        w0test = (s0*zres^(1/2)+z0*sres^(1/2))/(2^(1/2)*(sres^(1/2)*zres^(3/2)+zktsk*zres)^(1/2));
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
        atimeseta = a*eta;
        
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
    
    W = diag(w_l);
    for i=1:length(dims.q)
        W = blkdiag(W,scaling.q(i).W);
    end    

%     S = zeros(size(G,2));
% 
%     %LP-cones
%     for i=1:dims.l
%         S = S + G(i,:)'*G(i,:)/scaling.l.wl(i)^2;
%     end
% 
%     %SOCs
%     for i=1:length(dims.q)
%         %new scaling representation
%         wold0 = scaling.q(i).a;
%         wold1 = scaling.q(i).q;
%         w0 = (1/2*(wold0^2+wold1'*wold1+1))^(1/2);
%         w1 = (-1)/(2*w0)*(wold0*wold1+(eye(size(wold1,1))+1/(1+wold0)*wold1*wold1')*wold1);
%         w = [w0;w1];
%         beta = scaling.q(i).eta;
%         %compute new S
%         G_blockrow = G(dims.l+sum(dims.q(1:i-1))+1:dims.l+sum(dims.q(1:i)),:);
%         S = S + G_blockrow'*G_blockrow/beta^2;
%         %compute scalars&vectors for down- and upgrade
%         a = 2/beta^2;
%         b = a;
%         update.q(i).u = sqrt(a)*G_blockrow'*w;
%         sqrt(a)*G_blockrow'*w;
%         update.q(i).v = sqrt(b)*G_blockrow'*[1; zeros(dims.q(i)-1,1)];
%         sqrt(b)*G_blockrow'*[1; zeros(dims.q(i)-1,1)];
%     end
%     D = S + EPS*eye(size(S));
%     Ls = chol(D,'lower');
%     Dnew = Ls*Ls';
%     GtWinv2G = EPS*eye(size(G,2))+G'*W^(-2)*G;
%     nnzGtWinv2G = nnz(GtWinv2G);
%     Z = Ls\A';
%     M = A*(EPS*eye(size(G,2))+G'*W^(-2)*G)^(-1)*A'+EPS*eye(size(A,1));
%     Ms = chol(M,'lower');
    
    RHS = [bx;by;bz];
    EPS2 = 0;
    K = [EPS2*eye(size(A,2)) A' G'; A -EPS2*eye(size(A,1)) zeros(size(A,1),size(G,1)); G zeros(size(G,1),size(A,1)) -W^2];
    x_backslash = K\RHS;
    
end