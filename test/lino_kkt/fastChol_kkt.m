function [L_one,D] = fastChol_kkt(G,dims,scaling,EPS)
%Computes a cholesky factorization of G'*W(-2)*G using down- and updates

S = zeros(size(G,2));

%LP-cones
for i=1:dims.l
    S = S + G(i,:)'*G(i,:)/scaling.l.wl(i)^2;
end

%SOCs
for i=1:length(dims.q)
    %new scaling representation
    wold0 = scaling.q(i).a;
    wold1 = scaling.q(i).q;
    w0 = (1/2*(wold0^2+wold1'*wold1+1))^(1/2);
    w1 = (-1)/(2*w0)*(wold0*wold1+(eye(size(wold1,1))+1/(1+wold0)*wold1*wold1')*wold1);
    w = [w0;w1];
    beta = scaling.q(i).eta;
    %compute new S
    G_blockrow = G(dims.l+sum(dims.q(1:i-1))+1:dims.l+sum(dims.q(1:i)),:);
    S = S + G_blockrow'*G_blockrow/beta^2;
    %compute scalars&vectors for down- and upgrade
    a = 2/beta^2;
    b = a;
    update.q(i).u = sqrt(a)*G_blockrow'*w;
    update.q(i).v = sqrt(b)*G_blockrow'*[1; zeros(dims.q(i)-1,1)];    
end

% for i=1:length(dims.q)
%     S = S + update.q(i).u*update.q(i).u'-update.q(i).v*update.q(i).v';
% end

%S = S + eye(size(S))*EPS;
LD = ldlchol(sparse(S));
for i=1:length(dims.q)
    LD = ldlupdate(LD,sparse(update.q(i).u));
    LD = ldlupdate(LD,sparse(update.q(i).v),'-');
end
[L_one,D] = ldlsplit(LD);
L_one = L_one*sqrt(D);