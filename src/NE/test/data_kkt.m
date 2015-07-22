function[A,G,s,z,dims,bx,by,bz] = data_kkt(nEqConstr,nStates,nLP_SOC_Constr,condA,condG,DELTA1,DELTA2)
%Generate random data for KKT-System

% rng('shuffle');

if (condA<0 || condA>0 && condA<1)
    error('no valid condition number for A')
end
if (condG<0 || condG>0 && condG<1)
    error('no valid condition number for G')
end

n = nEqConstr;
m = nStates;
k = nLP_SOC_Constr;

assert(m>=n,'nStates must be larger or equal to nEqConstraints!');
assert(n+k>=m,'nEqConstr+nLP_SOC_Constr must be larger or equal to nStates, otherwise KKT-Matrix is singular!')

% %Maximal dimensions
% maxEq = 30;
% maxStates = 300;
% maxCone = 50;

% %random dimensions
% if condA~=0 %assure that A is not a vector or a scalar
%     n = randi(maxEq)+1;
% else
%     n = randi(maxEq);
% end
% m = randi(maxStates)+(n-1); %m >= n
% if condG~=0 %assert that G is not a vector or a scalar
%     k = randi(maxCone)+m-n+1;
% else
%     k = randi(maxCone)+m-n; %n+k >= m
% end
% 
% if k == 0;
%     k = 1;
% end

%random A and G
A = rand(n,m)*200-100; %random A with values [-10,10]
G = rand(k,m)*200-100; %random G with values [-10,10]

%A and G with specific condition-number
if condA~=0
    [UA,SA,VA] = svd(A);
    SA_new = zeros(n,m);
    SA_new(1,1) = condA; %max singular-value A
    SA_new(2,2) = 1; %min singular value A
    for o=3:n
        SA_new(o,o) = rand*(condA-1)+1;
    end
    A = sparse(UA*SA_new*VA');
end
if condG~=0
    [UG,SG,VG] = svd(G);
    SG_new = zeros(k,m);
    SG_new(1,1) = condG; %max singular-value G
    SG_new(2,2) = 1; %min singular value G
    if k<=m
        l = k;
    else
        l = m;
    end
    for o=3:l
        SG_new(o,o) = rand*(condG-1)+1;
    end
    G = sparse(UG*SG_new*VG');
end

%random SOCs
dims.l = randi(k);
while (k-dims.l)==1
    dims.l = randi(k);
end
if k==dims.l
    ncones = 0;
else
    ncones = randi(floor((k-dims.l)/2));
end

s = rand(dims.l,1)*10;
z = rand(dims.l,1)*10;
dims.q = zeros(ncones,1);
for o=1:ncones
    maxDim = k-dims.l-sum(dims.q(1:o))-2*(ncones-(o));
    if o~=ncones
        dims.q(o) = randi([2,maxDim]);
    else
        dims.q(o) = maxDim;
    end
    % s or z close to the cone-border?
    s_close = randi(2,1)-1; %s_close = {0,1}
    s_temp = zeros(dims.q(o),1);
    z_temp = zeros(dims.q(o),1);
    s_temp(2:end) = rand(size(s_temp,1)-1,1)*10;
    z_temp(2:end) = rand(size(z_temp,1)-1,1)*10;
    if s_close == 1
        s_temp(1) = norm(s_temp(2:end)) + DELTA1;
        z_temp(1) = norm(z_temp(2:end)) + DELTA2;
    else
        s_temp(1) = norm(s_temp(2:end)) + DELTA1;
        z_temp(1) = norm(z_temp(2:end)) + DELTA2;
    end
    assert(s_temp(1)-norm(s_temp(2:dims.q(o))) >= 0,'s not in SOC');
    assert(z_temp(1)-norm(z_temp(2:dims.q(o))) >= 0,'z not in SOC');
    s = [s; s_temp];
    z = [z; z_temp];        
end

% soc_length = 0;
% for o=1:ncones
%     soc_length = soc_length + dims.q(o);
% end
% soc_length = dims.l + soc_length;


%RHS with random values [-100,100]
bx = rand(m,1)*200-100;
by = rand(n,1)*200-100;
bz = rand(k,1)*200-100;