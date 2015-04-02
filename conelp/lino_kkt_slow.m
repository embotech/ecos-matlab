function [x,y,z,nitref] = lino_kkt_slow(varargin)
% Solve KKT system using the cholesky factorization

if nargin == 11 %init
    L = varargin{1};
    bx = varargin{2};
    by = varargin{3};
    bz = varargin{4};
    A = varargin{5};
    G = varargin{6};
    W = varargin{7};
    nItref = varargin{8};
    LINSYSACC = varargin{9};
    delta = varargin{10};
    Winv = varargin{11};
elseif nargin == 12
    L = varargin{1};
    bx = varargin{2};
    by = varargin{3};
    bz = varargin{4};
    A = varargin{5};
    G = varargin{6};
    scaling = varargin{7};
    dims = varargin{8};
    nItref = varargin{9};
    LINSYSACC = varargin{10};
    delta = varargin{11};
    Winv = varargin{12};
    
    W = diag(scaling.l.wl);
    for k = 1:length(dims.q)
         W = blkdiag(W, scaling.q(k).W);
    end
else
    error('lino_kkt_slow needs 11 or 12 input arguments')
end

n = length(bx);
p = length(by);
I_n = eye(n);
RHS = [bx;by;bz];

Winv2 = Winv'*Winv;
T = G'*Winv2;
Tbar = (delta*I_n+T*G)^(-1);

% Prepare RHS of normal equations
RHS_ne = A*Tbar*(bx+T*bz)-by;

% initial solve
v = conelp_forwardsub(L,RHS_ne);
dy = conelp_backwardsub(L',v);
dx = Tbar*(bx+T*bz-A'*dy);
dz = Winv2*(G*dx-bz);

% iterative refinement
bnorm = 1+norm(RHS,inf);
for i = 1:nItref

    % errors
    ex = bx - A'*dy - G'*dz;
    ey = by - A*dx;
    ez = bz - G*dx + W^2*dz;
    
        if(norm(ex,inf)/bnorm < LINSYSACC && ...
           norm(ey,inf)/bnorm < LINSYSACC && ...
           norm(ez,inf)/bnorm < LINSYSACC ), break; end
    
    % solve for correction
    RHS_error_ne = A*Tbar*(ex+T*ez)-ey;
    
    v = conelp_forwardsub(L,RHS_error_ne);
    ddy = conelp_backwardsub(L',v);
    ddx = Tbar*(ex+T*ez-A'*ddy);
    ddz = Winv2*(G*ddx-ez);

    dy = dy + ddy;
    dx = dx + ddx;
    dz = dz + ddz;   
end
    
   % copy out variables
    x = dx;
    y = dy;
    z = dz;
    nitref = i;
