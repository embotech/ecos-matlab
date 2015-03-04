function [x,y,z,nitref] = conelp_solve_chol(L,bx,by,bz,A,G,V,dims,nItref,LINSYSACC,delta)
% Solve KKT system using the cholesky factorization

if ( nargin == 10 ), delta = 0;
end

n = length(bx);
p = length(by);
I_n = eye(n);
RHS = [bx;by;bz];

T = G'*V^(-1);
Tbar = (delta*I_n+T*G)^(-1);

% rows to be added
Nstretch = 2;

% Prepare RHS of normal equation
bztilde = conelp_stretch(bz, dims, Nstretch);
RHS_ne = A*Tbar*(bx+T*bztilde)-by;

% initial solve
v = conelp_forwardsub(L,RHS_ne);
dy = conelp_backwardsub(L',v);
dx = Tbar*(bx+T*bztilde-A'*dy);
dz = V^(-1)*(G*dx-bztilde);

% iterative refinement
bnorm = 1+norm(RHS,inf);
for i = 1:nItref
    
    dz_unstretched = conelp_unstretch(dz,dims, Nstretch);
    G_unstretched = conelp_unstretch(G,dims,Nstretch);
    
    % errors
    ex = bx - A'*dy - G_unstretched'*dz_unstretched;
    ey = by - A*dx;
    eztilde = conelp_stretch(bz - G_unstretched*dx,dims,Nstretch) + V*dz;
    
        if(norm(ex,inf)/bnorm < LINSYSACC && ...
           norm(ey,inf)/bnorm < LINSYSACC && ...
           norm(eztilde,inf)/bnorm < LINSYSACC ), break; end
    
    % solve for correction
    RHS_error_ne = A*Tbar*(ex+T*eztilde)-ey;
    
    v = conelp_forwardsub(L,RHS_error_ne);
    ddy = conelp_backwardsub(L',v);
    ddx = Tbar*(ex+T*eztilde-A'*ddy);
    ddz = V^(-1)*(G*ddx-eztilde);

    dy = dy + ddy;
    dx = dx + ddx;
    dz = dz + ddz;   
end
    
   % copy out variables
    x = dx;
    y = dy;
    z = conelp_unstretch(dz,dims,Nstretch);
    nitref = i;

