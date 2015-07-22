%Run KKT-solvers
clear;
format long e;
N = 1; % Number of problems to solve
NORM = Inf; % Norm of solution to analyze
EPS = 5e-9; % Regularization on normal equations solver

%Problem data dimensions: nStates>=nEqConstr and
%nEqConstr+nLP_SOC_Constry>=nStates
nEqConstr = 70; %Number of equality constraints (size(A,1))
nStates = 100; %Number of states (size(A,2), size(G,2))
nLP_SOC_Constr = 80; %Number of LP and SOC constraints (size(G,1))

assert(nStates>=nEqConstr,'nStates must be larger or equal to nEqConstraints!');
assert(nEqConstr+nLP_SOC_Constr>=nStates,'nEqConstr+nLP_SOC_Constr must be larger or equal to nStates, otherwise KKT-Matrix is singular!')

%Problem data conditioning
DELTA1 = 1e-7; % how much in the cone are s and z (or how regular is W)
DELTA2 = 1;
condA = 10; % Condition number of A, enter 0 if no specific condition number desired
condG = 10; % Condition number of G, enter 0 if no specific condition number desired

for ntest=1:N
    [A,G,s,z,dims,bx,by,bz] = data_kkt(nEqConstr,nStates,nLP_SOC_Constr,condA,condG,DELTA1,DELTA2);
    A = 0;
    G = 0;
    while(sprank([A;G]) < nStates)
        A = 0;
        G = 0;
        while(sprank(A)~=nEqConstr)
            A = sprand(nEqConstr,nStates,0.2);
        end
        while(sprank(G)~=nLP_SOC_Constr)
            G = sprand(nLP_SOC_Constr,nStates,0.2);
        end
    end
    while(sprank([A;G]) < nStates)
        A = sprand(nEqConstr,nStates,0.2);
        G = sprand(nLP_SOC_Constr,nStates,0.2);
    end
    nlp = dims.l; %number of lp-cones
    if(dims.q)
        nsoc = length(dims.q); %number of socs
        socdims = int64(dims.q);
        % dumpdata_mex(A,G,s,z,nlp,nsoc,socdims,bx,by,bz,EPS);
        [x_backslash, K, RHS] = test_kkt_matlab(A,G,s,z,dims,EPS,bx,by,bz);
        [dx,dy,dz] = linokkt_mex(A,G,s,z,nlp,nsoc,socdims,EPS,bx,by,bz);
    else
        disp('LP-only!');
        % dumpdata_mex(A,G,s,z,nlp,0,0,bx,by,bz,EPS);
        [x_backslash, K, RHS] = test_kkt_matlab(A,G,s,z,dims,EPS,bx,by,bz);
        [dx,dy,dz] = linokkt_mex(A,G,s,z,nlp,0,0,EPS,bx,by,bz);
    end
    x = [dx;dy;dz];
    
    % Compute residuals and print worst
    cnorm = norm(K*x-RHS,Inf)
    mnorm = norm(K*x_backslash-RHS,Inf)
    
end