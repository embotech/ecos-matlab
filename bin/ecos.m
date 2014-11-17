function [x,y,info,s,z] = ecos(c,G,h,dims,varargin)
% ECOS - Embedded COnic Solver.
%
% Self-dual homogeneous embedding interior point implementation for optimization
% over linear or second-order cones. ECOS does not support semi-definite
% cones (feel free to contact the developers if you wish to add SDP support 
% to ECOS.
%
%   [x,y,info,s,z] = ECOS(c,G,h,dims) Solves a pair of primal and dual
%   cone programs
% 
%        minimize    c'*x
%        subject to  G*x + s = h
%                    s >= 0
%
%        maximize    -h'*z
%        subject to  G'*z + c = 0
%                    z >= 0.
%
%   The inequalities are with respect to a cone K defined as the Cartesian
%   product of N+1 cones:
% 
%        K = K_0 x K_1 x .... x K_N.
% 
%     The first cone, K_0, is the nonnegative orthant of dimension dims.l.
%     The next N cones are second order cones of dimension dims.q(1), ...,
%     dims.q(N), where a second order cone of dimension m is defined as
% 
%         { (u0, u1) in R x R^{m-1} | u0 >= ||u1||_2 }.
% 
%     INPUT arguments:
% 
%         c is a dense column vector of size n
% 
%         dims is a struct with the dimensions of the components of cone K.
%         It has two fields.
%         - dims.l, the dimension of the nonnegative orthant C_0, with l>=0.
%         - dims.q, a row vector of N integers with the dimensions of the 
%           second order cones K_1, ..., K_N. (N >= 0 and q(i) >= 3.)
% 
%         G is a sparse matrix of size (m,n), where
% 
%             m = dims.l + dims.q(1) + ... + dims.q(N).
%               = dims.l + sum(dims.q)
% 
%         Each column of G describes a vector
% 
%             v = ( v_0, v_1, ..., v_N )
% 
%         in V = R^dims.l x R^dims.q(1) x ... x R^dims.q(N)
%         stored as a column vector
% 
%             [ v_0; v_1; ...; v_N ].
% 
%         h is a dense column vector of size m, representing a vector in V,
%         in the same format as the columns of G.
%         
%
%
%   [x,y,info,s,z] = ECOS(c,G,h,dims,A,b) Solves a pair of primal and 
%   dual cone programs
% 
%        minimize    c'*x
%        subject to  G*x + s = h
%                    A*x = b
%                    s >= 0
%
%        maximize    -h'*z - b'*y
%        subject to  G'*z + A'*y + c = 0
%                    z >= 0.
%
%      where c,G,h,dims are defined as above, and A is a sparse matrix of 
%      size (p,n), and b is a dense matrix of size (p,1).
%       
%      It is assumed that rank(A) = p and rank([A; G]) = n.
%
% 
%   [x,y,info,s,z] = ECOS(c,G,h,dims,otps) and 
%   [x,y,info,s,z] = ECOS(c,G,h,dims,A,b,otps) are as above, with the struct 
%   otps used to control settings of the solver. The following fields can
%   be present (use ECOSOPTIMSET to obtain a default initialization):
%
%      .verbose - whether to inform on progess (default: true)
%      .feastol - stopping tolerance on infeasibilities (default: 1e-5)
%      .abstol  - stopping tolerance (default: 1e-6)
%      .reltol  - rel. stopping toletance (default: 1e-6)
%      .maxit   - maximum number of iterations (default: 30)
% 
%
% Details on ECOS can be found at http://embotech.com/ECOS and in the paper 
%    Alexander Domahidi, Eric Chu, Stephen Boyd. "ECOS: An Embedded
%    Conic Solver." In proceedings of European Control Conference (ECC), 
%    pp. 3071-3076, Zurich, Switzerland, July 2013."
%
%
% (c) A. Domahidi, ETH Zurich & embotech GmbH, Zurich, Switzerland, 2012-14.
%
%
% COPYING: ECOS is under GPLv3. For commercial licenses and professional
% support send an email to ecos@embotech.com.
%
% See also ECOSQP ECOSOPTIMSET ECOS_LICENSE

if (length(varargin) == 1)
    otps = varargin{1};
elseif (length(varargin) == 2)
    A = varargin{1};
    b = varargin{2};
elseif (length(varargin) == 3)
    A = varargin{1};
    b = varargin{2};
    otps = varargin{3};
end

if exist('otps','var')
    if (isfield(otps, 'bool_vars_idx') && ...
        (min(otps.bool_vars_idx(:)) < 1 || max(otps.bool_vars_idx(:)) > max(size(c))) )
        error('ecos:InvalidInput', 'otps.bool_vars_idx must be in [1,length(c)]');
    end

    if (isfield(otps, 'int_vars_idx') && ...
            (min(otps.int_vars_idx(:)) < 1 || max(otps.int_vars_idx(:)) > max(size(c))) )
        error('ecos:InvalidInput', 'otps.int_vars_idx must be in [1,length(c)]');
    end
end

if (nargin == 4)
    [x,y,info,s,z] = ecos_c(c,G,h,dims);
elseif (nargin == 5)
    [x,y,info,s,z] = ecos_c(c,G,h,dims,otps);
elseif (nargin == 6)
    [x,y,info,s,z] = ecos_c(c,G,h,dims,A,b);
elseif (nargin == 7)
    [x,y,info,s,z] = ecos_c(c,G,h,dims,A,b,otps);
else
    error('ecos:InvalidInput', 'Invalid call to ecos, please type "help ecos" for correct usage');
end

return