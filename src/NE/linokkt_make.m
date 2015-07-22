function ne_make
% Compile linokkt_mex

d = '' ;
if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims' ;
end
eval (sprintf ('mex %s -I../../ecos/include -I../../ecos/external/CHOLMOD/Include -I../../ecos/AMD/Include -I../../ecos/external/LDL/Include -I../../ecos/external/SuiteSparse_config -L../../ecos/external/CHOLMOD/Lib -L../../ecos/external/AMD/Lib -L../../ecos/external/COLAMD/Lib -L../../ecos/external/SuiteSparse_config -L../../ecos/external/SuiteSparse_config/xerbla -L../../ecos/external/CAMD/Lib -L../../ecos/external/CCOLAMD/Lib -lcholmod -lamd -lcolamd -lsuitesparseconfig -lcerbla -lcamd -lccolamd -lm -lgoto2 -llapack linokkt_mex.c ../../ecos/src/lino_kkt.c ../../ecos/src/cone.c ../../ecos/src/splamm.c ../../ecos/src/spla.c ../../ecos/src/timer.c', d));
fprintf('linokkt_mex successfully compiled.\n');

