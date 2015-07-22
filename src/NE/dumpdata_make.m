function dumpdata_make
% Compile dumpdata_mex

d = '' ;
if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims' ;
end
eval (sprintf ('mex %s -I../../ecos/include -I../../ecos/external/SuiteSparse_config test/dumpdata_mex.c ../../ecos/src/splamm.c ../../ecos/src/spla.c', d));
fprintf('dumpdata_mex successfully compiled.\n');

