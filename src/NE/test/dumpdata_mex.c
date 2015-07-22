/* Dump data to file data.h */ 

#include "mex.h"
#include "matrix.h"
#include "splamm.h"
#include <stdio.h>

/* The mex-function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if(nrhs == 11){
        /* Get inputs */
        mxArray* Amx = (mxArray*)prhs[0];
        mxArray* Gmx = (mxArray*)prhs[1];
        pfloat* s = (pfloat*)mxGetPr(prhs[2]);
        pfloat* z = (pfloat*)mxGetPr(prhs[3]);
        idxint nlp = (idxint)mxGetScalar(prhs[4]);
        idxint nsoc = (idxint)mxGetScalar(prhs[5]);
        idxint* socdims = mxGetData(prhs[6]);
        pfloat* bx = mxGetData(prhs[7]);
        pfloat* by = mxGetData(prhs[8]);
        pfloat* bz = mxGetData(prhs[9]);
        pfloat delta = (pfloat)mxGetScalar(prhs[10]);
        
        spmat* A = createSparseMatrix((mwIndex)mxGetM(Amx),(mwIndex)mxGetN(Amx),(mwIndex)mxGetNzmax(Amx),(mwIndex*)mxGetJc(Amx),(mwIndex*)mxGetIr(Amx),(pfloat*)mxGetPr(Amx));
        spmat* G = createSparseMatrix((mwIndex)mxGetM(Gmx),(mwIndex)mxGetN(Gmx),(mwIndex)mxGetNzmax(Gmx),(mwIndex*)mxGetJc(Gmx),(mwIndex*)mxGetIr(Gmx),(pfloat*)mxGetPr(Gmx));
        
        idxint i, j;
        FILE* fp = fopen("data_valgrind.h","w");
        
        /* Print header to file */
        fprintf(fp,"/* Data for valgrind */\n");
        
        /* Include */
        fprintf(fp, "#include \"lino_kkt.c\"\n\n");
        
        fprintf(fp,"/* nlp, nsoc, socdims and delta */\n");
        fprintf(fp,"idxint nlp = %i;\n",nlp);
        fprintf(fp,"idxint nsoc = %i;\n",nsoc);
        fprintf(fp,"idxint socdims[%i] = {",nsoc);
        if(!nsoc){
            fprintf(fp,"0};\n");
        }
        else{
            for(i = 0; i < nsoc; i++){
                if(i != nsoc-1)
                    fprintf(fp,"%i,",socdims[i]);
                else
                    fprintf(fp,"%i};\n",socdims[i]);
            }
        }
        fprintf(fp,"pfloat delta = 1e-7;\n"); 
    
        fprintf(fp,"/* Sparse matrix A */\n");
        fprintf(fp,"idxint Am = %i;\n",A->m);
        fprintf(fp,"idxint An = %i;\n",A->n);
        fprintf(fp,"idxint Annz = %i;\n",A->nnz);
        fprintf(fp,"idxint ajc[%i] = {",A->n+1);
        for(i = 0; i <= A->n; i++){
            if(i != A->n)
                fprintf(fp,"%i,",A->jc[i]);
            else
                fprintf(fp,"%i};\n",A->jc[i]);
        }
        
        fprintf(fp,"idxint air[%i] = {",A->nnz);
        for(i = 0; i < A->nnz; i++){
            if(i != A->nnz-1)
                fprintf(fp,"%i,",A->ir[i]);
            else
                fprintf(fp,"%i};\n",A->ir[i]);
        }
        
        fprintf(fp,"pfloat apr[%i] = {",A->nnz);
        for(i = 0; i < A->nnz; i++){
            if(i != A->nnz-1)
                fprintf(fp,"%f,",A->pr[i]);
            else
                fprintf(fp,"%f};\n",A->pr[i]);
        }
        
        fprintf(fp,"/* Sparse matrix G */\n");
        fprintf(fp,"idxint Gm = %i;\n",G->m);
        fprintf(fp,"idxint Gn = %i;\n",G->n);
        fprintf(fp,"idxint Gnnz = %i;\n",G->nnz);
        fprintf(fp,"idxint gjc[%i] = {",G->n+1);
        for(i = 0; i <= G->n; i++){
            if(i != G->n)
                fprintf(fp,"%i,",G->jc[i]);
            else
                fprintf(fp,"%i};\n",G->jc[i]);
        }
        
        fprintf(fp,"idxint gir[%i] = {",G->nnz);
        for(i = 0; i < G->nnz; i++){
            if(i != G->nnz-1)
                fprintf(fp,"%i,",G->ir[i]);
            else
                fprintf(fp,"%i};\n",G->ir[i]);
        }
        
        fprintf(fp,"pfloat gpr[%i] = {",G->nnz);
        for(i = 0; i < G->nnz; i++){
            if(i != G->nnz-1)
                fprintf(fp,"%f,",G->pr[i]);
            else
                fprintf(fp,"%f};\n",G->pr[i]);
        }
        
        fprintf(fp,"/* s */\n");
        fprintf(fp,"pfloat s[%i] = {",G->m);
        for(i = 0; i < G->m; i++){
            if(i != G->m-1)
                fprintf(fp,"%f,",s[i]);
            else
                fprintf(fp,"%f};\n",s[i]);
        }
        
        fprintf(fp,"/* z */\n");
        fprintf(fp,"pfloat z[%i] = {",G->m);
        for(i = 0; i < G->m; i++){
            if(i != G->m-1)
                fprintf(fp,"%f,",z[i]);
            else
                fprintf(fp,"%f};\n",z[i]);
        }
        
        fprintf(fp,"/* RHS */\n");
        fprintf(fp,"pfloat bx[%i] = {",A->n);
        for(i = 0; i < A->n; i++){
            if(i != A->n-1)
                fprintf(fp,"%f,",bx[i]);
            else
                fprintf(fp,"%f};\n",bx[i]);
        }
        fprintf(fp,"pfloat by[%i] = {",A->m);
        for(i = 0; i < A->m; i++){
            if(i != A->m-1)
                fprintf(fp,"%f,",by[i]);
            else
                fprintf(fp,"%f};\n",by[i]);
        }
        fprintf(fp,"pfloat bz[%i] = {",G->m);
        for(i = 0; i < G->m; i++){
            if(i != G->m-1)
                fprintf(fp,"%f,",bz[i]);
            else
                fprintf(fp,"%f};\n",bz[i]);
        }
        fclose(fp);   
        FREE(A);
        FREE(G);     
    }
}
