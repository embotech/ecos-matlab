#include "mex.h"
#include "matrix.h"
#include "lino_kkt.h"
#include <stdio.h>

/* The mex-function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if(nrhs == 11){
        idxint* socdims = NULL;
        /* Get inputs */
        mxArray* Amx = (mxArray*)prhs[0];
        mxArray* Gmx = (mxArray*)prhs[1];
        pfloat* s = (pfloat*)mxGetPr(prhs[2]);
        pfloat* z = (pfloat*)mxGetPr(prhs[3]);
        idxint nlp = (idxint)mxGetScalar(prhs[4]);
        idxint nsoc = (idxint)mxGetScalar(prhs[5]);
        if(nsoc){
            socdims = mxGetData(prhs[6]);
        }
        pfloat delta = (pfloat)mxGetScalar(prhs[7]);
        pfloat* bx = mxGetData(prhs[8]);
        pfloat* by = mxGetData(prhs[9]);
        pfloat* bz = mxGetData(prhs[10]);
        spmat* A = createSparseMatrix((mwIndex)mxGetM(Amx),(mwIndex)mxGetN(Amx),(mwIndex)mxGetJc(Amx)[mxGetN(Amx)],(mwIndex*)mxGetJc(Amx),(mwIndex*)mxGetIr(Amx),(pfloat*)mxGetPr(Amx));
        spmat* G = createSparseMatrix((mwIndex)mxGetM(Gmx),(mwIndex)mxGetN(Gmx),(mwIndex)mxGetJc(Gmx)[mxGetN(Gmx)],(mwIndex*)mxGetJc(Gmx),(mwIndex*)mxGetIr(Gmx),(pfloat*)mxGetPr(Gmx));
        idxint i, j, it;

        mxArray *dxmx, *dymx, *dzmx;
        mwIndex *dxjc, *dxir, *dyjc, *dyir, *dzjc, *dzir;
        pfloat *dxpr, *dypr, *dzpr;
        
         /* Allocate memory for Cone */
        cone* C = mxMalloc(sizeof(cone));
        C->lpc = mxMalloc(sizeof(lpcone));
        if(nsoc){
            C->soc = mxMalloc(nsoc*sizeof(socone));
        }
        /* mexPrintf("sizeof(C) = %i\nsizeof(lpc) = %i\nsizeof(soc) = %i\n",sizeof(*C),sizeof(C->lpc),nsoc*sizeof(C->soc)); */

        for(i = 0; i < nsoc; i++){
            C->soc[i].skbar = mxMalloc(socdims[i]*sizeof(pfloat));
            C->soc[i].zkbar = mxMalloc(socdims[i]*sizeof(pfloat));
            C->soc[i].q = mxMalloc((socdims[i]-1)*sizeof(pfloat));
        }    
        
        C->lpc->w = mxMalloc(nlp*sizeof(pfloat));
        C->lpc->v = mxMalloc(nlp*sizeof(pfloat));
        if(nsoc){
            C->nsoc = nsoc;
        }
        else{
            C->nsoc = 0;
        }
        C->lpc->p = nlp;

        for(i = 0; i < nsoc; i++){
            C->soc[i].p = socdims[i];
        }

        pfloat lambda[G->m];        
        if(updateScalings(C,s,z,lambda) == OUTSIDE_CONE){
            mexErrMsgTxt("Outside cone!");
        }
        /* Setup and factor */
        
        pfc* mypfc;
        if(nsoc){
            mypfc = neSetup(nlp,nsoc,socdims,G,A,delta);
        }
        else{
            mypfc = neSetup(nlp,0,0,G,A,delta);
        };
        NEfactor(mypfc,C);
        it = NEsolve(mypfc,C,bx,by,bz);
        /* Free memory of cone */
        for(i = 0; i < C->nsoc; i++){
           mxFree(C->soc[i].skbar);
           mxFree(C->soc[i].zkbar);
           mxFree(C->soc[i].q);
        }
        if(nsoc){
            mxFree(C->soc);           
        }
        mxFree(C->lpc->w);
        mxFree(C->lpc->v);
        mxFree(C->lpc);
        mxFree(C);
        
        /* copy solution */        
        dxmx = mxCreateSparse((mwSize)(A->n),1,(mwSize)(A->n),mxREAL);
        dxpr = mxGetPr(dxmx); 
        dxjc = mxGetJc(dxmx); 
        dxir = mxGetIr(dxmx);
        
        dxjc[0] = 0;
        dxjc[1] = (mwIndex)(A->n);

        for(i = 0; i < dxjc[1]; i++){
            dxpr[i] = (pfloat)mypfc->dx[i];
            dxir[i] = i;
        }
        
        dymx = mxCreateSparse((mwSize)(A->m),1,(mwSize)(A->m),mxREAL);
        dypr = mxGetPr(dymx); 
        dyjc = mxGetJc(dymx); 
        dyir = mxGetIr(dymx);
        
        dyjc[0] = 0;
        dyjc[1] = (mwIndex)(A->m);

        for(i = 0; i < dyjc[1]; i++){
            dypr[i] = (pfloat)mypfc->dy[i];
            dyir[i] = i;
        }
        
        dzmx = mxCreateSparse((mwSize)(G->m),1,(mwSize)(G->m),mxREAL);
        dzpr = mxGetPr(dzmx); 
        dzjc = mxGetJc(dzmx); 
        dzir = mxGetIr(dzmx);
        
        dzjc[0] = 0;
        dzjc[1] = (mwIndex)(G->m);

        for(i = 0; i < dzjc[1]; i++){
            dzpr[i] = (pfloat)mypfc->dz[i];
            dzir[i] = i;
        }
        
        
        if(nlhs == 3){
            plhs[0] = dxmx;
            plhs[1] = dymx;
            plhs[2] = dzmx;                 
        }
        
        if(nlhs > 3){
            plhs[0] = dxmx;
            plhs[1] = dymx;
            plhs[2] = dzmx;
            plhs[3] = mxCreateDoubleScalar(it);        
        }
        neCleanup(mypfc,nsoc,nlp);
    }
}