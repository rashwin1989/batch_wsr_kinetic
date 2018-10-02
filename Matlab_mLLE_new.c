/* $Revision: 1.8.6.2 $ */
/*=========================================================
 * Example for illustrating how to use pass complex data 
 * from MATLAB to C and back again
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2006 The MathWorks, Inc.
 *=======================================================*/
#include "mex.h"

/*#include "PR_EoS.h" */

/* The gateway routine. */
/* to compute fugacity coefficient
   input: 
          1) c:  mixing ratios for initial values
          2) P: pressure
          3) T: temperature
          4) Pc vector: critical pressure
          5) Tc vector: critical temperature
          6) w vector: acentric factor
          7) tk vector: type of k, binary interaction parameters
                        check subroutine calculate_kij in PR_EoS.f90 for details
          8) xa0 vector: molar concentration, original 
          9) xb0 vector: molar concentration, original
         10) xai vector: molar concentration, initial guess
         11) xbi vector: molar concentration, initial guess
   output:
          0) xa vector: molar concentration
          1) xb vector: molar concentration
          2) c1 vector: molar concentration
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *P, *T, *Pc, *Tc, *w, *xa, *xb, *xa0, *xb0, *c, *c1;
    mwSize  rows, cols, n;
    
    int *ntk;
    double *tk;

    int i, n_miscible;

    int bKijSet = 1;
    double  *kij_tmp;
    
    /* check for the proper number of arguments */
    if(nrhs != 9)
      mexErrMsgTxt("9 inputs required.");
    if(nlhs > 3)
      mexErrMsgTxt("Too many output arguments.");
    /*Check that both inputs are row vectors*/
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1  
      || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1 
      || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1 )
      mexErrMsgTxt("1st - 3rd inputs must be scalars.");

    rows=mxGetM(prhs[3]);
    cols=mxGetN(prhs[3]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("4th input must be row or column vectors.");

    /*
    if(n==rows) {
      if (mxGetM(prhs[0]) != n-1)
        mexErrMsgTxt("1st input must be consitent with the 5th.");
    }else{
      if (mxGetN(prhs[0]) != n-1)
        mexErrMsgTxt("1st input must be consitent with the 5th.");
    }
    */

    if(  mxGetM(prhs[4]) != rows || mxGetN(prhs[4]) != cols  
      || mxGetM(prhs[5]) != rows || mxGetN(prhs[5]) != cols  
      || mxGetM(prhs[6]) != rows || mxGetN(prhs[6]) != cols  
      || mxGetM(prhs[7]) != rows || mxGetN(prhs[7]) != cols  
      || mxGetM(prhs[8]) != rows || mxGetN(prhs[8]) != cols)
      mexErrMsgTxt("5th to 9th inputs must have the same shape.");

    /* get pointers to the inputs */
    c  = mxGetPr(prhs[0]);
    P  = mxGetPr(prhs[1]); 
    T  = mxGetPr(prhs[2]);
    Pc = mxGetPr(prhs[3]);
    Tc = mxGetPr(prhs[4]);
    w  = mxGetPr(prhs[5]);
    tk = mxGetPr(prhs[6]);
    xa0= mxGetPr(prhs[7]);
    xb0= mxGetPr(prhs[8]);

    ntk= (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) {
      ntk[i]=tk[i];
    }

    kij_tmp = (double*)malloc(sizeof(double)*n*n);

    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    xa = mxGetPr(plhs[0]);
    xb = mxGetPr(plhs[1]);
    c1 = mxGetPr(plhs[2]);

    /* call the C subroutine */
    calculate_kij_(&bKijSet,T,&n,Pc,Tc,w,ntk,kij_tmp);
    species_lle4_(P,T,&n,Pc,Tc,w,ntk,xa0,xb0,c,c1,xa,xb,&n_miscible);

    free(kij_tmp);
    free(ntk);
    return;
}
