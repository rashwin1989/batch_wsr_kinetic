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
          1) P: pressure
          2) T: temperature
          3) Pc vector: critical pressure
          4) Tc vector: critical temperature
          5) w vector: acentric factor
          6) tk vector: type of k, binary interaction parameters
                        check subroutine calculate_kij in PR_EoS.f90 for details
          7) MW vector
          8) x vector: molar fractions
   output:
          1) rho
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *P, *T, *Pc, *Tc, *w, *x, *MW, *coef_ab, *rho;
    mwSize  rows, cols, n;
    
    int *ntk;
    double *tk;

    int i;

    int bKijSet = 1;
    double  *kij_tmp;
    
    /* check for the proper number of arguments */
    if(nrhs != 8)
      mexErrMsgTxt("8 inputs required.");
    if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");
    /*Check that both inputs are row vectors*/
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1  
      || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1  )
      mexErrMsgTxt("1st - 2nd inputs must be scalars.");

    rows=mxGetM(prhs[2]);
    cols=mxGetN(prhs[2]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("3rd input must be row or column vectors.");

    if(  mxGetM(prhs[3]) != rows || mxGetN(prhs[3]) != cols  
      || mxGetM(prhs[4]) != rows || mxGetN(prhs[4]) != cols  
      || mxGetM(prhs[5]) != rows || mxGetN(prhs[5]) != cols  
      || mxGetM(prhs[6]) != rows || mxGetN(prhs[6]) != cols 
      || mxGetM(prhs[7]) != rows || mxGetN(prhs[7]) != cols )
      mexErrMsgTxt("4th to 8th inputs must have the same shape.");

    /* get pointers to the inputs */
    P  = mxGetPr(prhs[0]); 
    T  = mxGetPr(prhs[1]);
    Pc = mxGetPr(prhs[2]);
    Tc = mxGetPr(prhs[3]);
    w  = mxGetPr(prhs[4]);
    tk = mxGetPr(prhs[5]);
    MW = mxGetPr(prhs[6]);
    x  = mxGetPr(prhs[7]);

    coef_ab = (double*)malloc(sizeof(double)*n);

    ntk= (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) {
      ntk[i]=tk[i];
      coef_ab[i]=-1;
    }

    kij_tmp = (double*)malloc(sizeof(double)*n*n);

    /* create a new array and set the output pointer to it */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    rho = mxGetPr(plhs[0]);

    /*
    printf("P=%lf T=%lf\n", *P, *T);
    printf("Pc:");for (i=0; i<n; i++) printf("%f ", Pc[i]);printf("\n");
    printf("Tc:");for (i=0; i<n; i++) printf("%f ", Tc[i]);printf("\n");
    printf("w :");for (i=0; i<n; i++) printf("%f ", w[i]);printf("\n");
    printf("tk:");for (i=0; i<n; i++) printf("%d ",ntk[i]);printf("\n");
    printf("x :");for (i=0; i<n; i++) printf("%f ", x[i]);printf("\n");
    */

    /* call the C subroutine */
    calculate_kij_(&bKijSet,T,&n,Pc,Tc,w,ntk,kij_tmp);
    density_pr_eos_(P,T,&n,Pc,Tc,w,MW,x,ntk,coef_ab,rho);

    free(coef_ab);
    free(ntk);
    free(kij_tmp);
    return;
}
