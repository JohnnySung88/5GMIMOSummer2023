#include "mex.h"
#include <math.h>

#define L(x, y) L[(x) + (y)*n]
#define M(x, y) M[(x) + (y)*n]

/* The computational routine */
void myCholeskyC(double *M, double *L, mwSize n)
{

    mwSize i;
    mwSize j;
    mwSize k;
    mwSize l;

    /* multiply each element y by x */
    for (i = 0; i < n; i++)
    {

        double sum = 0;
        for (k = 0; k < n; k++)
        {
            sum += L(i, k) * L(i, k);
        } // end of for-loop k
        L(i, i) = sqrt(M(i, i) - sum);

        for (j = i + 1; j < n; j++)
        {
            double sum2 = 0;
            for (l = 0; l < n; l++)
            {
                sum2 += L(i, l) * L(j, l);
            } // end of for-loop l

            L(j, i) = (M(j, i) - sum2) / L(i, i);
        } // end of for-loop j
    }     // end of for-loop i
} // end of computational routine

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double multiplier; /* input scalar */
    double *inMatrix;  /* 1xN input matrix */
    size_t ncols;      /* size of matrix */
    double *outMatrix; /* output matrix */

    /* check for proper number of arguments */
    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Only one input required.");
    }
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required.");
    }
    /* make sure the input argument is type double */
    if (!mxIsDouble(prhs[0]) ||
        mxIsComplex(prhs[0]))
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "Input matrix must be type double.");
    }

    /* create a pointer to the real data in the input matrix  */
    inMatrix = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)ncols, (mwSize)ncols, mxREAL);

    /* set output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    myCholeskyC(inMatrix, outMatrix, (mwSize)ncols);
}