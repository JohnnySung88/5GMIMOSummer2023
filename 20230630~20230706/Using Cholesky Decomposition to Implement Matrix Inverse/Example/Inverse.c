#include "mex.h"
#include <lapack.h>
//NOTICE: Using command:"mex Inverse.c -llapack" to compile

void inverseMatrix(int n, double *matrix, double *invMatrix)
{
    int *IPIV, LWORK, INFO;
    double *WORK;

    IPIV = (int*)mxMalloc(n * sizeof(int));
    LWORK = n * n;
    WORK = (double*)mxMalloc(LWORK * sizeof(double));

    for (int i = 0; i < n * n; i++) {
        invMatrix[i] = matrix[i];
    }

    dgetrf_(&n, &n, invMatrix, &n, IPIV, &INFO);
    dgetri_(&n, invMatrix, &n, IPIV, WORK, &LWORK, &INFO);

    mxFree(IPIV);
    mxFree(WORK);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("Invalid number of input arguments. Expected 1 input.");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("Invalid number of output arguments. Expected 1 output.");
    }

    double *matrix = mxGetPr(prhs[0]);
    mwSize n = mxGetN(prhs[0]);

    mxArray *outMatrix = mxCreateDoubleMatrix(n, n, mxREAL);
    double *invMatrix = mxGetPr(outMatrix);

    inverseMatrix(n, matrix, invMatrix);

    plhs[0] = outMatrix;
}