#include "mex.h"
#include <lapack.h>
//NOTICE: Using command:"mex CDR.c -llapack" to compile

void inverseMatrix(int n, double *matrix, double *invMatrix)
{
    char uplo = 'U';  // Use upper triangular part of the matrix
    int info;

    // Perform Cholesky decomposition
    dpotrf_(&uplo, &n, matrix, &n, &info);

    if (info > 0) {
        mexErrMsgTxt("Cholesky decomposition failed. The matrix is not positive definite.");
    }

    // Compute inverse using Cholesky decomposition
    dpotri_(&uplo, &n, matrix, &n, &info);

    if (info > 0) {
        mexErrMsgTxt("Matrix inversion failed. The matrix is singular.");
    }

    // Compute the inverse matrix by copying the upper triangular part
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            invMatrix[i*n + j] = matrix[i*n + j];
            invMatrix[j*n + i] = matrix[i*n + j];
        }
    }
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
