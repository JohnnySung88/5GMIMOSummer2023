#include "mex.h"
#include <lapack.h>
#include <complex.h>
//NOTICE: Using command:"mex CD.c -llapack" to compile

typedef double complex Complex;

void inverseMatrixReal(int n, double *matrix, double *invMatrix) {
    char uplo = 'U';  
    int info;

    dpotrf_(&uplo, &n, matrix, &n, &info);
    if (info > 0) {
        mexErrMsgTxt("Cholesky decomposition failed. The matrix is not positive definite.");
    }

    dpotri_(&uplo, &n, matrix, &n, &info);
    if (info > 0) {
        mexErrMsgTxt("Matrix inversion failed. The matrix is singular.");
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            invMatrix[i*n + j] = matrix[i*n + j];
            invMatrix[j*n + i] = matrix[i*n + j];
        }
    }
}

void inverseMatrixComplex(int n, Complex *matrix, Complex *invMatrix) {
    char uplo = 'U';  
    int info;

    zpotrf_(&uplo, &n, matrix, &n, &info);
    if (info > 0) {
        mexErrMsgTxt("Cholesky decomposition failed. The matrix is not positive definite.");
    }

    zpotri_(&uplo, &n, matrix, &n, &info);
    if (info > 0) {
        mexErrMsgTxt("Matrix inversion failed. The matrix is singular.");
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            invMatrix[i*n + j] = matrix[i*n + j];
            invMatrix[j*n + i] = conj(matrix[i*n + j]);
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1) {
        mexErrMsgTxt("Invalid number of input arguments. Expected 1 input.");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("Invalid number of output arguments. Expected 1 output.");
    }

    mwSize n = mxGetN(prhs[0]);

    if (mxIsComplex(prhs[0])) {
        double *realPart = mxGetPr(prhs[0]);
        double *imagPart = mxGetPi(prhs[0]);

        Complex *matrix = mxMalloc(n * n * sizeof(Complex));
        for (int i = 0; i < n * n; i++) {
            matrix[i] = realPart[i] + I*imagPart[i];
        }

        mxArray *outMatrix = mxCreateDoubleMatrix(n, n, mxCOMPLEX);
        double *realInvPart = mxGetPr(outMatrix);
        double *imagInvPart = mxGetPi(outMatrix);

        Complex *invMatrix = mxMalloc(n * n * sizeof(Complex));
        inverseMatrixComplex(n, matrix, invMatrix);

        for (int i = 0; i < n * n; i++) {
            realInvPart[i] = creal(invMatrix[i]);
            imagInvPart[i] = cimag(invMatrix[i]);
        }

        mxFree(matrix);
        mxFree(invMatrix);
        plhs[0] = outMatrix;
    }
    else {
        double *matrix = mxGetPr(prhs[0]);
        mxArray *outMatrix = mxCreateDoubleMatrix(n, n, mxREAL);
        double *invMatrix = mxGetPr(outMatrix);

        inverseMatrixReal(n, matrix, invMatrix);
        plhs[0] = outMatrix;
    }
}
