#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

void make_L(int M, double complex *l , double complex *mat){
    double complex sum_2conj;
    for(int col = 0;col < M;col++){
        for(int row = 0;row < M;row++){
            if(row == col){
                sum_2conj = 0;
                for(int ram = 0; ram < col; ram++)
                    sum_2conj += l[row*M + ram] * conj(l[row*M + ram]);
                l[row*M + col] = csqrt(mat[row*M + col] - sum_2conj);
            }
            else if(col < row){
                sum_2conj = 0;
                for(int ram = 0; ram < col; ram++)
                    sum_2conj += l[row*M + ram] * conj(l[col*M + ram]);
                l[row*M + col] = (1/l[col*M + col]) * (mat[row*M + col] - sum_2conj);
            }
        }
    }
}

void invt_L(int M, double complex *l , double complex *L){
    double complex *R = (double complex*) malloc(M * M * sizeof(double complex));
    for(int row = 0;row < M;row++){
        for(int col = 0;col < M;col++){
            R[row*M + col] = L[row*M + col];
            if(row == col)
                l[row*M + col] = 1;
            else
                l[row*M + col] = 0;
        }
    }
    for(int row = 0;row < M;row++){
        for(int col = 0;col <= row;col++){
            if(row != col){
                for(int ram = 0;ram <= col;ram++)
                    l[row*M + ram] -= R[row*M + col] * l[col*M + ram];
                R[row*M + col] = 0;
            }
            else{
                for(int ram = 0;ram <= col;ram++)
                    l[row*M + ram] /= R[row*M + col];
                R[row*M + col] = 1;
            }
        }
    }
    free(R);
}

void M_conj(int M, double complex *a , double complex *l){
    double complex *H = (double complex*) malloc(M * M * sizeof(double complex));
    for(int row = 0;row < M;row++){
        for(int col = 0;col < M;col++){
            a[row*M + col] = 0;
            H[row*M + col] = conj(l[col*M + row]);
        }
    }
    for(int row = 0;row < M;row++)
        for(int col = 0;col < M;col++)
            for(int ram = 0;ram < M;ram++)
                a[row*M + col] += H[row*M + ram] * l[ram*M + col];
    free(H);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of input and output arguments */    
    if(nrhs!=1) {
        mexErrMsgIdAndTxt( "MATLAB:mexFunction:invalidNumInputs",
                "One input required.");
    } 
    if(nlhs!=1) {
        mexErrMsgIdAndTxt( "MATLAB:mexFunction:invalidNumOutputs",
                "One output required.");
    } 

    /* Input must be of complex type */
    if(!mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotComplex",
                          "Input must be complex.");
    }
    
    /* Get the size of the input matrix */
    int M = mxGetM(prhs[0]); //rows
    int N = mxGetN(prhs[0]); //columns

    /* We expect square matrices, check dimensions */
    if(M != N) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidInput",
                          "Input must be a square matrix.");
    }

    /* Create output matrix */
    plhs[0] = mxCreateDoubleMatrix(M, N, mxCOMPLEX);

    /* Create array of complex numbers for input and output matrices */
    double complex* inputMat = (double complex*) malloc(M * N * sizeof(double complex));
    double complex* outputMat = (double complex*) malloc(M * N * sizeof(double complex));
    double complex* L = (double complex*) malloc(M * N * sizeof(double complex));
    double complex* l = (double complex*) malloc(M * N * sizeof(double complex));

    double* inputReal = mxGetPr(prhs[0]);
    double* inputImag = mxGetPi(prhs[0]);
    double* outputReal = mxGetPr(plhs[0]);
    double* outputImag = mxGetPi(plhs[0]);

    /* Fill in the input matrix */
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            inputMat[j*M + i] = inputReal[j*M + i] + I * inputImag[j*M + i];
        }
    }

    make_L(M, L, inputMat);   
    invt_L(M, l, L);   
    M_conj(M, outputMat, l);   

    /* Write the result to the output matrix */
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            outputReal[j*M + i] = creal(outputMat[j*M + i]);
            outputImag[j*M + i] = cimag(outputMat[j*M + i]);
        }
    }

    /* Free allocated memory */
    free(inputMat);
    free(outputMat);
    free(L);
    free(l);
}
