#include "mex.h"
#include <complex.h>
#include <math.h>

#define	NUM_IN	prhs[0]
#define	M_IN	prhs[1]
#define X_OUT	plhs[0]

typedef struct {
    double *real;
    double *imag;
} Complex2;
Complex2 Ans;

int binaryToGray32(int num) {
    num = num ^ (num >> 1);
    return num;
}

int numToPos(int num, int size) {
    int ans;
    if (num > 0) { // Right side
        ans = size / 2 + (num - 1) / 2;
    } else { // Left side
        ans = size / 2 - 1 - (-num - 1) / 2;
    }
    return ans;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int col_NUM, row_NUM; 
    double *NUM_real, *NUM_imag;
    double *M;
    double *N;
    col_NUM = mxGetN(NUM_IN);
    row_NUM = mxGetM(NUM_IN);
    X_OUT = mxCreateDoubleMatrix(row_NUM, col_NUM, mxREAL);
    NUM_real = mxGetPr(NUM_IN);
    NUM_imag = mxGetPi(NUM_IN);
    M = mxGetPr(M_IN);
    double* Ans = mxGetPr(X_OUT);

    int k;
    double complex ram[col_NUM];

    int QAM = (int)M[0];
    int BIT = (int)(log(QAM) / log(2));
    int h_size = (int)pow(2, BIT / 2);
    int r_mask = (1 << (BIT / 2)) - 1;
    int l_mask = r_mask << (BIT / 2);
    int r_pos, l_pos;

    for (k = 0; k < col_NUM; k++) {
        // Real part
        l_pos = numToPos((int)round(NUM_real[k]), h_size);
        l_pos = binaryToGray32(l_pos) << (BIT / 2);
        // Imaginary part
        r_pos = numToPos((int)round(-NUM_imag[k]), h_size);
        r_pos = binaryToGray32(r_pos);
        // Output
        Ans[k] = l_pos + r_pos;
    }
}
