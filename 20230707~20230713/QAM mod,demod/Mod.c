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

unsigned int grayToBinary32(unsigned int num) {
    num = num ^ (num >> 16);
    num = num ^ (num >> 8);
    num = num ^ (num >> 4);
    num = num ^ (num >> 2);
    num = num ^ (num >> 1);
    return num;
}

int posToNum(int pos, int size) {
    int ans;
    if (pos >= size / 2) { // Right side
        if (pos == size / 2)
            ans = 1;
        else
            ans = 1 + 2 * (pos - size / 2);
    } else { // Left side
        if (pos == size / 2 - 1)
            ans = -1;
        else
            ans = -1 - 2 * (size / 2 - pos - 1);
    }
    return ans;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int col_NUM, row_NUM; // Added row_NUM declaration
    double *NUM;
    double *M;
    double *N;
    col_NUM = mxGetN(NUM_IN);
    row_NUM = mxGetM(NUM_IN);
    X_OUT = mxCreateDoubleMatrix(row_NUM, col_NUM, mxCOMPLEX);
    Ans.real = mxGetPr(X_OUT);
    Ans.imag = mxGetPi(X_OUT);
    NUM = mxGetPr(NUM_IN);
    M = mxGetPr(M_IN);

    int k;
    double complex ram[col_NUM];

    int QAM = (int)M[0];
    int BIT = (int)(log(QAM) / log(2));
    int mask = (int)QAM - 1;
    int h_size = (int)pow(2, BIT / 2);
    int r_mask = mask >> BIT / 2;
    int l_mask = r_mask << BIT / 2;
    int r_gray, l_gray;
    int r_pos, l_pos;

    for (k = 0; k < col_NUM; k++) {
        // Init
        ram[k] = 0;
        // Real part
        l_gray = ((int)NUM[k] & l_mask) >> BIT / 2;
        l_pos = grayToBinary32(l_gray);
        ram[k] += posToNum(l_pos, h_size);
        // Imaginary part
        r_gray = (int)NUM[k] & r_mask;
        r_pos = grayToBinary32(r_gray);
        ram[k] += posToNum(r_pos, h_size) * -I;
        // Output
        Ans.real[k] = creal(ram[k]);
        Ans.imag[k] = cimag(ram[k]);
    }
}