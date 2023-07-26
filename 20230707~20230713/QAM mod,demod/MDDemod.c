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

// 將二進位轉換為格雷碼
int binaryToGray32(int num) {
    num = num ^ (num >> 1);
    return num;
}

// 將數字轉換為位置
int numToPos(int num, int size) {
    int ans;
    if (num > 0) { // 在右側
        ans = size / 2 + (num - 1) / 2;
    } else { // 在左側
        ans = size / 2 - 1 - (-num - 1) / 2;
    }
    return ans;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int col_NUM, row_NUM; 
    double *NUM_real, *NUM_imag;
    double *M;
    col_NUM = mxGetN(NUM_IN);
    row_NUM = mxGetM(NUM_IN);
    X_OUT = mxCreateDoubleMatrix(row_NUM, col_NUM, mxREAL);
    NUM_real = mxGetPr(NUM_IN);
    NUM_imag = mxGetPi(NUM_IN);
    M = mxGetPr(M_IN);
    double* Ans = mxGetPr(X_OUT);

    int k;
    double complex ram[col_NUM];

    int QAM = (int)M[0]; // QAM 點個數
    int BIT = (int)(log(QAM) / log(2)); // QAM Bit數
    int h_size = (int)pow(2, BIT / 2); // QAM 星座圖的一半大小
    int r_mask = (1 << (BIT / 2)) - 1; // 右半部的位遮罩
    int l_mask = r_mask << (BIT / 2); // 左半部的位遮罩
    int r_pos, l_pos;

    double complex points[QAM]; // QAM 星座點陣列

    // 生成 QAM 星座點
    for (int i = 0; i < h_size; i++) {
        for (int j = 0; j < h_size; j++) {
            points[i*h_size + j] = (2.0 * i - h_size + 1) + I * (2.0 * j - h_size + 1);
        }
    }

    for (k = 0; k < col_NUM; k++) {
        double complex z = NUM_real[k] + I * NUM_imag[k];
        double min_dist = INFINITY;
        int min_index = 0;

        // 迭代所有星座點，找到與輸入信號距離最小的星座點
        for (int i = 0; i < QAM; ++i) {
            double dist = cabs(z - points[i]); // 計算歐幾里得距離
            if (dist < min_dist) {
                min_dist = dist;
                min_index = i;
            }
        }

        // 實部的格雷碼位置
        l_pos = numToPos((int)round(creal(points[min_index])), h_size);
        l_pos = binaryToGray32(l_pos) << (BIT / 2);
        // 虛部的格雷碼位置
        r_pos = numToPos((int)round(-cimag(points[min_index])), h_size);
        r_pos = binaryToGray32(r_pos);
        // 計算最終解碼結果
        Ans[k] = l_pos + r_pos;
    }
}