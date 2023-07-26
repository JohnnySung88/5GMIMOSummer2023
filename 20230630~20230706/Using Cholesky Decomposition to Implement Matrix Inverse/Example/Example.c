#include "mex.h"

typedef struct
{
    double *real;
    double *imag;
} Complex;
Complex A;
Complex Ainv;

void Inverse(int n, Complex A, Complex Ainv)
{
    /*
    L 和 Linv 變數與 A 和 Ainv 變數是不同的變數，它們分別用於不同的目的。這兩組變數在記憶體中佔用不同的空間，並且在函式執行期間可以獨立操作。
    在 Inverse 函式中，你可以根據你的需求，使用 L 和 Linv 變數進行計算反矩陣的相關操作。在函式結束時，你可以使用 free 函式來釋放 L 和 Linv 變數所使用的記憶體。
    總結起來，Complex 型別的變數 L 和 Linv 是區域變數，用於在 Inverse 函式中存儲計算反矩陣的中間結果。
    */
    Complex L, Linv;
    L.real = malloc(n * n * sizeof(double));
    L.imag = malloc(n * n * sizeof(double));
    // put your code here
    free(L.real);
    free(L.imag);
}

/*
這是 MEX 檔案的入口函式mexFunction，它接受四個參數：
nlhs 是輸出參數的數量， plhs 是指向輸出參數的指標陣列，
nrhs 是輸入參數的數量， prhs 是指向輸入參數的指標陣列。
nlhs 是 "Number of Left Hand Side" 的縮寫，它指定了在 MATLAB 中可以接收多少個輸出參數。在這個函式中，nlhs 用於確定要創建多少個輸出參數的指標陣列。
plhs 是 "Pointer to Left Hand Side" 的縮寫，它指定了要返回給 MATLAB 的輸出參數。在這個函式中，plhs 用於存儲計算結果，以供 MATLAB 使用。
nrhs 是 "Number of Right Hand Side" 的縮寫，它指定了在 MATLAB 中傳遞給這個 MEX 函式的輸入參數的數量。
prhs 是 "Pointer to Right Hand Side" 的縮寫，它指定了在 MATLAB 中傳遞給這個 MEX 函式的輸入參數。這些參數可以是數組、矩陣、結構等 MATLAB 中的任意類型。
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // 在函式內部，首先聲明了兩個整數變數 M 和 N，它們將用於存儲輸入複數矩陣 A 的維度。
    int M, N;

    // 接著，使用 mxGetPr 函式和 mxGetPi 函式分別獲取輸入複數矩陣 A 的實部和虛部，並將它們分別存儲到 A.real 和 A.imag 中。
    A.real = mxGetPr(prhs[0]);
    A.imag = mxGetPi(prhs[0]);

    // 然後，使用 mxCreateDoubleMatrix 函式創建一個 M 行 N 列的複數矩陣，並將其存儲到 plhs[0] 中。
    M = mxGetM(prhs[0]); // row
    N = mxGetN(prhs[0]); // col
    plhs[0] = mxCreateDoubleMatrix(M, N, mxCOMPLEX);

    // 最後，使用 mxGetPr 函式和 mxGetPi 函式分別獲取 Ainv 矩陣的實部和虛部，並將它們分別存儲到 Ainv.real 和 Ainv.imag 中，這樣就完成了將輸入複數矩陣轉換為輸出複數矩陣的操作。
    Ainv.real = mxGetPr(plhs[0]);
    Ainv.imag = mxGetPi(plhs[0]);
}
