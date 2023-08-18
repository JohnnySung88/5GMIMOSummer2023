#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
    double *real;
    double *imag;
} Complex;

void LDPC(int xa, int xb, int ha, int hb, Complex X, Complex H, Complex output) {
    Complex Eta, S, Lumda;
    Eta.real = malloc(ha*hb*sizeof(double));
    S.real = malloc(hb*sizeof(double)); 
    Lumda.real = malloc(xb*sizeof(double));

    int i, j, k, f, v;
    int iter = 5;
    int idx_pl, Ssgn;
    double Spsi, Qtemp, Qtemppsi, Qtempsgn;
    int *idx = calloc(hb, sizeof(int));
    int *Ssgn_temp = calloc(hb, sizeof(int));

    double Little_num = pow(0.1, 25);

    for(i = 0; i < xb; i++) {
        Lumda.real[i] = X.real[i];
    }

    memset(Eta.real, 0, ha*hb*sizeof(double));

    for(i = 0; i < iter; i++) {
        for(j = 0; j < ha; j++) {
            idx_pl = 0;
            for(k = 0; k < hb; k++) {
                if(H.real[j + k*324] == 1) {
                    idx[idx_pl] = k;
                    idx_pl++;
                }
            }
            
            Ssgn = 1;
            Spsi = 0;
            for(f = 0; f < idx_pl; f++) {
                S.real[f] = Lumda.real[idx[f]] - Eta.real[j + idx[f]*324];
                Ssgn_temp[f] = (S.real[f] > 0) ? 1 : -1;
                Ssgn *= Ssgn_temp[f];
                Spsi += -log(Little_num + tanh(fabs(S.real[f] * 0.5)));
            }
            
            for(v = 0; v < idx_pl; v++) {
                Qtemp = Lumda.real[idx[v]] - Eta.real[j + idx[v]*324];
                Qtemppsi = -log(Little_num + tanh(fabs(S.real[v]*0.5)));
                Qtempsgn = (Qtemp > 0) ? Ssgn : -Ssgn;
                Eta.real[j + idx[v]*324] = Qtempsgn * (-log(Little_num + tanh(fabs(Spsi - Qtemppsi) * 0.5)));
                Lumda.real[idx[v]] = Qtemp + Eta.real[j + idx[v]*324];
            }
        }
    }

    for(i = 0; i < ha; i++) {
        output.real[i] = (Lumda.real[i] > 0) ? 0 : 1;
    }

    free(Eta.real);
    free(S.real);
    free(Lumda.real);
    free(idx);
    free(Ssgn_temp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Complex X, H, Output;
    int Xa = mxGetM(prhs[0]);
    int Xb = mxGetN(prhs[0]);
    int Ha = mxGetM(prhs[1]);
    int Hb = mxGetN(prhs[1]);

    X.real = mxGetPr(prhs[0]);
    H.real = mxGetPr(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(1, Ha, mxCOMPLEX);
    Output.real = mxGetPr(plhs[0]);
    
    LDPC(Xa, Xb, Ha, Hb, X, H, Output);
}
