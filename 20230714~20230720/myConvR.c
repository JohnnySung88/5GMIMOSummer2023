#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variable declarations here */
    double *x, *h, *y;  /* pointers to input/output arrays */
    int m, n, i, j;

    /* Check for proper number of input and output arguments */    
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:mexfunc:invalidNumInputs",
                "Two input arguments required.");
    } 
    if (nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:mexfunc:maxlhs",
                "Too many output arguments.");
    }

    /* Check if the inputs are of the correct type */
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt( "MATLAB:mexfunc:inputNotDouble",
                "All input arguments must be of type double.");
    }

    /* Get the length of the sequences */
    m = mxGetNumberOfElements(prhs[0]);
    n = mxGetNumberOfElements(prhs[1]);

    /* Create output array y */
    plhs[0] = mxCreateDoubleMatrix(1, m+n-1, mxREAL);

    /* Get pointers to input arrays x and h, and output array y */
    x = mxGetPr(prhs[0]);
    h = mxGetPr(prhs[1]);
    y = mxGetPr(plhs[0]);

    /* Perform the convolution */
    for(i = 0; i < m+n-1; i++)
    {
        y[i] = 0;
        for(j = 0; j <= i; j++)
        {
            if(j < m && (i-j) < n) {
                y[i] += x[j]*h[i-j];
            }
        }
    }
}