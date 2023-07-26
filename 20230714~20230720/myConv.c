#include "mex.h"
#include "complex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variable declarations */
    double *xr, *xi, *hr, *hi;  /* pointers to input arrays */
    double *yr, *yi; /* pointers to output arrays */
    double complex *x, *h, *y; /* pointers to complex input/output arrays */
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

    /* Check if inputs are complex */
    bool isComplex = mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]);

    if (isComplex) {
        /* Create complex output array y */
        plhs[0] = mxCreateDoubleMatrix(1, m+n-1, mxCOMPLEX);

        /* Get pointers to input arrays x and h, and output array y */
        xr = mxGetPr(prhs[0]);
        xi = mxGetPi(prhs[0]);
        hr = mxGetPr(prhs[1]);
        hi = mxGetPi(prhs[1]);
        yr = mxGetPr(plhs[0]);
        yi = mxGetPi(plhs[0]);

        /* Allocate memory for complex input/output arrays */
        x = malloc(m * sizeof(double complex));
        h = malloc(n * sizeof(double complex));
        y = malloc((m+n-1) * sizeof(double complex));

        /* Construct complex input arrays */
        for(i = 0; i < m; i++)
            x[i] = xr[i] + I* (xi ? xi[i] : 0);

        for(i = 0; i < n; i++)
            h[i] = hr[i] + I* (hi ? hi[i] : 0);
    }
    else {
        /* Create real output array y */
        plhs[0] = mxCreateDoubleMatrix(1, m+n-1, mxREAL);

        /* Get pointers to input arrays x and h, and output array y */
        xr = mxGetPr(prhs[0]);
        hr = mxGetPr(prhs[1]);
        yr = mxGetPr(plhs[0]);

        /* Allocate memory for real input/output arrays */
        x = malloc(m * sizeof(double complex));
        h = malloc(n * sizeof(double complex));
        y = malloc((m+n-1) * sizeof(double complex));

        /* Construct real input arrays */
        for(i = 0; i < m; i++)
            x[i] = xr[i];

        for(i = 0; i < n; i++)
            h[i] = hr[i];
    }

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
        /* Copy the real and imaginary parts to output array */
        yr[i] = creal(y[i]);
        if (isComplex) {
            yi[i] = cimag(y[i]);
        }
    }

    /* Free allocated memory */
    free(x);
    free(h);
    free(y);
}
