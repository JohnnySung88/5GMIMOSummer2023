#include "mex.h"
#include "complex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variable declarations */
    double *xr, *xi, *hr, *hi; /* pointers to input arrays */
    double *yr, *yi;           /* pointers to output arrays */
    double complex *x, *h, *y; /* pointers to complex input/output arrays */
    int m, n, i, j, L;

    /* Check for proper number of input and output arguments */
    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("MATLAB:mexfunc:invalidNumInputs",
                          "Two input arguments required.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("MATLAB:mexfunc:maxlhs",
                          "Too many output arguments.");
    }

    /* Check if the inputs are of the correct type */
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
    {
        mexErrMsgIdAndTxt("MATLAB:mexfunc:inputNotDouble",
                          "All input arguments must be of type double.");
    }

    /* Get the length of the sequences */
    m = mxGetNumberOfElements(prhs[0]);
    n = mxGetNumberOfElements(prhs[1]);

    /* Check if inputs are complex */
    bool isComplex = mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]);

    L = m + n - 1; // length of output

    if (isComplex)
    {
        /* Create complex output array y */
        plhs[0] = mxCreateDoubleMatrix(1, L, mxCOMPLEX);

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
        y = malloc(L * sizeof(double complex));

        /* Construct complex input arrays */
        for (i = 0; i < m; i++)
            x[i] = xr[i] + I * (xi ? xi[i] : 0);

        for (i = 0; i < n; i++)
            h[i] = hr[i] + I * (hi ? hi[i] : 0);

        /* Perform the cross-correlation */
        for (i = -(m - 1); i < n; i++)
        {
            y[L - (i + m - 1) - 1] = 0;
            for (j = 0; j < m; j++)
            {
                if (i + j >= 0 && i + j < n)
                {
                    y[L - (i + m - 1) - 1] += x[j] * conj(h[i + j]); // corrected here, use conjugate of h
                }
            }
            /* Copy the real and imaginary parts to output array */
            yr[L - (i + m - 1) - 1] = creal(y[L - (i + m - 1) - 1]);
            yi[L - (i + m - 1) - 1] = cimag(y[L - (i + m - 1) - 1]);
        }
    }
    else
    {
        /* Create real output array y */
        plhs[0] = mxCreateDoubleMatrix(1, L, mxREAL);

        /* Get pointers to input arrays x and h, and output array y */
        xr = mxGetPr(prhs[0]);
        hr = mxGetPr(prhs[1]);
        yr = mxGetPr(plhs[0]);

        /* Allocate memory for real input/output arrays */
        x = malloc(m * sizeof(double complex));
        h = malloc(n * sizeof(double complex));
        y = malloc(L * sizeof(double complex));

        /* Construct real input arrays */
        for (i = 0; i < m; i++)
            x[i] = xr[i];

        for (i = 0; i < n; i++)
            h[i] = hr[i];

        /* Perform the cross-correlation */
        for (i = -(m - 1); i < n; i++)
        {
            y[L - (i + m - 1) - 1] = 0;
            for (j = 0; j < m; j++)
            {
                if (i + j >= 0 && i + j < n)
                {
                    y[L - (i + m - 1) - 1] += x[j] * h[i + j];
                }
            }
            /* Copy the real and imaginary parts to output array */
            yr[L - (i + m - 1) - 1] = creal(y[L - (i + m - 1) - 1]);
        }
    }

    /* Free allocated memory */
    free(x);
    free(h);
    free(y);
}
