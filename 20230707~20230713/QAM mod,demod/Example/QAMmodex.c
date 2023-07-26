#include "mex.h"
#include <math.h>

// Function to compute M-ary QAM modulation with Gray encoding
void myqammod(double *inputSignal, int inputSignalLength, int M, double *outputSignal)
{
    int bitsPerSymbol = log2(M); // Number of bits per symbol
    int grayCodedSymbol;
    double normalization = sqrt(M);
    int i, j, bit;

    for (i = 0; i < inputSignalLength; i++)
    {
        grayCodedSymbol = 0;

        // Gray encode the input
        for (j = 0; j < bitsPerSymbol; j++)
        {
            bit = ((int)inputSignal[i] >> j) & 1;                             // Get jth bit of the ith input
            grayCodedSymbol |= (bit ^ ((int)inputSignal[i] >> (j + 1))) << j; // Compute the Gray code
        }

        // Now compute the QAM symbol
        outputSignal[2 * i] = ((grayCodedSymbol % ((int)sqrt(M))) * 2 - M + 1) / normalization;     // I
        outputSignal[2 * i + 1] = ((grayCodedSymbol / ((int)sqrt(M))) * 2 - M + 1) / normalization; // Q
    }
}

// This function will be called from MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *inputSignal;  // input array
    int M;                // M-ary QAM value
    double *outputSignal; // output array
    int i, inputSignalLength;

    // Check for proper number of arguments
    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Two inputs required.");
    }

    // Check that both inputs are noncomplex double scalar values
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) ||
        !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) ||
        mxGetNumberOfElements(prhs[1]) != 1 ||
        (M = (int)mxGetScalar(prhs[1])) != 4 && M != 16 && M != 64)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                          "Input matrix must be noncomplex double and M must be 4, 16 or 64.");
    }

    // Get the length of the input array
    inputSignalLength = mxGetN(prhs[0]);

    // Create matrix for the return argument
    plhs[0] = mxCreateDoubleMatrix(1, 2 * inputSignalLength, mxREAL); // Double the length for I and Q

    // Assign pointers to each input and output
    inputSignal = mxGetPr(prhs[0]);
    outputSignal = mxGetPr(plhs[0]);

    // Call the QAM modulation function
    myqammod(inputSignal, inputSignalLength, M, outputSignal);
}
