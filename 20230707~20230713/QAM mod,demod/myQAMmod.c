#include "mex.h"
#include <math.h>
#include <stdlib.h>

// Function to compute Gray encoding
int gray_encode(int binary)
{
    return binary ^ (binary >> 1);
}

// Function to compute M-ary QAM modulation with Gray encoding
void myqammod(int *inputSignal, int inputSignalLength, int M, double *realPart, double *imagPart)
{
    int bitsPerSymbol = log2(M); // Number of bits per symbol
    int M_sqrt = sqrt(M);
    int grayCodedSymbol;
    int i;

    for (i = 0; i < inputSignalLength; i++)
    {
        // Gray encode the input
        grayCodedSymbol = gray_encode(inputSignal[i]);

        // Now compute the QAM symbol
        realPart[i] = (2.0 * (grayCodedSymbol % M_sqrt) - M_sqrt + 1);     // I
        imagPart[i] = (2.0 * (grayCodedSymbol / M_sqrt) - M_sqrt + 1); // Q
    }
}

// This function will be called from MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *inputSignal;  // input array
    int M;                // M-ary QAM value
    double *outputSignalReal, *outputSignalImag; // output array
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
    plhs[0] = mxCreateDoubleMatrix(1, inputSignalLength, mxCOMPLEX); // Create a complex array

    // Assign pointers to each input and output
    inputSignal = mxGetPr(prhs[0]);
    outputSignalReal = mxGetPr(plhs[0]);
    outputSignalImag = mxGetPi(plhs[0]);

    // Convert input signal from double to int
    int* inputSignalInt = malloc(inputSignalLength * sizeof(int));
    for (i = 0; i < inputSignalLength; i++) {
        inputSignalInt[i] = (int)inputSignal[i];
    }

    // Call the QAM modulation function
    myqammod(inputSignalInt, inputSignalLength, M, outputSignalReal, outputSignalImag);

    free(inputSignalInt);
}
