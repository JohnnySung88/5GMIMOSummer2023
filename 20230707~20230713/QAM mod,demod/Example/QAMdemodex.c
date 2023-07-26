#include "mex.h"
#include <math.h>

// Function to compute M-ary QAM demodulation with Gray decoding
void myqamdemod(double *inputSignal, int inputSignalLength, int M, double *outputSignal)
{
    int bitsPerSymbol = log2(M); // Number of bits per symbol
    int grayCodedSymbol;
    int i, j, bit, symbol;

    for (i = 0; i < inputSignalLength; i += 2)
    {
        // Convert I and Q back to Gray coded symbol
        grayCodedSymbol = ((int)round(inputSignal[i] * sqrt(M) + M - 1) / 2) +
                          ((int)round(inputSignal[i + 1] * sqrt(M) + M - 1) / 2) * ((int)sqrt(M));

        // Now decode the Gray code
        symbol = grayCodedSymbol ^ (grayCodedSymbol >> 1);

        // Assign the symbol to the output
        outputSignal[i / 2] = (double)symbol;
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
    plhs[0] = mxCreateDoubleMatrix(1, inputSignalLength / 2, mxREAL); // Half the length for demodulated symbols

    // Assign pointers to each input and output
    inputSignal = mxGetPr(prhs[0]);
    outputSignal = mxGetPr(plhs[0]);

    // Call the QAM demodulation function
    myqamdemod(inputSignal, inputSignalLength, M, outputSignal);
}
