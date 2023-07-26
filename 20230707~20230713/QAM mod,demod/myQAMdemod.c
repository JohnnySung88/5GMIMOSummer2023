#include "mex.h"
#include <math.h>
#include <stdlib.h>

// Function to compute Gray decoding
int gray_decode(int gray)
{
    int binary = 0;

    for (; gray; gray >>= 1)
    {
        binary ^= gray;
    }

    return binary;
}

// Function to compute M-ary QAM demodulation with Gray decoding
void myqamdemod(double *inputSignalReal, double *inputSignalImag, int inputSignalLength, int M, int *outputSignal)
{
    int bitsPerSymbol = log2(M); // Number of bits per symbol
    int M_sqrt = sqrt(M);
    int grayCodedSymbol;
    int i;

    for (i = 0; i < inputSignalLength; i++)
    {
        // compute the QAM symbol
        grayCodedSymbol = ((int)((inputSignalReal[i] + M_sqrt - 1) / 2.0)) % M_sqrt +
                          ((int)((inputSignalImag[i] + M_sqrt - 1) / 2.0)) * M_sqrt;

        // Gray decode the input
        outputSignal[i] = gray_decode(grayCodedSymbol);
    }
}

// This function will be called from MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *inputSignalReal, *inputSignalImag; // input array
    int M; // M-ary QAM value
    int *outputSignal; // output array
    int i, inputSignalLength;

    // Check for proper number of arguments
    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Two inputs required.");
    }

    // Check that both inputs are noncomplex double scalar values
    if (!mxIsDouble(prhs[0]) || !mxIsComplex(prhs[0]) ||
        mxGetNumberOfElements(prhs[1]) != 1 ||
        (M = (int)mxGetScalar(prhs[1])) != 4 && M != 16 && M != 64)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
                          "Input matrix must be noncomplex double and M must be 4, 16 or 64.");
    }

    // Get the length of the input array
    inputSignalLength = mxGetN(prhs[0]);

    // Create matrix for the return argument
    plhs[0] = mxCreateNumericMatrix(1, inputSignalLength, mxINT32_CLASS, mxREAL);

    // Assign pointers to each input and output
    inputSignalReal = mxGetPr(prhs[0]);
    inputSignalImag = mxGetPi(prhs[0]);
    outputSignal = (int*)mxGetData(plhs[0]);

    // Call the QAM demodulation function
    myqamdemod(inputSignalReal, inputSignalImag, inputSignalLength, M, outputSignal);
}
