#include "mex.h"
#include <complex.h>
#include <math.h>
#include "matrix.h" 

#define	 Y_IN		prhs[0]
#define	 H_IN		prhs[1]
#define TX_IN		prhs[2]
#define LMMSE_OUT	plhs[0]

typedef struct{
    double *real;
    double *imag;
	int		row;
	int		col;
}Complex2;
Complex2 LMMSE,Y,H;

int i,j,M,Tx,row,col,ram;

int YPos(int sc,int slot,int row){
	return sc + 1644*slot + 1644*560*row;
}
int HPos(int sc,int slot,int row,int col){
	return sc + 1644*slot + 1644*560*row + 1644*560*Tx*col;
}
void matrix_muilt(int Arow,int Acol,double complex A[Arow][Acol],int Brow,int Bcol,double complex B[Brow][Bcol],double complex Ans[Arow][Bcol]){
	int i,j,k;
	double complex ram;
	for(i=0;i<Arow;++i)
		for(k=0;k<Acol;++k){
			ram = A[i][k];
			for(j=0;j<Bcol;++j)
				Ans[i][j] += ram*B[k][j];
		}
}
void herm   (double complex output[Tx][Tx],double complex input[Tx][Tx]){
    for(i=0;i<Tx;i++){
        for(j=0;j<Tx;j++)
            output[i][j] = conj(input[j][i]);
    }
}
void multip1(double complex output[Tx]    ,double complex input1[Tx][Tx],double complex input2[Tx]){
	for(i=0;i<Tx;i++)
	    output[i] = 0;
	for(i=0;i<Tx;i++)
	    for(ram = 0;ram<Tx;ram++)
            output[i] += input1[i][ram] * input2[ram];
}
void multip2(double complex output[Tx][Tx],double complex input1[Tx][Tx],double complex input2[Tx][Tx]){
    for(i=0;i<Tx;i++)
        for(j=0;j<Tx;j++)
			output[i][j] = 0;
	for(i=0;i<Tx;i++)
        for(j=0;j<Tx;j++)
            for(ram=0;ram<Tx;ram++)
                output[i][j] += input1[i][ram] * input2[ram][j];
}
void matadd (double complex output[Tx][Tx],double complex input1[Tx][Tx],double complex input2[Tx][Tx]){
	for(i=0;i<Tx;i++)
		for(j=0;j<Tx;j++)
			output[i][j] = input1[i][j] + input2[i][j];
}
void inv   (double complex output[2][2],double complex input[2][2]){
	double complex det = input[0][0] * input[1][1] - input[1][0] * input[0][1];
	output[0][0] = input[1][1]/det;
	output[0][1] =-input[0][1]/det;
	output[1][0] =-input[0][1]/det;
	output[1][1] = input[0][0]/det;
}
void culcal(){
	int SC,slot;
	int i,j;
	double complex NoM[Tx][Tx];
	double complex h  [Tx][Tx];
	double complex hh [Tx][Tx];
	double complex mul[Tx][Tx];
	double complex inh[Tx][Tx];
	double complex y  [Tx]	  ;
	double complex x  [Tx]	  ;
	double complex tep[Tx]	  ;
	double complex te2[Tx][Tx];
	//init No
	for(i=0;i<Tx;i++)
		for(j=0;j<Tx;j++){
			if(i==j)
				NoM[i][j] = 1;
			else
				NoM[i][j] = 0;
		}
			
	//start calcul
	for(SC=0;SC<1644;SC++){
		for(slot=0;slot<560;slot++){
			for(i=0;i<Tx;i++){
				y[i]		= Y.real[YPos(SC,slot,i  )] + I*Y.imag[YPos(SC,slot,i)];
				for(j=0;j<Tx;j++)
					h[i][j]	= H.real[HPos(SC,slot,i,j)] + I*H.imag[HPos(SC,slot,i,j)];
			}
			herm	(hh	,h);
			multip2	(mul,hh	,h	);
			matadd	(te2,mul,NoM);
			inv		(inh,te2	);
			multip1 (tep,hh	,y	);
			multip1 ( x ,inh,tep);
			for(i=0;i<Tx;i++){
				LMMSE.real[YPos(SC,slot,i)] = creal(x[i]);
				LMMSE.imag[YPos(SC,slot,i)] = cimag(x[i]);
			}
		}
	}
}

//Y=1644*560*4 ; H = 1644*560*4*4; ZF=1644*560*4
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	Y.real  = mxGetPr(  Y_IN  );
	Y.imag  = mxGetPi(  Y_IN  );
	Y.row	= mxGetM (  Y_IN  );
	Y.col	= mxGetN (  Y_IN  );
	H.real  = mxGetPr(  H_IN  );
	H.imag  = mxGetPi(  H_IN  );
	H.row	= mxGetM (  H_IN  );
	H.col	= mxGetN (  H_IN  );
	Tx 		=*mxGetPr( TX_IN  );
	M 		=*mxGetPr( TX_IN  );
	mwSize dims[3] = {1644,560,Tx};
	LMMSE_OUT	= mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxCOMPLEX);
	LMMSE.real = mxGetPr( LMMSE_OUT );
	LMMSE.imag = mxGetPi( LMMSE_OUT );
	if(Tx==2)
		culcal();
	else
		mexErrMsgTxt("only support in 2x2 MIMO");
}