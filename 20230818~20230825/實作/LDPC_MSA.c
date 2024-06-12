#include "mex.h"
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <limits.h>
#include "matrix.h" 

#define	LLR_IN		prhs[0]
#define ITER_IN		prhs[1]
#define Hrow_IN		prhs[2]//row table
#define RsizeIN		prhs[3]//row table size
#define LDPC_num_IN prhs[4]
#define C_attenIN   prhs[5]
#define LLR2_OUT	plhs[0]

typedef struct{
    float *real;
}Complex2;
Complex2 LLRm,Hrow,Hrow_size,LLR2;
float ITER_num;
float C_atten;
float LDPC_num;

int H_row_master[ 648][ 8] ={};
int H_row_master_size[ 648]={};
float  LLR[1296][5464]={};
float    R[ 648][1296]={};
float    S[   8]	   ={};

int LLRPos(int row,int col){	//1296*5480
	return row + 1296*col; 
}
float sign(float ram){
    return (float)((ram>0)-(ram<0));
}
void init_LLR2(){
	int i,j;
	for(i=0;i<1296;i++)
		for(j=0;j<LDPC_num;j++)
			LLR2.real[LLRPos(i,j)] = LLRm.real[LLRPos(i,j)];
}
void init_Hpos(){
	int i,j;
	for(i=0;i<648;i++){
		H_row_master_size[i] = Hrow_size.real[i];
		for(j=0;j<8;j++)
			H_row_master[i][j] = Hrow.real[i+j*648] - 1;
	}
}
void init_R(){
	int row,col;
	for(row=0;row<648;row++)
		for(col=0;col<1296;col++)
			R[row][col] = 0;
}
void init_S(int row,int index){
	int i;
	for(i=0;i<H_row_master_size[row];i++){
		S[i] 	= LLR2.real[ LLRPos(H_row_master[row][i] ,index) ] - R[row][ H_row_master[row][i] ];
	}
}
void update(int index){
	int row,col,kcol,i;
	float product,mini,ram;
	for(row=0;row<648;row++){
		init_S(row,index);
		for(col=0;col<H_row_master_size[row];col++){
			product = 1;
			ram    = S[col];
			S[col] = DBL_MAX;
			mini   = DBL_MAX;
			for(kcol = 0;kcol<H_row_master_size[row];kcol++){
				product = product * sign( S[kcol] );
				if( fabs(S[kcol])<mini )
					mini = fabs(S[kcol]);
			}
			S[col] = ram;
			R[row][ H_row_master[row][col] ] = product*mini*C_atten;
			LLR2.real[ LLRPos( H_row_master[row][col] , index ) ] = S[col] + R[row][H_row_master[row][col]];
		}
	}
}
void culcal(){
	int index,ITER;
	init_Hpos();
	for(index=0;index<LDPC_num;index++){
		init_R();
		for(ITER=0;ITER < (int)ITER_num ;ITER++)//疊代
			update(index);
	}
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
	/* Ensure the input is single precision */
    if (!mxIsSingle(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notSingle",
                          "Input matrix must be type single.");
    }
	LLRm.real 		= mxGetSingles( LLR_IN      );
	ITER_num  		=*mxGetSingles( ITER_IN     );
	Hrow.real 		= mxGetSingles( Hrow_IN     );
	Hrow_size.real 	= mxGetSingles( RsizeIN 	   );
	LDPC_num	    =*mxGetSingles( LDPC_num_IN );
	C_atten			=*mxGetSingles( C_attenIN   );
	mwSize dims[2] 	= {1296,LDPC_num};
	LLR2_OUT = mxCreateNumericArray(2,dims,mxSINGLE_CLASS,mxREAL);
	LLR2.real= mxGetSingles( LLR2_OUT );
	init_LLR2();
	culcal();
}