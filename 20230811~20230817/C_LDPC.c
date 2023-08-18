#include"mex.h"
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
typedef struct{
double *real;
double *imag;
}Complex;

void LDPC(int xa,int xb,int ha,int hb,Complex X,Complex H,Complex output){
    Complex Eta,S,Lumda;
    
    Eta.real=malloc(ha*hb*sizeof(double));   //新增動態記憶體
//     idx.real=malloc(8*sizeof(double));
    S.real=malloc(8*sizeof(double));
    Lumda.real=malloc(xb*sizeof(double));
    int i,j,k,m,n,p,q,f,v;
    int iter;
    int idx_pl,Ssgn;
    double Spsi,Qtemp,Qtemppsi,Qtempsgn;
    int idx[8]={0};
    int Ssgn_temp[8]={0};
    double Little_num;
    iter=5;
    Little_num=1;
    for(i=0;i<25;i++)
    {Little_num*=0.1;}
    for(i=0;i<xb;i++)
    {Lumda.real[i]=X.real[i];}
    
    for(i=0;i<ha*hb;i++) //Eta初始化
        Eta.real[i]=0;
    for(i=0;i<8;i++) //idx初始化
        idx[i]=0;
    for(i=0;i<8;i++) //S初始化
        S.real[i]=0;
    
    
    for(i=0;i<iter;i++){ //疊代START
        
        for(j=0;j<ha;j++){ //LDPC.par_bits 
            
            idx_pl=0;            
            for(k=0;k<hb;k++){  //idx
                if(H.real[j+k*324]==1)
                {idx[idx_pl]=k;idx_pl++;}}//CHECK
            
            Ssgn=1;Spsi=0;
            for(f=0;f<idx_pl;f++){
                S.real[f]=Lumda.real[idx[f]]-Eta.real[j+idx[f]*324];//<<<<<<<<<<
                if(S.real[f]>0)
                {Ssgn_temp[f]=1;}
                else
                {Ssgn_temp[f]=-1;}
                Ssgn*=Ssgn_temp[f];
                Spsi+=-log(Little_num+tanh(fabs(S.real[f]*0.5)));}
            
            for(v=0;v<idx_pl;v++){
                Qtemp=Lumda.real[idx[v]]-Eta.real[j+idx[v]*324];//<<<<<<<<<<<<<<<
                Qtemppsi=-log(Little_num+tanh(fabs(S.real[v]*0.5)));
                if(Qtemp>0)
                    Qtempsgn=Ssgn;
                else
                    Qtempsgn=-Ssgn;
                Eta.real[j+idx[v]*324]=Qtempsgn*(-log(Little_num+tanh(fabs(Spsi-Qtemppsi)*0.5)));//<<<<<<<<<<<<<
                Lumda.real[idx[v]]=Qtemp+Eta.real[j+idx[v]*324];//<<<<<<<<<<<<<
            }
            
        }}
    for(i=0;i<ha;i++){
        if(Lumda.real[i]>0)
            output.real[i]=0;
        else
            output.real[i]=1;
    }
//     printf("%d\n",xa);//1
//     printf("%d\n",xb);//648
//     printf("%d\n",ha);//324
//     printf("%d\n",hb);//648
    
//     printf("%d\n",idx_pl);
//     printf("%d\n",Ssgn);
//     printf("%f\n",Qtemp);//<<<<<<<<<
//     printf("%f\n",Qtempsgn);
//     printf("%f\n",Eta.real[0]);//有問題啊你
//     printf("%f\n",Eta.real[32]);
//     printf("%f\n",Eta.real[75]);
//     printf("%f\n",Eta.real[106]);
//     printf("%f\n",Eta.real[112]);
    
//     printf("%f\n",Eta.real[324]);
//     printf("%f\n",Eta.real[324]);
//     printf("%f\n",Eta.real[324]);
//     printf("%f\n",Eta.real[324]);
//     printf("%f\n",Lumda.real[0]);
//     printf("%f\n",Lumda.real[1]);
//     printf("%f\n",Lumda.real[2]);
//     printf("%f\n",Lumda.real[3]);
//     printf("%f\n",Lumda.real[4]);
//     printf("%f\n",Lumda.real[5]);
//     printf("%f\n",Lumda.real[6]);
//     printf("%f\n",Lumda.real[638]);
//     printf("%f\n",Lumda.real[639]);
//     printf("%f\n",Lumda.real[640]);
//     printf("%f\n",Lumda.real[641]);
//     printf("%f\n",Lumda.real[642]);
//     printf("%f\n",Lumda.real[643]);
//     printf("%f\n",Lumda.real[644]);
//     printf("%f\n",Lumda.real[645]);
//     printf("%f\n",Lumda.real[646]);
//     printf("%f\n",Lumda.real[647]);
//     printf("%f\n",H.real[0]);
//     printf("%d\n",idx[0]);
//     printf("%d\n",idx[1]);
//     printf("%d\n",idx[2]);
//     printf("%d\n",idx[3]);
//     printf("%d\n",idx[4]);
//     printf("%d\n",idx[5]);
//     printf("%d\n",idx[6]);
//     printf("%f\n",S.real[0]);//<<<<<<<<<<<<<<<<
//     printf("%f\n",S.real[1]);
//     printf("%f\n",S.real[2]);
//     printf("%f\n",S.real[5]);
//     printf("%f\n",S.real[6]);//Check
//     printf("%d\n",Ssgn); //Check
//     printf("%f\n",Little_num);
//     printf("%f\n",Spsi);
//     printf("%d\n",idx_pl);
    free(Eta.real);         //釋放記憶體
    free(S.real);
    free(Lumda.real);

}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    Complex X;
    Complex H;
    Complex Output;
    int Xa,Xb,Ha,Hb;
//     double *Maxpos;
     
    Xa=mxGetM(prhs[0]);   
    Xb=mxGetN(prhs[0]);
    Ha=mxGetM(prhs[1]);
    Hb=mxGetN(prhs[1]);
    
    X.real=mxGetPr(prhs[0]);   
    H.real=mxGetPr(prhs[1]);
    
        
    plhs[0]=mxCreateDoubleMatrix(1,Ha,mxCOMPLEX);
    Output.real=mxGetPr(plhs[0]);
//     Acorr.imag=mxGetPi(plhs[0]);
//     plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
//     Maxpos=mxGetPr(plhs[0]);
   LDPC(Xa,Xb,Ha,Hb,X,H,Output);    
}
