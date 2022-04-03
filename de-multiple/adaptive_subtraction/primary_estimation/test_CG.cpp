#include "iostream.h"
#include "math.h"
#include "stdio.h"
#include "fstream.h"
#include "stdlib.h"
#include "alloc.c"
    
extern "C" {
void real_CG (float *A,float *x1,float *x2,float *f,float *P1,float *P2,float *R1,float *R2,float *AP,int *np,float *err);
void real_CG_ (float *A,float *x1,float *x2,float *f,float *P1,float *P2,float *R1,float *R2,float *AP,int *np,float *err);
void real_cg (float *A,float *x1,float *x2,float *f,float *P1,float *P2,float *R1,float *R2,float *AP,int *np,float *err);
void real_cg_ (float *A,float *x1,float *x2,float *f,float *P1,float *P2,float *R1,float *R2,float *AP,int *np,float *err);
}   

int main()
{
    int i,j,win;
    win=2;
    float err=0.001;
    float **A;
    A=alloc2float(win,win);
    float **AT;
    AT=alloc2float(win,win);  

    float *f=new float [win];
    float *b=new float [win];

    float *P1=new float [win];
    float *P2=new float [win];
    float *R1=new float [win];
    float *R2=new float [win];
    float *X1=new float [win];
    float *X2=new float [win];
    float *AP=new float [win];

    A[0][0]=2;
    A[0][1]=7;
    A[1][0]=7;
    A[1][1]=2;

    cout<<"Matrix A is==="<<endl;
    for(i=0;i<win;i++)
      {
         for(j=0;j<win;j++)
            cout<<A[i][j]<<"  ";
         cout<<endl;
      }

    for(i=0;i<win;i++)
       for(j=0;j<win;j++)
          AT[i][j]=A[i][j];

    b[0]=1;
    b[1]=1;
   
 
 
    cout<<"Matrix b is==="<<endl;
    for(i=0;i<win;i++)
       cout<<b[i]<<"  ";
    cout<<endl;
   
    real_cg_ (A[0], X1, X2, b, P1, P2, R1, R2, AP, &win, &err);

    cout<<"Solution is==="<<endl;
    cout<<X2[0]<<"  "<<X2[1]<<endl;

    return 0;
}
