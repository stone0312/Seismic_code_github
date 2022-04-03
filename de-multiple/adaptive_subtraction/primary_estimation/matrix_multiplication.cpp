using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#define pai 3.14159265

int complex_matrix_multiply_vector(complex<float> *a, complex<float>**b, complex<float>*c, int m, int n)
{
   int ix,iz;
   for(ix=0;ix<m;ix++)
      a[ix]=(0.0,0.0);

   for(ix=0;ix<m;ix++)
     {
       for(iz=0;iz<n;iz++)
         a[ix]+=b[ix][iz]*c[iz];
     }

   return 0;

}

int main()
{
   int ir1,ir2,ix,iz;

   int n=3;
   complex<float> **a;
   a=alloc2complex(n,n);
   complex<float> *b;
   b=alloc1complex(n);
   complex<float> *x;
   x=alloc1complex(n);

   complex<float> **a1;
   a1=alloc2complex(n,n);
   complex<float> *b1;
   b1=alloc1complex(n);
   complex<float> *x1;
   x1=alloc1complex(n);

   for(ix=0;ix<n;ix++)
     for(iz=0;iz<n;iz++)
       a[ix][iz]=(0.0,0.0);

   a[0][0].real()=1; 
   a[0][0].imag()=1; 
   a[1][0].real()=2; 
   a[1][0].imag()=1; 
   a[1][1].real()=3; 
   a[1][1].imag()=1; 
   a[2][0].real()=5; 
   a[2][0].imag()=1; 
   a[2][1].real()=6; 
   a[2][1].imag()=1; 
   a[2][2].real()=7; 
   a[2][2].imag()=1; 

   b[0].real()=2;
   b[0].imag()=1;
   b[1].real()=5;
   b[1].imag()=1;
   b[2].real()=1;
   b[2].imag()=1;

   complex_matrix_multiply_vector(x,a,b,n,n);

   cout<<"True Vector b is===="<<endl;
   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<x[ix]<<endl;

   b1[0].real()=1;
   b1[0].imag()=0;
   b1[1].real()=1;
   b1[1].imag()=0;
   b1[2].real()=1;
   b1[2].imag()=0;

   for(ir1=0;ir1<n;ir1++)
     {
       for(ir2=0;ir2<n;ir2++)
         a1[ir2][ir1]=a[ir2][ir1]*b[ir1];
     }

    complex_matrix_multiply_vector(x1,a1,b1,n,n);

   cout<<"Calculated Vector b1 is===="<<endl;
   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<x1[ix]<<endl;

   return 0;

}























