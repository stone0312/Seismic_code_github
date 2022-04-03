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

//subroutine to solve the lower triangular linear equations with chasing-after algorithm
int matrix_pursuing(complex <float> **a, complex <float> *b, complex <float> *x, int n, float lamda1)
{
  int ix,ix1,iz;
  complex<float> tmp;
 
  for(ix=0;ix<n;ix++)
    for(iz=0;iz<n;iz++) 
       a[ix][iz].real()+=lamda1;

  for(iz=0;iz<n;iz++)
     x[iz]=(0.0,0.0);

  for(iz=0;iz<n;iz++)
    {
       tmp=(0.0,0.0);
       for(ix=0;ix<iz;ix++)
         tmp+=a[iz][ix]*x[ix];

       x[iz]=(b[iz]-tmp)/a[iz][iz];

    } 

  return 0;  

}

int main()
{
   int ix,iz;

   int n=3;
   complex<float> **a;
   a=alloc2complex(n,n);
   complex<float> *b;
   b=alloc1complex(n);
   complex<float> *x;
   x=alloc1complex(n);

   complex<float> *b1;
   b1=alloc1complex(n);

   float lamda1=0.01;

   for(ix=0;ix<n;ix++)
     for(iz=0;iz<n;iz++)
       a[ix][iz]=(0.0,0.0);

   a[0][0].real()=1; 
   a[0][0].imag()=2; 
   a[1][0].real()=2; 
   a[1][0].imag()=1; 
   a[1][1].real()=2; 
   a[1][1].imag()=3; 
   a[2][0].real()=1; 
   a[2][0].imag()=3; 
   a[2][1].real()=3; 
   a[2][1].imag()=4; 
   a[2][2].real()=4; 
   a[2][2].imag()=5; 

   b[0].real()=1;
   b[0].imag()=0;
   b[1].real()=1;
   b[1].imag()=0;
   b[2].real()=1;
   b[2].imag()=0;

   x[0]=b[0]/a[0][0];
 
   cout<<x[0]<<endl;

   matrix_pursuing(a,b,x,n,lamda1); 

   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<x[ix]<<endl;

   complex_matrix_multiply_vector(b1,a,x,n,n);

   cout<<"True Vector b is===="<<endl;
   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<b[ix]<<endl;

   cout<<"Calculated Vector b1 is===="<<endl;
   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<b1[ix]<<endl;

   return 0;

}























