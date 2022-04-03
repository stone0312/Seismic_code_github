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

int complex_matrix_multiply_matrix(complex<float> **a, complex<float>**b, complex<float>**c, int m, int n, int p)
{
   int ix,iy,iz;
   for(ix=0;ix<m;ix++)
     for(iy=0;iy<p;iy++)
       a[ix][iy]=(0.0,0.0);

   for(ix=0;ix<m;ix++)
    {
      for(iy=0;iy<p;iy++)
        {
          for(iz=0;iz<n;iz++)
            a[ix][iy]+=b[ix][iz]*c[iz][iy];
        }
    }

  return 0;

}



int complex_vector_dotproduct(complex<float> *p1,complex<float> *ap,complex<float> *apcong,int np,float p1ap)
{  
   int ip;
   complex<float> p1ap1_tmp;
   p1ap1_tmp=(0.0,0.0);

   for(ip=0;ip<np;ip++)
     {
        apcong[ip].real()=ap[ip].real();
        apcong[ip].imag()=-ap[ip].imag();
     }

   for(ip=0;ip<np;ip++)
     p1ap1_tmp+=p1[ip]*apcong[ip];  

   p1ap=p1ap1_tmp.real();


   return 0;
}

//subroutine to solve the lower triangular linear equations with CG algorithm :AX2=f
int complex_CG(complex <float> **a, complex <float> *x1, complex <float> *x2,complex <float> *f,complex <float> *p1,complex <float> *p2,complex <float> *r1,complex <float> *r2,complex <float> *ap, complex <float> *apcong,int np, float err)
{
  int ip,ix,iz;
  
  float alfa1,alfa2,beta1,beta2;
  float r1norm1,r2norm1,r1norm,p1ap,r2norm;
  complex<float> p1ap1_tmp;
 
  float norm2min=0.000001;
  int iter_max=1000;

  for(ip=0;ip<np;ip++)
    {
      x1[ip].real()=0.0;
      x1[ip].imag()=0.0;
      x2[ip].real()=0.0;
      x2[ip].imag()=0.0;
      r1[ip].real()=0.0;
      r1[ip].imag()=0.0;
      r2[ip].real()=0.0;
      r2[ip].imag()=0.0;
      p1[ip].real()=0.0;
      p1[ip].imag()=0.0;
      p2[ip].real()=0.0;
      p2[ip].imag()=0.0;
    } 
 
  for(ip=0;ip<np;ip++)
    {
      r1[ip]=f[ip]; 
      p1[ip]=r1[ip]; 
    }

  int iteration=1;

  while(iteration<iter_max) 
    {
       r1norm1=0.0;
       for(ip=0;ip<np;ip++)
         r1norm1+=(r1[ip].real()*r1[ip].real()+r1[ip].imag()*r1[ip].imag());
       r1norm=r1norm1;
       
       complex_matrix_multiply_vector(ap,a,p1,np,np);        


       p1ap1_tmp=(0.0,0.0);

       for(ip=0;ip<np;ip++)
       {
         apcong[ip].real()=ap[ip].real();
         apcong[ip].imag()=-ap[ip].imag();
       }

      for(ip=0;ip<np;ip++)
        p1ap1_tmp+=p1[ip]*apcong[ip];  
      p1ap=p1ap1_tmp.real();


       if(fabs(p1ap)<norm2min)
         iteration=iter_max;

       else
        {
          alfa1=r1norm/p1ap;

          cout<<iteration<<" , "<<alfa1<<endl;

          for(ip=0;ip<np;ip++)
           {
             x2[ip].real()=x1[ip].real()+alfa1*p1[ip].real();
             x2[ip].imag()=x1[ip].imag()+alfa1*p1[ip].imag();
           }

          for(ip=0;ip<np;ip++)
           {
             r2[ip].real()=r1[ip].real()-alfa1*ap[ip].real();
             r2[ip].imag()=r1[ip].imag()-alfa1*ap[ip].imag();
           }

          r2norm1=0.0;
          for(ip=0;ip<np;ip++)
            r2norm1+=(r2[ip].real()*r2[ip].real()+r2[ip].imag()*r2[ip].imag());  
          r2norm=r2norm1;

          cout<<"Iteration Time===="<<iteration<<" , Error is==== "<<sqrt(r2norm)<<endl;

          beta2=r2norm/r1norm;     
          for(ip=0;ip<np;ip++)
            {
              p2[ip].real()=r2[ip].real()+beta2*p1[ip].real();
              p2[ip].imag()=r2[ip].imag()+beta2*p1[ip].imag();
            }
          for(ip=0;ip<np;ip++)
            {
              x1[ip]=x2[ip];
              r1[ip]=r2[ip];
              p1[ip]=p2[ip];
            } 

          if(sqrt(r2norm)<err)
             iteration=iter_max;
          else
             iteration+=1;

 
        }

          cout<<"222222====="<<iteration<<endl;
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
   complex<float> *x1;
   x1=alloc1complex(n);
   complex<float> *x2;
   x2=alloc1complex(n);
   complex<float> *p1;
   p1=alloc1complex(n);
   complex<float> *p2;
   p2=alloc1complex(n);
   complex<float> *r1;
   r1=alloc1complex(n);
   complex<float> *r2;
   r2=alloc1complex(n);
   complex<float> *ap;
   ap=alloc1complex(n);
   complex<float> *apcong;
   apcong=alloc1complex(n);

   complex<float> **at;
   at=alloc2complex(n,n);
   complex<float> **ata;
   ata=alloc2complex(n,n);
   complex<float> *atb;
   atb=alloc1complex(n);
   complex<float> *b1;
   b1=alloc1complex(n);

   float lamda1=0.01;

   for(ix=0;ix<n;ix++)
     for(iz=0;iz<n;iz++)
       a[ix][iz]=(0.0,0.0);

  a[0][0].real()=1.0; 
   a[0][0].imag()=1.0; 
   a[0][1].real()=3.0; 
   a[0][1].imag()=-2.0; 
   a[0][2].real()=2.0; 
   a[0][2].imag()=-2.0; 
   a[1][0].real()=3.0; 
   a[1][0].imag()=2.0; 
   a[1][1].real()=3.0; 
   a[1][1].imag()=2.0; 
   a[1][2].real()=5.0; 
   a[1][2].imag()=-2.0; 
   a[2][0].real()=2.0; 
   a[2][0].imag()=2.0; 
   a[2][1].real()=5.0; 
   a[2][1].imag()=2.0; 
   a[2][2].real()=4.0; 
   a[2][2].imag()=7.0; 

   for(ix=0;ix<n;ix++)
     a[ix][ix].real()+=0.1;

   for(ix=0;ix<n;ix++)
    {
      for(iz=0;iz<n;iz++)
        cout<<a[ix][iz]<<"   ";
      cout<<endl;
    }

   b[0].real()=6;
   b[0].imag()=-4;
   b[1].real()=22;
   b[1].imag()=0;
   b[2].real()=33;
   b[2].imag()=12;


   float err=0.00001;

   for(ix=0;ix<n;ix++)
     for(iz=0;iz<n;iz++)
      {
        at[ix][iz].real()=a[iz][ix].real();
        at[ix][iz].imag()=-a[iz][ix].imag();
      }

   complex_matrix_multiply_matrix(ata,at,a,n,n,n);

   complex_matrix_multiply_vector(atb,at,b,n,n);

   complex_CG(ata, x1, x2,atb,p1,p2,r1,r2,ap, apcong, n,err); 

   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<x2[ix]<<endl;

   complex_matrix_multiply_vector(b1,a,x2,n,n);

   cout<<"True Vector b is===="<<endl;
   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<b[ix]<<endl;

   cout<<"Calculated Vector b1 is===="<<endl;
   for(ix=0;ix<n;ix++)
      cout<<ix+1<<"==== "<<b1[ix]<<endl;

   return 0;

}























