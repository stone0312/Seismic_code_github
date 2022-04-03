using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#define pai 3.14159265

//subroutine a=b*c, a is m dimision vector, b is m*n matrix, c is n dimision vector
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

//a=b*c m*n,n*p--->m*p
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
int complex_CG(complex <float> **a, complex <float> *x1, complex <float> *x2,complex <float> *f,complex <float> *p1,complex <float> *p2,complex <float> *r1,complex <float> *r2,complex <float> *ap, complex <float> *apcong,int np, float err)
{
  int ip,ix,iz;

  float alfa1,alfa2,beta1,beta2;
  float r1norm1,r2norm1,r1norm,p1ap,r2norm;
  complex<float> p1ap1_tmp;

  float norm2min=0.00000001;
  int iter_max=30;

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

    }

  return 0;

}

int main()
{
    int m,n,im,in;
    m=3;
    n=9;    
 
    float err=0.00001;

    complex<float> **a;
    a=alloc2complex(n,m);
    complex<float> **at;
    at=alloc2complex(m,n);
    complex<float> **ata;
    ata=alloc2complex(n,n);
    complex<float> *atb;
    atb=alloc1complex(n);
    complex<float> *b;
    b=alloc1complex(m);
    complex<float> *b1;
    b1=alloc1complex(m);
 
   complex<float> *solu;
   solu=alloc1complex(n);
   complex<float> *x1;
   x1=alloc1complex(n);
   complex<float> *p11;
   p11=alloc1complex(n);
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
 
   for(im=0;im<m;im++)
     for(in=0;in<n;in++)
        a[im][in]=(0.0,0.0);

   a[0][0].real()=1.0;
   a[0][0].imag()=0.0;
   a[0][1].real()=2.0;
   a[0][1].imag()=1.0;
   a[0][2].real()=-2.0;
   a[0][2].imag()=3.0;
   a[1][3].real()=2.0;
   a[1][3].imag()=-2.0;
   a[1][4].real()=4.0;
   a[1][4].imag()=5.0;
   a[1][5].real()=3.0;
   a[1][5].imag()=6.0;
   a[2][6].real()=-7.0;
   a[2][6].imag()=8.0;
   a[2][7].real()=4.0;
   a[2][7].imag()=5.0;
   a[2][8].real()=1.0;
   a[2][8].imag()=-8.0;

   b[0].real()=1.0;
   b[0].imag()=2.0;
   b[1].real()=8.0;
   b[1].imag()=9.0;
   b[2].real()=11.0;
   b[2].imag()=0.0;

   cout<<"Matrix A is===="<<endl;
   for(im=0;im<m;im++)
    {
     for(in=0;in<n;in++)
       cout<<a[im][in]<<"   ";
     cout<<endl;
    }

   cout<<"Vector b is===="<<endl;
   for(im=0;im<m;im++)
       cout<<b[im]<<"   ";
   cout<<endl;

   for(in=0;in<n;in++)
     for(im=0;im<m;im++)
      {
        at[in][im].real()=a[im][in].real();
        at[in][im].imag()=-a[im][in].imag();
      }

   complex_matrix_multiply_vector(atb, at, b, n, m);

   complex_matrix_multiply_matrix(ata, at, a, n, m, n);

/*
   ata[0][0].real()=1.0;
   ata[0][0].imag()=0.0;
   ata[0][1].real()=4.0;
   ata[0][1].imag()=5.0;
   ata[1][0].real()=4.0;
   ata[1][0].imag()=-5.0;
   ata[1][1].real()=1.0;
   ata[1][1].imag()=0.0;

   atb[0].real()=1.0;
   atb[0].imag()=3.0;
   atb[1].real()=1.0;
   atb[1].imag()=0.0;
*/

   cout<<"Matrix ATA is===="<<endl;
   for(im=0;im<n;im++)
    {
     for(in=0;in<n;in++)
       cout<<ata[im][in]<<"   ";
     cout<<endl;
    }

   cout<<"Vector ATB is===="<<endl;
   for(im=0;im<n;im++)
       cout<<atb[im]<<"   ";
   cout<<endl;

   complex_CG(ata, x1, solu, atb, p11, p2, r1, r2, ap, apcong, n, err); 

  for(in=0;in<n;in++)
     cout<<in<<"   "<<solu[in]<<endl; 
  cout<<"================================================================"<<endl;

   complex_matrix_multiply_vector(b1, a, solu, m, n);
  
  for(in=0;in<m;in++)
     cout<<in<<"   "<<b1[in]<<endl; 

  return 0;

}











