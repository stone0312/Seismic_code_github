#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"

int main()
{
  int ns,lt,ltt;
  ns=2500;
  lt=4000;
  ltt=2*lt;

  int ix,it;
 
  complex<float> **u;
  u=alloc2complex(ltt,ns);

  complex<float> **u1;
  u1=alloc2complex(ltt,ns);

  for(ix=0;ix<ns;ix++)
    for(it=0;it<ltt;it++)
     {
       u[ix][it].real()=0.0;
       u[ix][it].imag()=0.0;
       u1[ix][it].real()=0.0;
       u1[ix][it].imag()=0.0;
     }
   
  ifstream swq1;
  swq1.open("/data2/swq/kx_omega_domain_green_function_real.dat",ios::binary);

  ifstream swq2;
  swq2.open("/data2/swq/kx_omega_domain_green_funciton_imag.dat",ios::binary);

  ofstream swq3;
  swq3.open("/data2/swq/x_omega_domain_green_function_real1.dat",ios::binary);

  ofstream swq4;
  swq4.open("/data2/swq/x_omega_domain_green_function_imag1.dat",ios::binary);
  

  for(ix=0;ix<ns;ix++)
     {
        for(it=0;it<ltt;it++)
           {
               swq1.read((char*)&u[ix][it].real(),sizeof(u[ix][it].real()));
               swq2.read((char*)&u[ix][it].imag(),sizeof(u[ix][it].imag()));
           }
     } 
    
    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);

    p2=fftwf_plan_dft_1d(ns,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

    for(it=0;it<ltt;it++)
      {
         if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(ix=0;ix<ns;ix++)
                 {
                    in2[ix][0]=u[ix][it].real();
                    in2[ix][1]=u[ix][it].imag();
                 }
           }

        fftwf_execute(p2);

        for(ix=0;ix<ns;ix++)
          {
             u1[ix][it].real()=out2[ix][0];
             u1[ix][it].imag()=out2[ix][1];
          }
      }

   for(ix=0;ix<ns;ix++)
     for(it=0;it<ltt;it++)
       {
          swq3.write((char*)&u1[ix][it].real(),sizeof(u1[ix][it].real()));
          swq4.write((char*)&u1[ix][it].imag(),sizeof(u1[ix][it].imag()));
       } 

    return 0;

}






































