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
   ns=200;
   lt=4000;
   ltt=2*lt;

   int is,it;

   complex<float> **u;
   u=alloc2complex(ltt,ns);

   for(is=0;is<ns;is++)
      for(it=0;it<ltt;it++)
        u[is][it]=(0.0,0.0);

   float **amp;
   amp=alloc2float(ltt,ns);

   complex<float> **u1;
   u1=alloc2complex(ltt,ns);

   for(is=0;is<ns;is++)
      for(it=0;it<ltt;it++)
        u1[is][it]=0.0;

    fftwf_complex *in3,*out3,*in4,*out4;
    fftwf_plan p3,p4;

    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p3=fftwf_plan_dft_1d(ns,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(ltt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

    ifstream swq1;
    swq1.open("/data1/swq/srme/100th_fft_predicted_model_real_positive_2lt.dat",ios::binary);

    ifstream swq11;
    swq11.open("/data1/swq/srme/100th_fft_predicted_model_imag_positive_2lt.dat",ios::binary);

    ofstream swq2;
    swq2.open("/data1/swq/srme/test_2dfft.dat",ios::binary);
    
    for(is=0;is<ns;is++)
       {
          swq1.seekg(is*4,ios::beg);
          swq1.read((char*)&u[is][0].real(),sizeof(u[is][0].real()));
          swq11.seekg(is*4,ios::beg);
          swq11.read((char*)&u[is][0].imag(),sizeof(u[is][0].imag()));

          for(it=1;it<ltt/2+1;it++)
             {
               swq1.seekg((ns-1)*4,ios::cur);
               swq11.seekg((ns-1)*4,ios::cur);
               swq1.read((char*)&u[is][it].real(),sizeof(u[is][it].real()));
               swq11.read((char*)&u[is][it].imag(),sizeof(u[is][it].imag()));
             }
       }
    swq1.close();
    swq11.close();
 
    for(is=0;is<ns;is++)
      for(it=ltt/2+1;it<ltt;it++)
         {
            u[is][it].real()=u[ns-is-1][ltt-it].real();
            u[is][it].imag()=-u[ns-is-1][ltt-it].imag();
         }
/*
    for(it=0;it<ltt;it++)
       cout<<it<<"   "<<u[50][it].real()<<endl;
    return 0;
    ofstream swq22;
    swq22.open("/data1/swq/srme/test_read_real.dat",ios::binary);
    for(is=0;is<ns;is++)
      for(it=0;it<ltt;it++)
        {
           swq22.write((char*)&u[is][it].real(),sizeof(u[is][it].real()));
        }
    return 0;
*/
    for(it=0;it<ltt;it++)
      {
          if((in3==NULL)||(out3==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(is=0;is<ns;is++)
                  {
                     in3[is][0]=u[is][it].real();
                     in3[is][1]=u[is][it].imag();
                  }
               fftwf_execute(p3);
    
               for(is=0;is<ns;is++)
                  {
                     u[is][it].real()=out3[is][0];
                     u[is][it].imag()=out3[is][1];
                  } 
            }
      }
     for(is=0;is<ns;is++)
       {
          if((in4==NULL)||(out4==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(it=0;it<ltt;it++)
                  {
                     in4[it][0]=u[is][it].real();
                     in4[it][1]=u[is][it].imag();
                  }
               fftwf_execute(p4);

               for(it=0;it<ltt;it++)
                  {
                     u1[is][it].real()=out4[it][0];
                     u1[is][it].imag()=out4[it][1];
                  }
            }
       }

    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq2.write((char*)&u1[is][it].real(),sizeof(u1[is][it].real()));
    swq2.close();

    cout<<"ALL DONE!"<<endl;
   
    return 0;


}



















































