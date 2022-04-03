#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"

int main()
{
  int ns,nr,nsr,lt,ltt;
  ns=1;
  nr=1;
  nsr=ns*nr;
  lt=4000;
  ltt=1*lt;

  int ix,ir,it;
 
  complex<float> *u1;
  u1=alloc1complex(ltt);

    for(it=0;it<ltt;it++)
      u1[it]=(0.0,0.0); 

  float *mf1;
  mf1=alloc1float(ltt);

  float *amp;
  amp=alloc1float(ltt);

  ofstream swq4;
  swq4.open("/home/swq/usrmge/smooth20_1st_trace_deghosting_100m_lamda100.dat",ios::binary);

  ofstream swq5;
  swq5.open("/home/swq/usrmge/smooth20_1st_trace_deghosting_100m_lamda100_amp.dat",ios::binary);


    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

  ifstream swq2;
  swq2.open("/home/swq/usrmge/smooth20_1st_trace_deghosting_100m_lamda100_real.dat",ios::binary);
  if(!swq2)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 

  ifstream swq22;
  swq22.open("/home/swq/usrmge/smooth20_1st_trace_deghosting_100m_lamda100_imag.dat",ios::binary);
  if(!swq22)
    {
       cout<<"cannot open file"<<endl;
       abort();
    }

  for(ix=0;ix<nsr;ix++)
    {
      for(it=0;it<ltt;it++)
       {
            swq2.read((char*)&u1[it].real(),sizeof(u1[it].real()));
            swq22.read((char*)&u1[it].imag(),sizeof(u1[it].imag()));
       }
    


      for(it=0;it<ltt;it++)
        {
           amp[it]=sqrt(u1[it].real()*u1[it].real()+u1[it].imag()*u1[it].imag()); 
           swq5.write((char*)&amp[it],sizeof(amp[it]));          
        }

     
         if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(it=0;it<ltt;it++)
                 {
                    in2[it][0]=u1[it].real();
                    in2[it][1]=u1[it].imag();
                 }
           }

        fftwf_execute(p2);

        for(it=0;it<ltt;it++)
           mf1[it]=out2[it][0]/ltt;

         for(it=0;it<ltt;it++)
            swq4.write((char*)&mf1[it],sizeof(mf1[it]));

         if((ix+1)%1000==0)
          cout<<ix+1<<" trace ifftw done!"<<endl;
    }
     swq2.close();
     swq4.close();
     swq22.close();

     return 0;

}






































