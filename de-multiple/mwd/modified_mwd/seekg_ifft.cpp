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
  ns=479;
  nr=479;
  nsr=ns*nr;
  lt=3000;
  ltt=2*lt;

  int is,ir,it;
 
  complex<float> *u1;
  u1=alloc1complex(ltt);

  for(it=0;it<ltt;it++)
      u1[it]=(0.0,0.0); 

  float *mf1;
  mf1=alloc1float(ltt);

  ifstream swq2;
  swq2.open("/data2/swq/zj_data/479shot_fft_predicted_model_real_positive_impulse.dat",ios::binary);
  if(!swq2)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 

  ifstream swq22;
  swq22.open("/data2/swq/zj_data/479shot_fft_predicted_model_imag_positive_impulse.dat",ios::binary);
  if(!swq22)
    {
       cout<<"cannot open file"<<endl;
       abort();
    }

  ofstream swq4;
  swq4.open("/data2/swq/zj_data/check_multishots1.dat",ios::binary);
 
  fftwf_complex *in2,*out2;
  fftwf_plan p2;
  in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

  p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

   for(is=0;is<nsr;is++)
     {
        for(it=0;it<ltt;it++)
            u1[it]=(0.0,0.0);
        for(it=0;it<lt;it++)
            mf1[it]=0.0;

        swq2.seekg(is*4,ios::beg);
        swq2.read((char*)&u1[0].real(),sizeof(u1[0].real()));
        swq22.seekg(is*4,ios::beg);
        swq22.read((char*)&u1[0].imag(),sizeof(u1[0].imag()));

//        for(it=1;it<lt+1;it++)
        for(it=1;it<900;it++)
          {
            swq2.seekg((nsr-1)*4,ios::cur);
            swq22.seekg((nsr-1)*4,ios::cur);
            swq2.read((char*)&u1[it].real(),sizeof(u1[it].real()));
            swq22.read((char*)&u1[it].imag(),sizeof(u1[it].imag()));
          }
        for(it=ltt-899;it<ltt;it++)
          {
            u1[it].real()=u1[ltt-it].real();
            u1[it].imag()=-u1[ltt-it].imag();
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

         for(it=0;it<lt;it++)
            swq4.write((char*)&mf1[it],sizeof(mf1[it]));
        cout<<is<<" trace done!"<<endl;
    }
     swq2.close();
     swq4.close();
     swq22.close();

     return 0;

}






































