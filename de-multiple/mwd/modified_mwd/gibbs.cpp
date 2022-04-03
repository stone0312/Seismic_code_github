#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define cut 500
int main()
{
  int ns,nr,nsr,lt,ltt;
  ns=479;
  nr=479;
  nsr=ns*nr;
  lt=3000;
  ltt=2*lt;

  int ix,ir,it;
 
  complex<float> **u1;
  u1=alloc2complex(ltt,ns);

  for(ix=0;ix<ns;ix++)
    for(it=0;it<ltt;it++)
      u1[ix][it]=(0.0,0.0); 

  float *mf1;
  mf1=alloc1float(ltt);

  ifstream swq2;
  swq2.open("/data2/swq/zj_data/479shot_raw_real.dat",ios::binary);
  if(!swq2)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 

  ifstream swq22;
  swq22.open("/data2/swq/zj_data/479shot_raw_imag.dat",ios::binary);
  if(!swq22)
    {
       cout<<"cannot open file"<<endl;
       abort();
    }

  for(ix=0;ix<ns;ix++)
    for(it=0;it<ltt;it++)
       {
            swq2.read((char*)&u1[ix][it].real(),sizeof(u1[ix][it].real()));
            swq22.read((char*)&u1[ix][it].imag(),sizeof(u1[ix][it].imag()));
       }

  for(ix=0;ix<ns;ix++)
    for(it=cut;it<ltt-cut;it++)
        u1[ix][it]=(0.0,0.0);

  ofstream swq4;
  swq4.open("/data2/swq/zj_data/test_gibbs_1st_shot_cut500.dat",ios::binary);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
   
   for(ix=0;ix<ns;ix++)
     {
          for(it=0;it<ltt;it++)
           swq2.read((char*)&u1[ir][it],sizeof(u1[ir][it]));
         
         if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(it=0;it<ltt;it++)
                 {
                    in2[it][0]=u1[ix][it].real();
                    in2[it][1]=u1[ix][it].imag();
                 }
           }

        fftwf_execute(p2);

        for(it=0;it<ltt;it++)
           mf1[it]=out2[it][0];

         for(it=0;it<lt;it++)
            swq4.write((char*)&mf1[it],sizeof(mf1[it]));
        cout<<ix<<" trace done!"<<endl;
    }
     swq2.close();
     swq4.close();
     swq22.close();

     return 0;

}






































