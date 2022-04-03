
using namespace std;
#include <iostream>
#include "math.h"
#include <fstream>
#include "stdlib.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"

int main()
{
  int ns,nr,nsr,lt,ltt;
  ns=200;
  nr=200;
  nsr=ns*nr;
  lt=1000;
  ltt=2*lt;

  nsr=nsr*30;

  int ix,ir,it;
 
  float *u1;
  u1=alloc1float(ltt);

  complex <float> *mf1;
  mf1=alloc1complex(ltt);
  
  float *amp;
  amp=alloc1float(ltt);

  ifstream swq2;
  swq2.open("/data/swq/shallow_mod_for_paper/crg_200shots_sr0_wd100m_raw_pri_ftp_w40.dat",ios::binary);
  if(!swq2)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 
  ofstream swq4;
  swq4.open("/data/swq/shallow_mod_for_paper/crg_green_tau_p_frequency_domain_real.dat",ios::binary);
 
  ofstream swq44;
  swq44.open("/data/swq/shallow_mod_for_paper/crg_green_tau_p_frequency_domain_imag.dat",ios::binary);
/*    
  ofstream swq444;
  swq444.open("/home/swq/20130902/fft_zcj_amp.dat",ios::binary);
*/
    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

   
   for(ix=0;ix<nsr;ix++)
     {
          for(it=0;it<ltt;it++)
             u1[it]=0.0;

          for(it=0;it<lt;it++)
           swq2.read((char*)&u1[it],sizeof(u1[it]));
         
         if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(it=0;it<ltt;it++)
                 {
                    in2[it][0]=u1[it];
                    in2[it][1]=0.0;
                 }
           }

        fftwf_execute(p2);

        for(it=0;it<ltt;it++)
          {
             mf1[it].real()=out2[it][0]/ltt;
             mf1[it].imag()=out2[it][1]/ltt;
          }

        for(it=0;it<ltt;it++)
            amp[it]=sqrt(pow(mf1[it].real(),2)+pow(mf1[it].imag(),2));

         for(it=0;it<ltt;it++)
          {
            swq4.write((char*)&mf1[it].real(),sizeof(mf1[it].real()));
            swq44.write((char*)&mf1[it].imag(),sizeof(mf1[it].imag()));
//            swq444.write((char*)&amp[it],sizeof(amp[it]));
          }
        if(ix%100==0)
            cout<<ix+1<<" trace done!"<<endl;  
    }

        cout<<"all done!"<<endl;

     swq2.close();
     swq4.close();
     swq44.close();
//     swq444.close();

     return 0;

}






































