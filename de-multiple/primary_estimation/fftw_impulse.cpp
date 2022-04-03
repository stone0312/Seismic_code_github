using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"

#define pai 3.14159265


int main()
{
  int lt;
  lt=100;

  int it;
 
  float *u1;
  u1=alloc1float(lt);

  for(it=0;it<lt;it++)
    u1[it]=0.0;
  u1[50]=1.0;

  complex <float> *uf1;
  uf1=alloc1complex(lt);
  
  float *amp;
  amp=alloc1float(lt);

    fftwf_complex *in1,*out1;
    fftwf_plan p1;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

    p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
  
    if((in1==NULL)||(out1==NULL))
            cout<<"memory insufficient"<<endl;
    else
      {
        for(it=0;it<lt;it++)
          {
             in1[it][0]=u1[it];
             in1[it][1]=0.0;
          }
      }

    fftwf_execute(p1);
  
        for(it=0;it<lt;it++)
          {
             uf1[it].real()=out1[it][0];
             uf1[it].imag()=out1[it][1];
          }

        for(it=0;it<lt;it++)
            amp[it]=sqrt(pow(uf1[it].real(),2)+pow(uf1[it].imag(),2));


        for(it=0;it<lt;it++)
           cout<<it<<" ==== Real Part is ===="<<uf1[it].real()<<" ==== Imag Part is ===="<<uf1[it].imag()<<" ==== Amp is ===="<<amp[it]<<endl;


        cout<<"all done!"<<endl;

     return 0;

}






































