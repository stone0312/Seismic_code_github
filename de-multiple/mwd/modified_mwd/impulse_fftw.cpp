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
  ns=200;
  nr=200;
  nsr=ns*nr;
  lt=4000;
  ltt=2*lt;

  int ix,ir,it;
 
    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

    float *a;
    a=alloc1float(ltt);

    complex<float> *aa;
    aa=alloc1complex(ltt);

    for(it=0;it<ltt;it++) 
       a[it]=0.0;
    a[0]=1.0;

    for(it=0;it<ltt;it++)  
       {
          in2[it][0]=a[it];
          in2[it][1]=0.0;
       }

    fftwf_execute(p2);

    for(it=0;it<ltt;it++)
       {
          aa[it].real()=out2[it][0];
          aa[it].imag()=out2[it][1];
       }
     
    for(it=0;it<ltt;it++)
       cout<<it<<"  "<<aa[it]<<endl;



     return 0;

}






































