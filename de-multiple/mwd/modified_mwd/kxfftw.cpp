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
  ns=2500;
  nr=1;
  nsr=ns*nr;
  lt=4000;
  ltt=2*lt;

  int ix,it;

  float **u;
  u=alloc2float(lt,ns);
 
  complex <float> *uk;
  uk=alloc1complex(lt);
 
  complex<float> **u1;
  u1=alloc2complex(lt,ns);

  ifstream swq1;
  swq1.open("/data1/swq/srme_data/orig_shots.dat",ios::binary);
  if(!swq1)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 
  for(ix=0;ix<ns;ix++)
     for(it=0;it<lt;it++)
       swq1.read((char*)&u[ix][it],sizeof(u[ix][it]));
  swq1.close();


  ofstream swq2;
  swq2.open("/data1/swq/srme_data/test_kxfftw.dat",ios::binary);
  if(!swq2)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 

    fftwf_complex *in1,*out1,*in2,*out2;
    fftwf_plan p1, p2, p3,p4;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);

    p1=fftwf_plan_dft_1d(ns,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ns,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);


   for(it=0;it<lt;it++)
     {
        if((in1==NULL)||(out1==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(ix=0;ix<ns;ix++)
                 {
                    in1[ix][0]=u[ix][it];
                    in1[ix][1]=0.0;
                 }
           }

        fftwf_execute(p1);

        for(ix=0;ix<ns;ix++)
          {
             uk[ix].real()=out1[ix][0];
             uk[ix].imag()=out1[ix][1];
          }
 
         if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(ix=0;ix<ns;ix++)
                 {
                    in2[ix][0]=uk[ix].real();
                    in2[ix][1]=uk[ix].imag();
                 }
           }

        fftwf_execute(p2);

        for(ix=0;ix<ns;ix++)
          {
             u1[ix][it].real()=out2[ix][0];
             u1[ix][it].imag()=out2[ix][1];
          }

        cout<<it<<" time slice done!"<<endl;
     }

     for(ix=0;ix<ns;ix++)
        for(it=0;it<lt;it++)
           swq2.write((char*)&u1[ix][it].real(),sizeof(u1[ix][it].real()));
     swq2.close();
  

     return 0;

}






































