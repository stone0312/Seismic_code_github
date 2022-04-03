
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
  char fn1[256],fn2[256],fn3[256];
  int ns, nr,lt, flag;
  float thre;

  ifstream swq;
  swq.open("hilbert_envelop.par");
  swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>flag>>thre;
  swq.close();


  float *u1;
  u1=alloc1float(lt);

  float *u2;
  u2=alloc1float(lt);

  float *u3;
  u3=alloc1float(lt);
  

  ifstream swq1;
    swq1.open(fn1,ios::binary);
    if(!swq1)
      {
        cout<<"cannot open "<<fn1<<endl;
        return 0;
      }

  ofstream swq2;
    swq2.open(fn2,ios::binary);
    if(!swq2)
      {
        cout<<"cannot open "<<fn2<<endl;
        return 0;
      }

  ofstream swq3;
    swq3.open(fn3,ios::binary);
    if(!swq3)
      {
        cout<<"cannot open "<<fn3<<endl;
        return 0;
      }

  int is, ir, it;

    fftwf_complex *trace_t,*trace_w;
    fftwf_plan p1, p2;
    trace_t=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    trace_w=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

    p1=fftwf_plan_dft_1d(lt,trace_t,trace_w,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(lt,trace_w,trace_t,FFTW_BACKWARD,FFTW_MEASURE);

 
  for(is=0;is<ns*nr;is++)
    {
      for(it=0;it<lt;it++)
        swq1.read((char*)&u1[it],sizeof(u1[it]));

       for(it=0;it<lt;it++)
          {
             trace_t[it][0]=u1[it];
             trace_t[it][1]=0.0;
          }

       fftwf_execute(p1);
 
                   
       for(it=0;it<lt/2;it++)
         {
             trace_w[it][0]*=2.0;
             trace_w[it][1]*=2.0;
         }

       for(it=lt/2;it<lt;it++)
         {
             trace_w[it][0]=0.0;
             trace_w[it][1]=0.0;
         }

       fftwf_execute(p2);
       
       for(it=0;it<lt;it++)
         u2[it]=sqrt(pow(trace_t[it][0],2)+pow(trace_t[it][1],2))/(float)lt;

      for(it=0;it<lt;it++)
        swq2.write((char*)&u2[it],sizeof(u2[it]));     

       for(it=0;it<lt;it++)
        {
          u3[it]=u1[it];

         if(u2[it]<thre)
           u3[it]=0.0;
        }

      for(it=0;it<lt;it++)
        swq3.write((char*)&u3[it],sizeof(u3[it]));     

     cout<<is+1<<" trace done..."<<endl;
 
   }

  cout<<"All Done!"<<endl;

  return 0;  
}















































