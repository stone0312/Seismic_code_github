#include "iostream.h"
#include "fstream.h"
#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "stdlib.h"
#include "fftw3.h"
#define pai 3.14159265

using namespace std;

int main()
{
   int lt;
   float dt,rz,v;

   lt=7001;
   dt=0.5;
   rz=100.0;
   v=1500.0;
 
   int ltt;
   ltt=1*lt;
 
   float lamda;
   lamda=100;
  
   int it,tmp;
 
   float degh_amp_max;
 
   float *u;
   u=alloc1float(ltt);

   complex<float> *uf;
   uf=alloc1complex(ltt);

   complex<float> *ps;
   ps=alloc1complex(ltt);

   complex<float> *dg;
   dg=alloc1complex(ltt);

   complex<float> *udgf;
   udgf=alloc1complex(ltt);
  
   float *udg;
   udg=alloc1float(ltt);

   float *omega;
   omega=alloc1float(ltt);
   for(it=0;it<ltt/2+1;it++)
     omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
     omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

  float *degh_amp;
  degh_amp=alloc1float(ltt);

  float *amp_orig;
  amp_orig=alloc1float(ltt);

  float *filt;
  filt=alloc1float(ltt);

  fftwf_complex *in1,*out1,*in2,*out2;
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

  fftwf_plan p1,p2;

  p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
  p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);



  ifstream swq1;
  swq1.open("./1st_trace_primary_ghost_100m_7001.dat",ios::binary);
  if(!swq1)
     {
         cout<<"cannot open 1st_trace_primary_ghost_100m_7001.dat"<<endl;
         abort();
     }

  ofstream swq2;
  swq2.open("./orig_amp_1st_trace_100m_7001.dat",ios::binary);

  ofstream swq3;
  swq3.open("./deghost_operator_amp_1st_trace_100m_7001.dat",ios::binary);

  ofstream swq4;
  swq4.open("./deghost_amp_1st_trace_100m_lamda100_7001.dat",ios::binary);
 
  ofstream swq5;
  swq5.open("./1st_trace_100m_deghosting_lamda100_7001.dat",ios::binary);

  for(it=0;it<ltt;it++)
    u[it]=0.0;

  for(it=0;it<ltt;it++)
    swq1.read((char*)&u[it],sizeof(u[it]));
  swq1.close();      

  for(it=0;it<ltt;it++)
    {
       in1[it][0]=u[it];
       in1[it][1]=0.0;
    }

   fftwf_execute(p1);

   for(it=0;it<ltt;it++)
     {
        uf[it].real()=out1[it][0]/ltt;
        uf[it].imag()=out1[it][1]/ltt;
     }

   for(it=0;it<ltt;it++)
     {
       amp_orig[it]=sqrt(pow(uf[it].real(),2)+pow(uf[it].imag(),2));    
       swq2.write((char*)&amp_orig[it],sizeof(amp_orig[it]));
     }
   swq2.close();

   for(it=0;it<ltt;it++)
     {
         ps[it].real()=0.0;
         ps[it].imag()=2*rz*omega[it]/v;
 
         ps[it]=exp(ps[it]);

         dg[it].real()=1-ps[it].real();
         dg[it].imag()=ps[it].imag();         

         degh_amp[it]=sqrt(pow(dg[it].real(),2)+pow(dg[it].imag(),2)); 
         swq3.write((char*)&degh_amp[it],sizeof(degh_amp[it]));
     }
    swq3.close();
/*
   degh_amp_max=0.0;
   for(it=0;it<ltt;it++)
     {
        if(degh_amp[it]>degh_amp_max)
           degh_amp_max=degh_amp[it];
        else
           degh_amp_max=degh_amp_max;           
     }
   cout<<"maximum amp of deghosting operator is===="<<degh_amp_max<<endl;

    tmp=0; 
//    for(it=0;it<lt;it++)
//      cout<<it<<"  "<<degh_amp[it]/degh_amp_max<<endl;
//        cout<<it<<"  "<<degh_amp[it]<<endl;

    for(it=0;it<ltt;it++)
     {  
        if(degh_amp[it]/degh_amp_max>0.1)
          udgf[it]=uf[it]/dg[it];
        else
          {
             tmp+=1;

             dg[it].real()+=lamda;
             udgf[it]=uf[it]/dg[it];
          }


         filt[it]=sqrt(pow(udgf[it].real(),2)+pow(udgf[it].imag(),2));
         swq4.write((char*)&filt[it],sizeof(filt[it]));
     } 
    
    cout<<tmp<<endl; 
    swq4.close();
*/
    udgf[0]=(0.0,0.0);

    for(it=1;it<ltt;it++)
      udgf[it]=uf[it]/dg[it];


    for(it=0;it<ltt;it++)
     {
       filt[it]=sqrt(pow(udgf[it].real(),2)+pow(udgf[it].imag(),2));
       swq4.write((char*)&filt[it],sizeof(filt[it]));
     }
    swq4.close();
  

    for(it=0;it<ltt;it++)
    {
       in2[it][0]=udgf[it].real();
       in2[it][1]=udgf[it].imag();
    }

   fftwf_execute(p2);
   
   for(it=0;it<ltt;it++)
     udg[it]=out2[it][0];
   
   for(it=0;it<ltt;it++)
     swq5.write((char*)&udg[it],sizeof(udg[it]));
   swq5.close();

   return 0;

}



































