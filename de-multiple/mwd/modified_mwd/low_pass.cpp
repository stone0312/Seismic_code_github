#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

int main()
{
  int ns,nr,nsr,lt;
  float dt,fbeg,fend;
  ns=300;
  nr=479;
  nsr=ns*nr;
  lt=3000;
  dt=2.0;  //ms
  fbeg=60.0;
  fend=100.0;

  int ix,ir,it,ifbeg,ifend,itmp;

  ifbeg=int(fbeg*dt*lt/1000);
  ifend=int(fend*dt*lt/1000);
/*  
  cout<<ifbeg<<"  "<<ifend<<endl;
  return 0 ;
*/
  float *hann;
  hann=alloc1float(ifend-ifbeg+1);
  for(it=0;it<ifend-ifbeg+1;it++)
    hann[it]=cos(0.5*pai*it/(ifend-ifbeg));     
/* 
  for(it=0;it<ifend-ifbeg+1;it++)
    cout<<it<<"  "<<hann[it]<<endl;
  return 0;
*/

  float *u1;
  u1=alloc1float(lt);

  complex<float> *u2;
  u2=alloc1complex(lt);

  complex<float> *u3;
  u3=alloc1complex(lt);

  float *u;
  u=alloc1float(lt);

  ifstream swq2;
  swq2.open("/data2/swq/100th_shot_wlrm_zj_dz0125m.dat",ios::binary);
  if(!swq2)
    {
       cout<<"cannot open file"<<endl;
       abort();
    } 
    
  ofstream swq3;
  swq3.open("/data2/swq/lowpass_100th_shot_wlrm_zj_dz0125m.dat",ios::binary);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in3,*out3;
    fftwf_plan p3;
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    p3=fftwf_plan_dft_1d(lt,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
   
        for(ir=0;ir<nr;ir++)
         {
          for(it=0;it<lt;it++)
           swq2.read((char*)&u1[it],sizeof(u1[it]));
         
         if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(it=0;it<lt;it++)
                 {
                    in2[it][0]=u1[it];
                    in2[it][1]=0.0;
                 }
           }

        fftwf_execute(p2);

        for(it=0;it<lt;it++)
          {
             u2[it].real()=out2[it][0];
             u2[it].imag()=out2[it][1];
          }

        for(it=0;it<ifbeg;it++)
          u3[it]=u2[it];
        for(it=ifbeg;it<ifend+1;it++)
          {
            itmp=it-ifbeg;
//            u3[it]=u2[it]*hann[itmp];
            u3[it].real()=0.0;
            u3[it].imag()=0.0;
          }
        for(it=ifend+1;it<lt/2+1;it++)
          {
            u3[it].real()=0.0;
            u3[it].imag()=0.0;
          } 
        for(it=lt/2+1;it<lt;it++)
          {
             itmp=lt-it;
             u3[it].real()=u3[itmp].real(); 
             u3[it].imag()=-u3[itmp].imag();
          }


       if((in3==NULL)||(out3==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(it=0;it<lt;it++)
                 {
                    in3[it][0]=u3[it].real();
                    in3[it][1]=u3[it].imag();
                 }
           }

        fftwf_execute(p3);

        for(it=0;it<lt;it++)
          u[it]=out3[it][0]/lt;

         for(it=0;it<lt;it++)
            swq3.write((char*)&u[it],sizeof(u[it]));
        cout<<ir<<" trace done!"<<endl;
         }

     swq2.close();
     swq3.close();

     return 0;

}






































