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

using namespace std;

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn11[256],fn22[256],fn33[256],fn44[256],fn5[256],fn8[256];
   int nshot,nr,lt,shot,trace;
   float ds,dr,dt,fmax,fmin,st_beg,st_end;
   int ltt;
   int ifmax,ifmin;
   int is,ir,it,ix,ircrg,ix1,mcg_count,tmp;
 
   ifstream swq;
   swq.open("convolution_prediction_obc.par");
   if(!swq)
      {
         cout<<"cannot open convolution_prediction_obc.par"<<endl;
         abort();
      }

   swq>>fn1>>fn11>>fn2>>fn22>>fn3>>fn33>>fn4>>fn44>>fn5>>fn8>>nshot>>nr>>ds>>dr>>lt>>dt>>fmax>>fmin>>shot>>trace>>st_beg>>st_end;
   swq.close();

   cout<<"Fna of CSG Source-side Up-wavefield real part is===="<<fn1<<endl;
   cout<<"Fna of CSG Source-side Up-wavefield imaginary part is===="<<fn11<<endl;
   cout<<"Fna of CSG Source-side Down-wavefield real part is===="<<fn2<<endl;
   cout<<"Fna of CSG Source-side Down-wavefield imaginary part is===="<<fn22<<endl;
   cout<<"Fna of CRG Receiver-side Up-wavefield real part is===="<<fn3<<endl;
   cout<<"Fna of CRG Receiver-side Up-wavefield imaginary part is===="<<fn33<<endl;
   cout<<"Fna of CRG Receiver-side Down-wavefield real part is===="<<fn4<<endl;
   cout<<"Fna of CRG Receiver-side Down-wavefield imaginary part is===="<<fn44<<endl;
   cout<<"Fna of Predicted Multiple is===="<<fn5<<endl;
   cout<<"Fna of Predicted MCG is===="<<fn8<<endl;
   cout<<"No. of Shots is===="<<nshot<<endl; 
   cout<<"No. of Traces per Shot is===="<<nr<<endl;
   cout<<"Source Interval is===="<<ds<<"m"<<endl; 
   cout<<"Receiver Interval is===="<<dr<<"m"<<endl; 
   cout<<"lt and dt are===="<<lt<<" , "<<dt<<"ms"<<endl;
   cout<<"Maximum and Minimum Frequency is===="<<fmax<<"Hz , "<<fmin<<"Hz"<<endl;
   cout<<"The Output MCG of Shot No. is===="<<shot<<endl;
   cout<<"The Output MCG of Trace No. is===="<<trace<<endl;
   cout<<"The qualitive ratio is===="<<st_beg<<" , "<<st_end<<endl; 
   
   ltt=2*lt;

   float *omega;
   omega=alloc1float(ltt);
  
//calculate the omega and kx
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   ifmin=int(fmin*dt*ltt/1000);
   ifmax=int(fmax*dt*ltt/1000)+1;

   cout<<"Frequency slices to be calculated is==== "<<ifmax-ifmin+1<<endl;

   complex<float> **su;
   su=alloc2complex(ltt,nr);

   complex<float> **sd;
   sd=alloc2complex(ltt,nr);

   complex<float> **ru;
   ru=alloc2complex(ltt,nr);

   complex<float> **rd;
   rd=alloc2complex(ltt,nr);

 
   complex<float> **mcg;
   mcg=alloc2complex(ltt,nr);
   
   float **mcgt;
   mcgt=alloc2float(lt,nr);

   float **srmt;
   srmt=alloc2float(lt,nr);

   float mem;
   mem=0.0;
   
   mem+=(ltt*nr*2*6+ltt+lt)*4;
   mem=mem/1024/1024/4;
   cout<<"Memory Needed to be allocated is===="<<mem<<"MB"<<endl; 
 
//   cout<<"=========="<<endl;

   fftwf_complex *in1,*out1;
   fftwf_plan p1;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_ESTIMATE);

//   cout<<"=========="<<endl;

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
      {
         cout<<"cannot open "<<fn1<<endl;
         abort();
      } 
   ifstream swq11;
   swq11.open(fn11,ios::binary);
   if(!swq11)
      {
         cout<<"cannot open "<<fn11<<endl;
         abort();
      }
   
   ifstream swq2;
   swq2.open(fn2,ios::binary);
   if(!swq2)
      {
         cout<<"cannot open "<<fn2<<endl;
         abort();
      } 
   ifstream swq22;
   swq22.open(fn22,ios::binary);
   if(!swq22)
      {
         cout<<"cannot open "<<fn22<<endl;
         abort();
      } 

   ifstream swq3;
   swq3.open(fn3,ios::binary);
   if(!swq3)
      {
         cout<<"cannot open "<<fn3<<endl;
         abort();
      } 
   ifstream swq33;
   swq33.open(fn33,ios::binary);
   if(!swq33)
      {
         cout<<"cannot open "<<fn33<<endl;
         abort();
      } 

   ifstream swq4;
   swq4.open(fn4,ios::binary);
   if(!swq4)
      {
         cout<<"cannot open "<<fn4<<endl;
         abort();
      } 
   ifstream swq44;
   swq44.open(fn44,ios::binary);
   if(!swq44)
      {
         cout<<"cannot open "<<fn44<<endl;
         abort();
      }
   
   ofstream swq5;
   swq5.open(fn5,ios::binary);
   if(!swq5)
      {
         cout<<"cannot open "<<fn5<<endl;
         abort();
      }

   ofstream swq8;
   swq8.open(fn8,ios::binary);
   if(!swq8)
      {
         cout<<"cannot open "<<fn8<<endl;
         abort();
      }
   cout<<"====SRM Prediction in CSG domain Starts===="<<endl;
   
 
   for(is=0;is<nshot;is++) 
     {
        cout<<is+1<<" shot prediction starts..."<<endl;

        for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++)
            srmt[ir][it]=0.0;

        swq1.seekg(0,ios::beg);
        swq11.seekg(0,ios::beg);
        swq2.seekg(0,ios::beg);
        swq22.seekg(0,ios::beg);
        for(ix=0;ix<is;ix++)  
           {
             swq1.seekg(nr*ltt*4,ios::cur); 
             swq11.seekg(nr*ltt*4,ios::cur); 
             swq2.seekg(nr*ltt*4,ios::cur); 
             swq22.seekg(nr*ltt*4,ios::cur); 
           }
        for(ix=0;ix<nr;ix++)
          for(it=0;it<ltt;it++)
           {
              swq1.read((char*)&(su[ix][it].real()),sizeof(su[ix][it].real()));
              swq11.read((char*)&(su[ix][it].imag()),sizeof(su[ix][it].imag()));
              swq2.read((char*)&(sd[ix][it].real()),sizeof(sd[ix][it].real()));
              swq22.read((char*)&(sd[ix][it].imag()),sizeof(sd[ix][it].imag()));
           } 

         for(ir=0;ir<nr;ir++)
           {
              swq3.seekg(0,ios::beg);
              swq33.seekg(0,ios::beg);
              swq4.seekg(0,ios::beg);
              swq44.seekg(0,ios::beg);

              for(ix=0;ix<ir;ix++)
               {
                 swq3.seekg(nshot*ltt*4,ios::cur);
                 swq33.seekg(nshot*ltt*4,ios::cur);
                 swq4.seekg(nshot*ltt*4,ios::cur);
                 swq44.seekg(nshot*ltt*4,ios::cur);
               }

             for(ix=0;ix<nr;ix++)
                for(it=0;it<ltt;it++)                
                  mcg[ix][it]=(0.0,0.0);  

             for(ix=0;ix<nr;ix++)
               for(it=0;it<ltt;it++)
                 {
                   swq3.read((char*)&(ru[ix][it].real()),sizeof(ru[ix][it].real()));
                   swq33.read((char*)&(ru[ix][it].imag()),sizeof(ru[ix][it].imag()));
                   swq4.read((char*)&(rd[ix][it].real()),sizeof(rd[ix][it].real()));
                   swq44.read((char*)&(rd[ix][it].imag()),sizeof(rd[ix][it].imag()));
                 }

//start to convolve the CSG and CRG
            if(is<ir) 
             {  
              for(ix=is;ix<ir;ix++)
                {
                    for(it=ifmin;it<ifmax;it++)
                     {
                       mcg[ix][it]=su[ix][it]*rd[ix][it]+sd[ix][it]*ru[ix][it];
//                       mcg[ix][it]=su[ix][it]*rd[ix][it];
//                       mcg[ix][it]=sd[ix][it]*rd[ix][it];

                       mcg[ix][it].real()/=(ir-is);
                       mcg[ix][it].imag()/=(ir-is);
                     }
                }
             }
            else if(is==ir) 
             {  
                    for(it=ifmin;it<ifmax;it++)
                     {
                       mcg[is][it]=su[is][it]*rd[is][it]+sd[is][it]*ru[is][it];
//                       mcg[is][it]=su[is][it]*rd[is][it];
//                       mcg[is][it]=sd[is][it]*rd[is][it];
                     }
             }
            else
             {  
              for(ix=ir;ix<is;ix++)
                {
                    for(it=ifmin;it<ifmax;it++)
                     {
                       mcg[ix][it]=su[ix][it]*rd[ix][it]+sd[ix][it]*ru[ix][it];
//                       mcg[ix][it]=su[ix][it]*rd[ix][it];
//                       mcg[ix][it]=sd[ix][it]*rd[ix][it];

                       mcg[ix][it].real()/=(is-ir);
                       mcg[ix][it].imag()/=(is-ir);
                     }
                }
             }

              for(ix=0;ix<nr;ix++)
                 for(it=ltt/2+1;it<ltt;it++)
                    {
                      mcg[ix][it].real()=mcg[ix][ltt-it].real();
                      mcg[ix][it].imag()=-mcg[ix][ltt-it].imag();
                    }
            

              for(ix=0;ix<nr;ix++)
                {
                   for(it=0;it<ltt;it++)
                     {
                       in1[it][0]=mcg[ix][it].real();
                       in1[it][1]=mcg[ix][it].imag();
                     }
                    fftwf_execute(p1);

                    for(it=0;it<lt;it++)
                       mcgt[ix][it]=out1[it][0];
                }

                 for(ix=0;ix<nr;ix++)
                   for(it=0;it<lt;it++)
                    swq8.write((char*)&(mcgt[ix][it]),sizeof(mcgt[ix][it]));
/*
              if(is==shot) 
                {
                 for(ix=0;ix<nr;ix++)
                   for(it=0;it<lt;it++)
                    swq8.write((char*)&(mcgt[ix][it]),sizeof(mcgt[ix][it]));
                }
*/
              for(it=0;it<lt;it++)
                { 
                  for(ix=0;ix<nr;ix++)
                     srmt[ir][it]+=mcgt[ix][it];
                }

              for(it=0;it<lt;it++)
                 swq5.write((char*)&(srmt[ir][it]),sizeof(srmt[ir][it]));

           }

        cout<<is+1<<" shot prediction done..."<<endl;
 
     }

   cout<<"ALL DONE!"<<endl; 
 
      return 0;
}























































































































