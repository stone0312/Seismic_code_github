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
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256];
  int ns,nr,np,lt,ishot,itrace;
  float dt,fmin,fmax,min,max;

  ifstream swq;
  swq.open("usmg_tp.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>ns>>nr>>np>>lt>>dt>>fmin>>fmax>>min>>max>>ishot>>itrace;
  swq.close();

  cout<<"Fna of CSG_UP Real is===="<<fn1<<endl; 
  cout<<"Fna of CSG_UP Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG_DOWN Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG_DOWN Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<" , "<<dt<<endl;
  cout<<"minimum and maximum frequency to be calculated are===="<<fmin<<" , "<<fmax<<endl;

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1,ip; 
  ltt=2*lt;
 
  float *omega;
  omega=alloc1float(ltt);
  for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
  for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));
 
  ifmin=int(fmin*dt*ltt/1000);
  ifmax=int(fmax*dt*ltt/1000)+1;

  cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;
 
  complex<float> **usf;
  usf=alloc2complex(ltt,nr*np);
  complex<float> **urf;
  urf=alloc2complex(ltt,ns*np);
  complex<float> **urf1;
  urf1=alloc2complex(ltt,ns*np);

  complex<float> **mcgpf;
  mcgpf=alloc2complex(ltt,nr*np);
  float **mcgpt;
  mcgpt=alloc2float(lt,nr*np);
  float **mcgt;
  mcgt=alloc2float(lt,nr);

  float **mt;
  mt=alloc2float(lt,nr);

  float **urpt1;
  urpt1=alloc2float(lt,nr*np);

  fftwf_complex *in1,*out1;
  fftwf_plan p1;
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);
 
  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl; 

  ifstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;

  ifstream swq3;
  swq3.open(fn3,ios::binary);
  if(!swq3)
       cout<<"cannot open "<<fn3<<endl;

  ifstream swq4;
  swq4.open(fn4,ios::binary);
  if(!swq4)
       cout<<"cannot open "<<fn4<<endl;

  ofstream swq5;
  swq5.open(fn5,ios::binary);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;

  ofstream swq6;
  swq6.open(fn6,ios::binary);
  if(!swq6)
       cout<<"cannot open "<<fn6<<endl;

  ofstream swq7;
  swq7.open(fn7,ios::binary);
  if(!swq7)
       cout<<"cannot open "<<fn7<<endl;
/*
  ofstream swqtmp;
  swqtmp.open("./check_tp_read_new.dat",ios::binary);
  if(!swqtmp)
       cout<<"cannot open check_tp_read.dat"<<endl;
*/

//  for(is=0;is<ns;is++)
  for(is=10;is<11;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
           mt[ir][it]=0.0;

      for(ir=0;ir<nr*np;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real())); 
           swq2.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
         }
       
           swq3.seekg(0,ios::beg);
           swq4.seekg(0,ios::beg);

           for(ir=0;ir<nr;ir++)
              {
                for(ir1=0;ir1<nr;ir1++)
                 for(it=0;it<lt;it++)
                     mcgt[ir1][it]=0.0;

                for(ir1=0;ir1<nr*np;ir1++)
                  for(it=0;it<ltt;it++)
                    {
                       mcgpf[ir1][it].real()=0.0;
                       mcgpf[ir1][it].imag()=0.0;
                    }

                for(ir1=0;ir1<nr*np;ir1++)
                 for(it1=0;it1<ltt;it1++)
                  {
                    swq3.read((char*)&urf[ir1][it1].real(),sizeof(urf[ir1][it1].real()));
                    swq4.read((char*)&urf[ir1][it1].imag(),sizeof(urf[ir1][it1].imag()));
                  }

                for(it1=0;it1<ltt;it1++)
                  {
                     for(ir1=0;ir1<nr;ir1++)
                       {
                          for(ip=ir1*np;ip<(ir1+1)*np;ip++)
                            urf1[ip][it1]=urf[2*ir1*np+np-1-ip][it1];
                       }
                  }

                for(it=ifmin;it<ifmax;it++)
                 {
                   for(ir1=0;ir1<nr*np;ir1++)
                     mcgpf[ir1][it]=usf[ir1][it]*urf1[ir1][it]; 
                 } 

                for(ir1=0;ir1<nr*np;ir1++)
                 for(it=ltt/2+1;it<ltt;it++)
                   {   
                     mcgpf[ir1][it].real()=mcgpf[ir1][ltt-it].real();   
                     mcgpf[ir1][it].imag()=-mcgpf[ir1][ltt-it].imag();   
                   }
 
                 for(ir1=0;ir1<nr*np;ir1++)
                  {
      	            for(it=0;it<ltt;it++)
                      {
                        in1[it][0]=mcgpf[ir1][it].real();
                        in1[it][1]=mcgpf[ir1][it].imag();
                      }

                    fftwf_execute(p1);

                    for(it=0;it<lt;it++)
                      mcgpt[ir1][it]=out1[it][0];
                  }

                 for(it=0;it<lt;it++)
                  {
                   for(ir1=0;ir1<nr;ir1++)
                     {
                       for(ip=ir1*np;ip<(ir1+1)*np;ip++)
                         mcgt[ir1][it]+=mcgpt[ip][it];  
                     }
                  }

                if(is==ishot&&ir==itrace)
                  {
                     for(ir1=0;ir1<nr*np;ir1++)
                       for(it=0;it<lt;it++)
                          swq6.write((char*)&mcgpt[ir1][it],sizeof(mcgpt[ir1][it]));

                     for(ir1=0;ir1<nr;ir1++)
                       for(it=0;it<lt;it++)
                          swq7.write((char*)&mcgt[ir1][it],sizeof(mcgt[ir1][it]));
                  }

                 for(it=0;it<lt;it++)
                   {
                      for(ir1=0;ir1<nr;ir1++)
                         mt[ir][it]+=mcgt[ir1][it];
                   }

                 for(it=0;it<lt;it++)
        	   swq5.write((char*)&mt[ir][it],sizeof(mt[ir][it]));

                  cout<<"   "<<is+1<<" Shot, "<<ir+1<<" Trace Forward Prediction Done ..."<<endl;


              }             
/*
      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq5.write((char*)&mt[ir][it],sizeof(mt[ir][it]));
*/

      cout<<is+1<<" Shot Forward Prediction Done!"<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











