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

int green_function( complex <float> *sfk, complex <float> **gfk,  complex <float> **gf, float **g, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, int ifmin, int ifmax, float wd)
{
    float kz;
    complex <float> ps1;
    complex <float> ps2;
    complex <float> ps3;
    complex <float> ps4;
    complex <float> ps5;
    complex <float> ps6;
    complex <float> ps7;

    int ir,is,it;

//zeroing the arrays
    for(ir=0;ir<nrb;ir++) 
       for(it=0;it<ltt;it++)
         {
            gfk[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
          {
            gf[ir][it]=(0.0,0.0);
            g[ir][it]=0.0;
          }

    fftwf_complex *in3,*out3;
    fftwf_plan p3;
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    p3=fftwf_plan_dft_1d(nrb,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);

    fftwf_complex *in4,*out4;
    fftwf_plan p4;
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p4=fftwf_plan_dft_1d(ltt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

     for(ir=0;ir<nrb;ir++)
        {
           for(it=0;it<ltt/2+1;it++)
              {
                 if(it<ifmin)
                   {
                      gfk[ir][it].real()=0.0;
                      gfk[ir][it].imag()=0.0;
                   }
                 else if(it>ifmax&&it<ltt/2+1)
                    {
                      gfk[ir][it].real()=0.0;
                      gfk[ir][it].imag()=0.0;
                   }
                 else
                   {
                      if(pow(omega[it]/v,2)-pow(kx[ir],2)<0.0)
                        {
                           gfk[ir][it].real()=0.0;
                           gfk[ir][it].imag()=0.0;
                        }
                      else if(pow(omega[it]/v,2)-pow(kx[ir],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                        {
                             kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));

                             ps1.real()=0.0;               
                             ps1.imag()=-2*wd*kz;
                             ps1=exp(ps1);

                             gfk[ir][it]=sfk[ir]*ps1; 
                        }
                       else
                        {
                           gfk[ir][it].real()=0.0;
                           gfk[ir][it].imag()=0.0;
                        } 
                   }

                }
        }


//ifft

     cout<<"IFFT Begins..."<<endl;
 
       for(it=0;it<ltt/2+1;it++)
         {   
           for(ir=0;ir<nrb;ir++)
              {   
                 in3[ir][0]=gfk[ir][it].real();
                 in3[ir][1]=gfk[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               gf[ir][it].real()=out3[ir+boul][0];
               gf[ir][it].imag()=out3[ir+boul][1];
             }   
         }   

     for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
          {
             gf[ir][it].real()=gf[ir][ltt-it].real();
             gf[ir][it].imag()=-gf[ir][ltt-it].imag();
          }

        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=gf[ir][it].real();
                 in4[it][1]=gf[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             g[ir][it]=out4[it][0];
         }


    return 0;

}

int main()
{
   char fn1[256],fn2[256],fn3[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb, n;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda,wd;
 
   ifstream swq;
   swq.open("green_function.par");
   swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>wd;
   swq.close();

   cout<<"Fna of Green's Function in x-f domain Real Parts is===="<<fn1<<endl; 
   cout<<"Fna of Green's Function in x-f domain Imaginary Parts is===="<<fn2<<endl; 
   cout<<"Fna of Green's Function in x-t domain is===="<<fn3<<endl; 
   cout<<"No. of shots is===="<<ns<<endl;
   cout<<"No. of traces per shot is===="<<nr<<endl;
   cout<<"lt and dt is===="<<lt<<" , "<<dt<<"ms"<<endl;
   cout<<"Source and Receiver Interval are===="<<ds<<"m , "<<dr<<"m"<<endl;
   cout<<"Source and Receiver Depth are===="<<sdep<<"m , "<<rdep<<"m"<<endl;
   cout<<"Water Velocity is===="<<v<<endl;
   cout<<"Maximum Theta is===="<<theta<<endl;
   cout<<"Left and Right Boundary are===="<<boul<<" , "<<bour<<endl;  
   cout<<"Maximum and Minimum Frequency are===="<<fmax<<"Hz , "<<fmin<<"Hz"<<endl;
   cout<<"Water Bottom Depth is===="<<wd<<endl;

   ltt=2*lt;

   nrb=boul+nr+bour;

   int it,is,ir;

   float *omega;
   omega=alloc1float(ltt);
   float *kx;
   kx=alloc1float(nrb);

   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   ifmin=int(fmin*dt*ltt/1000);
   ifmax=int(fmax*dt*ltt/1000)+1;

   cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

   for(is=0;is<nrb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dr*nrb);
   for(is=nrb/2+1;is<nrb;is++)
      kx[is]=-2*pai*1.0/float(2*dr)+2*pai*float(is-nrb/2)/float(dr*nrb);
   
   complex<float> *sf;
   sf=alloc1complex(nrb); 

   complex<float> *sfk;
   sfk=alloc1complex(nrb);

   complex<float> **gfk;
   gfk=alloc2complex(ltt,nrb);

   complex<float> **gf;
   gf=alloc2complex(ltt,nr);

   float **g;
   g=alloc2float(ltt,nr);

   fftwf_complex *in2,*out2;
   fftwf_plan p2;
   in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
   out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
   p2=fftwf_plan_dft_1d(nrb,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

  ofstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl;

  ofstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;

  ofstream swq3;
  swq3.open(fn3,ios::binary);
  if(!swq3)
       cout<<"cannot open "<<fn3<<endl;
 
  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nrb;ir++)
          sf[ir]=(0.0,0.0);

       sf[is+boul].real()=1.0;
       sf[is+boul].imag()=0.0;
        
       for(ir=0;ir<nrb;ir++)
          {
            in2[ir][0]=sf[ir].real();
            in2[ir][1]=sf[ir].imag();
          }
       
        fftwf_execute(p2);

        for(ir=0;ir<nrb;ir++)
         {
           sfk[ir].real()=out2[ir][0]/nrb;
           sfk[ir].imag()=out2[ir][1]/nrb;
         }

       green_function(sfk, gfk, gf, g, omega, kx, sdep, rdep, ds, dr, nr, lt, ltt, boul, bour, nrb, v, theta, ifmin, ifmax, wd);

       for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
          { 
           swq1.write((char*)&gf[ir][it].real(),sizeof(gf[ir][it].real()));
           swq2.write((char*)&gf[ir][it].imag(),sizeof(gf[ir][it].imag()));
          }
  
       for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
           swq3.write((char*)&g[ir][it],sizeof(g[ir][it]));

        cout<<is+1<<" shot done!"<<endl;
    }
  
        return 0; 

}
































