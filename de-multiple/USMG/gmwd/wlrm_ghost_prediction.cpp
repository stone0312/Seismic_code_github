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

int deghosting_plus_extrapolation(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfm, complex <float> ** ufm, float **um, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, int ifmin, int ifmax, float wd)
{
    float kz;
    complex <float> ps1;

    int ir,is, it;

//zeroing the arrays
    for(ir=0;ir<nrb;ir++) 
       for(it=0;it<ltt;it++)
         {
            u[ir][it]=0.0; 
            uf[ir][it]=(0.0,0.0);
            ukf[ir][it]=(0.0,0.0);
            ukfm[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
            ufm[ir][it]=(0.0,0.0);

     for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
            um[ir][it]=0.0;

    for(ir=boul;ir<nr+boul;ir++)
      for(it=0;it<lt;it++)
           u[ir][it]=uraw[ir-boul][it];
 
    fftwf_complex *in1,*out1;
    fftwf_plan p1;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    p2=fftwf_plan_dft_1d(nrb,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

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
       for(it=0;it<ltt;it++)
          uf[ir][it]=(0.0,0.0); 

    for(ir=0;ir<nrb;ir++)
       {
         for(it=0;it<ltt;it++)
           {
              in1[it][0]=u[ir][it];
              in1[it][1]=0.0;
           }

          fftwf_execute(p1); 

          for(it=0;it<ltt;it++)
            {
              uf[ir][it].real()=out1[it][0]/ltt;
              uf[ir][it].imag()=out1[it][1]/ltt;
            }
       }  

    for(it=0;it<ltt;it++)
       {
           for(ir=0;ir<nrb;ir++)
              {
                 in2[ir][0]=uf[ir][it].real();
                 in2[ir][1]=uf[ir][it].imag();
              } 
           fftwf_execute(p2); 
          

           for(ir=0;ir<nrb;ir++)
             {
               ukf[ir][it].real()=out2[ir][0]/nrb;
               ukf[ir][it].imag()=out2[ir][1]/nrb;
             }
       }

     cout<<"Prodiction Begins..."<<endl;

     for(ir=0;ir<nrb;ir++)
        {
           for(it=0;it<ltt/2+1;it++)
              {
                 if(it<ifmin)
                   {
                      ukfm[ir][it].real()=0.0;
                      ukfm[ir][it].imag()=0.0;
                   }
                 else if(it>ifmax&&it<ltt/2+1)
                    {
                      ukfm[ir][it].real()=0.0;
                      ukfm[ir][it].imag()=0.0;
                   }
                 else
                   {
                      if(pow(omega[it]/v,2)-pow(kx[ir],2)<0.0)
                        {
                           ukfm[ir][it].real()=0.0;
                           ukfm[ir][it].imag()=0.0;
                        }
                      else if(pow(omega[it]/v,2)-pow(kx[ir],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                        {
                            kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));
                             ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                             ps1.imag()=-2*wd*kz;
                             ps1=exp(ps1);

                             ukfm[ir][it]=ukf[ir][it]*ps1;   //upward field

                        }
                       else
                        {
                           ukfm[ir][it].real()=0.0;
                           ukfm[ir][it].imag()=0.0;
                        } 
                   }
//                 if(it%1000==0)
//                    cout<<it<<" Frequency Slices Deghosting Done..."<<endl;
                }

            for(it=ltt/2+1;it<ltt;it++)
              {
                 ukfm[ir][it].real()=ukfm[ir][ltt-it].real();
                 ukfm[ir][it].imag()=-ukfm[ir][ltt-it].imag();
              }
        }

//ifft

     cout<<"IFFT Begins..."<<endl;
 
       for(it=0;it<ltt;it++)
         {   
           for(ir=0;ir<nrb;ir++)
              {   
                 in3[ir][0]=ukf[ir][it].real();
                 in3[ir][1]=ukf[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               ufm[ir][it].real()=out3[ir+boul][0];
               ufm[ir][it].imag()=out3[ir+boul][1];
             }   
         }   


        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufm[ir][it].real();
                 in4[it][1]=ufm[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             um[ir][it]=-out4[it][0];
         }


    return 0;

}

int main()
{
   char fn1[256],fn2[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb, n;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda,wd;
 
   ifstream swq;
   swq.open("wlrm_prediction.par");
   swq>>fn1>>fn2>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>lamda>>wd>>n;
   swq.close();

   cout<<"Fna of input shot gather is===="<<fn1<<endl; 
   cout<<"Fna of predicted multiple is===="<<fn2<<endl; 

   cout<<"No. of shots is===="<<ns<<endl;
   cout<<"No. of traces per shot is===="<<nr<<endl;
   cout<<"lt and dt is===="<<lt<<" , "<<dt<<"ms"<<endl;
   cout<<"Source and Receiver Interval are===="<<ds<<"m , "<<dr<<"m"<<endl;
   cout<<"Source and Receiver Depth are===="<<sdep<<"m , "<<rdep<<"m"<<endl;
   cout<<"Water Velocity is===="<<v<<endl;
   cout<<"Maximum Theta is===="<<theta<<endl;
   cout<<"Left and Right Boundary are===="<<boul<<" , "<<bour<<endl;  
   cout<<"Maximum and Minimum Frequency are===="<<fmax<<"Hz , "<<fmin<<"Hz"<<endl;
   cout<<"Regulation Parameter is===="<<lamda<<endl;

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

   float **uraw;
   uraw=alloc2float(lt,nr);

   float **u;
   u=alloc2float(ltt,nrb);

   complex<float> **uf;
   uf=alloc2complex(ltt,nrb); 

   complex<float> **ukf;
   ukf=alloc2complex(ltt,nrb);

   complex<float> **ukfm;
   ukfm=alloc2complex(ltt,nrb);

   complex<float> **ufm;
   ufm=alloc2complex(ltt,nr);

   float **um;
   um=alloc2float(lt,nr);

  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl;

  ofstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;
 
  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
            swq1.read((char*)&uraw[ir][it],sizeof(uraw[ir][it]));

       deghosting_plus_extrapolation(uraw, u, uf, ukf, ukfm, ufm, um, omega, kx,  sdep, rdep, ds, dr, nr, lt, ltt, boul, bour, nrb, v, theta, ifmin, ifmax, wd);

       for(ir=n;ir<nr;ir++)
        for(it=0;it<lt;it++)
           swq2.write((char*)&um[ir][it],sizeof(um[ir][it]));

        cout<<is+1<<" shot done!"<<endl;
    }
  
        return 0; 

}
































