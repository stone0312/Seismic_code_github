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

int source_deghosting_plus_wlrm(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu, complex <float> ** ukfd,  complex <float> **ufu, complex <float> **ufd,  float **uu, float **ud,  complex <float> **ukfuw, complex <float> ** ukfuwg,  complex <float> **ufuw, complex <float> **ufuwg, float **uuwlrm, float **uuwlrmg,  float *omega, float *kx, float sdep,float wd, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, float ifmin, float ifmax, float lamda)
{
   float kz;
   complex <float> ps1;
   complex <float> ps2;
   complex <float> ps3;
   complex <float> ps4;
   complex <float> ps5;
   complex <float> ps6;
  
   int ir,is,it; 

   float amp;

    for(ir=0;ir<nrb;ir++)
       for(it=0;it<ltt;it++)
         {
            u[ir][it]=0.0;
            uf[ir][it]=(0.0,0.0);
            ukf[ir][it]=(0.0,0.0);
            ukfu[ir][it]=(0.0,0.0);
            ukfd[ir][it]=(0.0,0.0);
            ukfuw[ir][it]=(0.0,0.0);
            ukfuwg[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
         {
            ufu[ir][it]=(0.0,0.0);
            ufd[ir][it]=(0.0,0.0);
            ufuw[ir][it]=(0.0,0.0);
            ufuwg[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
         {
            uu[ir][it]=0.0;
            ud[ir][it]=0.0;
            uuwlrm[ir][it]=0.0;
            uuwlrmg[ir][it]=0.0;
         }
   
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

    for(ir=boul;ir<nr+boul;ir++)
       for(it=0;it<lt;it++)
          u[ir][it]=uraw[ir-boul][it];

    for(ir=0;ir<nr;ir++)
       {
         for(it=0;it<ltt;it++)
           {
              in1[it][0]=u[ir+boul][it];
              in1[it][1]=0.0;
           }

          fftwf_execute(p1);

          for(it=0;it<ltt;it++)
            {
              uf[ir+boul][it].real()=out1[it][0]/ltt;
              uf[ir+boul][it].imag()=out1[it][1]/ltt;
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

/*
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
               uf[ir][it].real()=out3[ir+boul][0];
               uf[ir][it].imag()=out3[ir+boul][1];
             }
         }

       for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=uf[ir][it].real();
                 in4[it][1]=uf[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             uu[ir][it]=out4[it][0];
         }

       return 0;
*/

     cout<<"Deghosting Begins..."<<endl;

     for(ir=0;ir<nrb;ir++)
//     for(ir=0;ir<nrb/2+1;ir++)
        {
           for(it=0;it<ltt/2+1;it++)
              {
                 if(it<ifmin)
                   {
                      ukfu[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfd[ir][it].real()=0.0;//downward wavefield
                      ukfd[ir][it].imag()=0.0;
                      ukfuw[ir][it].real()=0.0;//upward wavefield
                      ukfuw[ir][it].imag()=0.0;//upward wavefield
                      ukfuwg[ir][it].real()=0.0;
                      ukfuwg[ir][it].imag()=0.0;
                   }
                 else if(it>ifmax&&it<ltt/2+1)
                    {
                      ukfu[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfd[ir][it].real()=0.0;//downward wavefield
                      ukfd[ir][it].imag()=0.0;
                      ukfuw[ir][it].real()=0.0;
                      ukfuw[ir][it].imag()=0.0;
                      ukfuwg[ir][it].real()=0.0;
                      ukfuwg[ir][it].imag()=0.0;
                   }
                 else
                   {
                      if(pow(omega[it]/v,2)-pow(kx[ir],2)<0.0)
                        {
                           ukfu[ir][it].real()=0.0;
                           ukfu[ir][it].imag()=0.0;
                           ukfd[ir][it].real()=0.0;
                           ukfd[ir][it].imag()=0.0;
                           ukfuw[ir][it].real()=0.0;
                           ukfuw[ir][it].imag()=0.0;
                           ukfuwg[ir][it].real()=0.0;
                           ukfuwg[ir][it].imag()=0.0;
                        }
                      else if(pow(omega[it]/v,2)-pow(kx[ir],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                        {
                            kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));
                             ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                             ps1.imag()=-2*sdep*kz;
                             ps1=exp(ps1);
                             ps2.real()=1.0-ps1.real()+lamda;
                             ps2.imag()=-ps1.imag();
 
                             ukfu[ir][it]=ukf[ir][it]/ps2;   //upward field
/*
                             ps3.real()=0.0;               //phase-shift extrapolating operator to achieve downward field
                             ps3.imag()=2*sdep*kz;
                             ps3=exp(ps3);
                             ps4.real()=1.0-ps3.real()+lamda;
                             ps4.imag()=-ps3.imag();
*/                             
                             ukfd[ir][it]=ukf[ir][it]*ps1;  //downward field

                             ukfd[ir][it].real()*=-1.0;
                             ukfd[ir][it].imag()*=-1.0;

                             ps5.real()=0.0;
                             ps5.imag()=-2*wd*kz;
                             ps5=exp(ps5);
                             ukfuw[ir][it]=ukfu[ir][it]*ps5;

                             ps6.real()=0.0;
                             ps6.imag()=-2*(wd+sdep)*kz;
                             ps6=exp(ps6);
                             ukfuwg[ir][it]=ukfu[ir][it]*ps6;

                        }
                       else
                        {
                           ukfu[ir][it].real()=0.0;
                           ukfu[ir][it].imag()=0.0;
                           ukfd[ir][it].real()=0.0;
                           ukfd[ir][it].imag()=0.0;
                           ukfuw[ir][it].real()=0.0;
                           ukfuw[ir][it].imag()=0.0;
                           ukfuwg[ir][it].real()=0.0;
                           ukfuwg[ir][it].imag()=0.0;
                        }
                   }

                }

            for(it=ltt/2+1;it<ltt;it++)
              {
                 ukfu[ir][it].real()=ukfu[ir][ltt-it].real();
                 ukfu[ir][it].imag()=-ukfu[ir][ltt-it].imag();
                 ukfd[ir][it].real()=ukfd[ir][ltt-it].real();
                 ukfd[ir][it].imag()=-ukfd[ir][ltt-it].imag();
                 ukfuw[ir][it].real()=ukfuw[ir][ltt-it].real();
                 ukfuw[ir][it].imag()=-ukfuw[ir][ltt-it].imag();
                 ukfuwg[ir][it].real()=ukfuwg[ir][ltt-it].real();
                 ukfuwg[ir][it].imag()=-ukfuwg[ir][ltt-it].imag();
              }
        }

/*
      for(ir=nrb/2+1;ir<nrb;ir++)
       for(it=0;it<ltt;it++)
              {
                  ukfu[ir][it].real()=ukfu[nrb-ir][it].real();
                  ukfu[ir][it].imag()=-ukfu[nrb-ir][it].imag();
                  ukfd[ir][it].real()=ukfd[nrb-ir][it].real();
                  ukfd[ir][it].imag()=-ukfd[nrb-ir][it].imag();
                  ukfuw[ir][it].real()=ukfuw[nrb-ir][it].real();
                  ukfuw[ir][it].imag()=-ukfuw[nrb-ir][it].imag();
                  ukfuwg[ir][it].real()=ukfuwg[nrb-ir][it].real();
                  ukfuwg[ir][it].imag()=-ukfuwg[nrb-ir][it].imag();

              }
*/

     
       cout<<"IFFT Begins..."<<endl;

       for(it=0;it<ltt;it++)
         {   
           for(ir=0;ir<nrb;ir++)
              {   
                 in3[ir][0]=ukfu[ir][it].real();
                 in3[ir][1]=ukfu[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               ufu[ir][it].real()=out3[ir+boul][0];
               ufu[ir][it].imag()=out3[ir+boul][1];
             }   
         }   

       for(ir=0;ir<nr;ir++)
         {   
           for(it=0;it<ltt;it++)
              {   
                 in4[it][0]=ufu[ir][it].real();
                 in4[it][1]=ufu[ir][it].imag();
              }   
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             uu[ir][it]=out4[it][0];
         }   

       for(it=0;it<ltt;it++)
         {
           for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukfd[ir][it].real();
                 in3[ir][1]=ukfd[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               ufd[ir][it].real()=out3[ir+boul][0];
               ufd[ir][it].imag()=out3[ir+boul][1];
             }   
         }   
   
       for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufd[ir][it].real();
                 in4[it][1]=ufd[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             ud[ir][it]=-out4[it][0];
         }

       for(it=0;it<ltt;it++)
         {
           for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukfuw[ir][it].real();
                 in3[ir][1]=ukfuw[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               ufuw[ir][it].real()=out3[ir+boul][0];
               ufuw[ir][it].imag()=out3[ir+boul][1];
             }   
         }   
   
       for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufuw[ir][it].real();
                 in4[it][1]=ufuw[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             uuwlrm[ir][it]=-out4[it][0];
         }

       for(it=0;it<ltt;it++)
         {
           for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukfuwg[ir][it].real();
                 in3[ir][1]=ukfuwg[ir][it].imag();
              }
           fftwf_execute(p3);

           for(ir=0;ir<nr;ir++)
             {
               ufuwg[ir][it].real()=out3[ir+boul][0];
               ufuwg[ir][it].imag()=out3[ir+boul][1];
             }
         }

       for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufuwg[ir][it].real();
                 in4[it][1]=ufuwg[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             uuwlrmg[ir][it]=out4[it][0];
         }
 
    return 0;


}

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda;
   float wd;
 
   ifstream swq;
   swq.open("wlrm_plus_ghost_vb.par");
   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>lamda>>wd;
   swq.close();

   cout<<"Fna of input shot gather is===="<<fn1<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting is===="<<fn2<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting is===="<<fn3<<endl; 
   cout<<"Fna of source-side wlrm is===="<<fn4<<endl; 
   cout<<"Fna of source-side wlrm and ghost is===="<<fn5<<endl; 
   cout<<"Fna of ghost is===="<<fn6<<endl; 
   cout<<"Fna of wlrm is===="<<fn7<<endl; 
   cout<<"Fna of wlrm and ghost is===="<<fn8<<endl; 

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
   cout<<"Water Depth is===="<<wd<<endl;

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

   cout<<"ifmax, ifmin"<<ifmax<<" , "<<ifmin<<endl;

   cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

   for(is=0;is<nrb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dr*nrb);
   for(is=nrb/2+1;is<nrb;is++)
      kx[is]=-2*pai*1.0/float(2*dr)+2*pai*float(is-nrb/2)/float(dr*nrb);

//int source_deghosting_plus_wlrm(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu, complex <float> ** ukfd,  complex <float> **ufu, complex <float> **ufd,  float **uu, float **ud,  complex <float> **ukfuw, complex <float> ** ukfuwg,  complex <float> **ufuw, complex <float> **ufuwg, float **uuwlrm, float **uuwlrmg,  float *omega, float *kx, float sdep,float wd, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, float ifmin, float ifmax, float lamda)

   float **uraw;
   uraw=alloc2float(lt,nr);

   float **u;
   u=alloc2float(ltt,nrb);

   complex<float> **uf;
   uf=alloc2complex(ltt,nrb); 

   complex<float> **ukf;
   ukf=alloc2complex(ltt,nrb);

   complex<float> **ukfu;
   ukfu=alloc2complex(ltt,nrb);

   complex<float> **ukfd;
   ukfd=alloc2complex(ltt,nrb);

    complex<float> **ufu;
   ufu=alloc2complex(ltt,nr);

   complex<float> **ufd;
   ufd=alloc2complex(ltt,nr);

    complex<float> **ukfuw;
   ukfuw=alloc2complex(ltt,nrb);

   complex<float> **ukfuwg;
   ukfuwg=alloc2complex(ltt,nrb);

    complex<float> **ufuw;
   ufuw=alloc2complex(ltt,nr);

   complex<float> **ufuwg;
   ufuwg=alloc2complex(ltt,nr);
 
   float **uu;
   uu=alloc2float(lt,nr);

   float **ud;
   ud=alloc2float(lt,nr);

   float **uuwlrm;
   uuwlrm=alloc2float(lt,nr);

   float **uuwlrmg;
   uuwlrmg=alloc2float(lt,nr);
 
  ifstream swq1;
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

  ofstream swq4;
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

  ofstream swq8;
  swq8.open(fn8,ios::binary);
  if(!swq8)
       cout<<"cannot open "<<fn8<<endl;

  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
            swq1.read((char*)&uraw[ir][it],sizeof(uraw[ir][it]));

//int source_deghosting_plus_wlrm(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu, complex <float> ** ukfd,  complex <float> **ufu, complex <float> **ufd,  float **uu, float **ud,  complex <float> **ukfuw, complex <float> ** ukfuwg,  complex <float> **ufuw, complex <float> **ufuwg, float **uuwlrm, float **uuwlrmg,  float *omega, float *kx, float sdep,float wd, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, float ifmin, float ifmax, float lamda)

       source_deghosting_plus_wlrm(uraw, u, uf, ukf, ukfu,  ukfd, ufu, ufd, uu, ud, ukfuw,  ukfuwg,  ufuw, ufuwg, uuwlrm, uuwlrmg,  omega, kx, sdep, wd,  ds,  dr,  nr,  lt,  ltt,  boul, bour,   nrb,  v, theta, ifmin,  ifmax, lamda);

       for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
         {   
           swq2.write((char*)&uu[ir][it],sizeof(uu[ir][it]));
           swq3.write((char*)&ud[ir][it],sizeof(ud[ir][it]));
           swq4.write((char*)&uuwlrm[ir][it],sizeof(uuwlrm[ir][it]));
           swq5.write((char*)&uuwlrmg[ir][it],sizeof(uuwlrmg[ir][it]));
         }

        cout<<is+1<<" shot source deghosting done..."<<endl;

        cout<<is+1<<" shot done!"<<endl;
    }
  
        return 0; 

}
































