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

int deghosting_plus_extrapolation(float **urawu, float **uu, complex <float> **ufu, complex <float> **ukfu,float **urawd, float **ud, complex <float> **ufd, complex <float> **ukfd, complex <float> **ukfm, complex <float> ** ufm, float **um, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, int ifmin, int ifmax, float wd, int n)
{
    float kz;
    complex <float> ps1;
    complex <float> ps2;
    complex <float> ps3;

    int ir,is,it;

//zeroing the arrays
    for(ir=0;ir<nrb;ir++) 
       for(it=0;it<ltt;it++)
         {
            uu[ir][it]=0.0; 
            ufu[ir][it]=(0.0,0.0);
            ukfu[ir][it]=(0.0,0.0);
            ud[ir][it]=0.0;
            ufd[ir][it]=(0.0,0.0);
            ukfd[ir][it]=(0.0,0.0);
            ukfm[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
            ufm[ir][it]=(0.0,0.0);

     for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
            um[ir][it]=0.0;

    for(ir=boul+n;ir<nr+boul+n;ir++)
      for(it=0;it<lt;it++)
        {
           uu[ir][it]=urawu[ir-boul-n][it];
           ud[ir][it]=urawd[ir-boul-n][it];
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
         {
           ufu[ir][it]=(0.0,0.0); 
           ufd[ir][it]=(0.0,0.0); 
         }

    for(ir=0;ir<nrb;ir++)
       {
         for(it=0;it<ltt;it++)
           {
              in1[it][0]=uu[ir][it];
              in1[it][1]=0.0;
           }

          fftwf_execute(p1); 

          for(it=0;it<ltt;it++)
            {
              ufu[ir][it].real()=out1[it][0]/ltt;
              ufu[ir][it].imag()=out1[it][1]/ltt;
            }
       }  

    for(it=0;it<ltt;it++)
       {
           for(ir=0;ir<nrb;ir++)
              {
                 in2[ir][0]=ufu[ir][it].real();
                 in2[ir][1]=ufu[ir][it].imag();
              } 
           fftwf_execute(p2); 
          

           for(ir=0;ir<nrb;ir++)
             {
               ukfu[ir][it].real()=out2[ir][0]/nrb;
               ukfu[ir][it].imag()=out2[ir][1]/nrb;
             }
       }


    for(ir=0;ir<nrb;ir++)
       {
         for(it=0;it<ltt;it++)
           {
              in1[it][0]=ud[ir][it];
              in1[it][1]=0.0;
           }

          fftwf_execute(p1);

          for(it=0;it<ltt;it++)
            {
              ufd[ir][it].real()=out1[it][0]/ltt;
              ufd[ir][it].imag()=out1[it][1]/ltt;
            }
       }

    for(it=0;it<ltt;it++)
       {
           for(ir=0;ir<nrb;ir++)
              {
                 in2[ir][0]=ufd[ir][it].real();
                 in2[ir][1]=ufd[ir][it].imag();
              }
           fftwf_execute(p2);


           for(ir=0;ir<nrb;ir++)
             {
               ukfd[ir][it].real()=out2[ir][0]/nrb;
               ukfd[ir][it].imag()=out2[ir][1]/nrb;
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

                             ps2.real()=0.0;
                             ps2.imag()=-2*(wd-rdep)*kz;
                             ps2=exp(ps2);

                             ps3.real()=0.0;
                             ps3.imag()=-2*(wd+rdep)*kz;
                             ps3=exp(ps3);

                             ukfm[ir][it]=ukfd[ir][it]*ps2-ukfd[ir][it]*ps1+ukfu[ir][it]*ps3-ukfu[ir][it]*ps1; //1 

//                             ukfm[ir][it]=ukfd[ir][it]*ps2-ukfd[ir][it]*ps1;                       //2
                                                  
//                             ukfm[ir][it]=ukfu[ir][it]*ps3-ukfu[ir][it]*ps1;                       //3

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
                 in3[ir][0]=ukfm[ir][it].real();
                 in3[ir][1]=ukfm[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               ufm[ir][it].real()=out3[ir+boul+n][0];
               ufm[ir][it].imag()=out3[ir+boul+n][1];
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
             um[ir][it]=out4[it][0];
         }


    return 0;

}

int main()
{
   char fn1[256],fn11[256],fn2[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb, n;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda,wd, min_off;
 
   ifstream swq;
   swq.open("wlrm_plus_ghost_prediction.par");
   swq>>fn1>>fn11>>fn2>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>lamda>>wd>>min_off;
   swq.close();

   cout<<"Fna of input upwave shot gather is===="<<fn1<<endl; 
   cout<<"Fna of input downwave shot gather is===="<<fn11<<endl; 
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

   n=int(min_off/dr+0.5);
   nrb=boul+n+nr+bour;

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

   float **urawu;
   urawu=alloc2float(lt,nr);

   float **urawd;
   urawd=alloc2float(lt,nr);
   
   float **uu;
   uu=alloc2float(ltt,nrb);

   complex<float> **ufu;
   ufu=alloc2complex(ltt,nrb); 

   complex<float> **ukfu;
   ukfu=alloc2complex(ltt,nrb);

   float **ud;
   ud=alloc2float(ltt,nrb);

   complex<float> **ufd;
   ufd=alloc2complex(ltt,nrb); 

   complex<float> **ukfd;
   ukfd=alloc2complex(ltt,nrb);

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

  ifstream swq11;
  swq11.open(fn11,ios::binary);
  if(!swq11)
       cout<<"cannot open "<<fn11<<endl;

  ofstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;
 
  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
            swq1.read((char*)&urawu[ir][it],sizeof(urawu[ir][it]));

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
            swq11.read((char*)&urawd[ir][it],sizeof(urawd[ir][it]));

       deghosting_plus_extrapolation(urawu, uu, ufu, ukfu, urawd, ud, ufd, ukfd,  ukfm, ufm, um, omega, kx,  sdep, rdep, ds, dr, nr, lt, ltt, boul, bour, nrb, v, theta, ifmin, ifmax, wd, n );

       for(ir=n;ir<nr;ir++)
        for(it=0;it<lt;it++)
           swq2.write((char*)&um[ir][it],sizeof(um[ir][it]));

        cout<<is+1<<" shot done!"<<endl;
    }
  
        return 0; 

}
































