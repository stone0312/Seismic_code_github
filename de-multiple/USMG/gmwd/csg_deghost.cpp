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

int deghosting_plus_extrapolation(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu,  complex <float> **ukfg, complex <float> **ukfm, complex<float> **ukfmg, complex <float> **ufu, complex <float> **ufg, complex <float> **ufm, complex<float>**ufmg,  float **uu, float **ug, float **um, float **umg, float **mg, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, float ifmin, float ifmax, float lamda, int flag, float wd)
{
    float kz;
    complex <float> ps1;
    complex <float> ps2;
    complex <float> ps3;
    complex <float> ps4;
    complex <float> ps5;
    complex <float> ps6;

    int ir,is, it;
    int tmp;
    tmp=0;

//zeroing the arrays
    for(ir=0;ir<nrb;ir++) 
       for(it=0;it<ltt;it++)
         {
            u[ir][it]=0.0; 
                      uf[ir][it].real()=0.0;//upward wavefield
                      uf[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfu[ir][it].imag()=0.0;
                      ukfg[ir][it].real()=0.0;//downward wavefield
                      ukfg[ir][it].imag()=0.0;
                      ukfm[ir][it].real()=0.0;//wlrm
                      ukfm[ir][it].imag()=0.0;
                      ukfmg[ir][it].real()=0.0;//ghosts of wlrm
                      ukfmg[ir][it].imag()=0.0;
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
         {
                      ufu[ir][it].real()=0.0;//upward wavefield
                      ufu[ir][it].imag()=0.0;
                      ufg[ir][it].real()=0.0;//downward wavefield
                      ufg[ir][it].imag()=0.0;
                      ufm[ir][it].real()=0.0;//wlrm
                      ufm[ir][it].imag()=0.0;
                      ufmg[ir][it].real()=0.0;//ghosts of wlrm
                      ufmg[ir][it].imag()=0.0;
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
         {
            uu[ir][it]=0.0;
            ug[ir][it]=0.0;
            um[ir][it]=0.0;
            umg[ir][it]=0.0;
         }

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
        for(ir=nrb/2+1;ir<nrb;ir++) 
           for(it=0;it<ltt;it++)
             ukf[ir][it].imag()=-ukf[nrb-ir][it].imag();
        for(it=0;it<ltt;it++)
         {
           for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukf[ir][it].real();
                 in3[ir][1]=ukf[ir][it].imag();
//                 in3[ir][1]=0.0;
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
//                 in4[it][1]=0.0;
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             uu[ir][it]=out4[it][0];
         }
    ofstream swq22;
    swq22.open("test.dat",ios::binary);
    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
           swq22.write((char*)&uu[ir][it],sizeof(uu[ir][it]));

    return 0;
*/
   if(flag==1)
    {//1111
     cout<<"Deghosting Begins..."<<endl;

     for(ir=0;ir<nrb;ir++)
        {
           for(it=0;it<ltt;it++)
              {
                      if(pow(omega[it]/v,2)-pow(kx[ir],2)<0.0)
                        {
                           ukfu[ir][it].real()=0.0;
                           ukfu[ir][it].imag()=0.0;
                           ukfg[ir][it].real()=0.0;
                           ukfg[ir][it].imag()=0.0;
                           ukfm[ir][it].real()=0.0;//wlrm
                           ukfm[ir][it].imag()=0.0;
                           ukfmg[ir][it].real()=0.0;//ghosts of wlrm
                           ukfmg[ir][it].imag()=0.0;
                        }
                      else
                        {
                             kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));
                             ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                             ps1.imag()=-2*rdep*kz;
                             ps1=exp(ps1);
                             ps2.real()=1.0-ps1.real()+lamda;
                             ps2.imag()=-ps1.imag();

                             ps3.real()=0.0;               //phase-shift extrapolating operator to achieve downward field
                             ps3.imag()=2*rdep*kz;
                             ps3=exp(ps3);
                             ps4.real()=1.0-ps3.real()+lamda;
                             ps4.imag()=-ps3.imag();

                             ukfu[ir][it]=ukf[ir][it]/ps2;   //upward field
                             ukfg[ir][it]=ukf[ir][it]/ps4;  //downward field
                       }
                }

         }

//ifft

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


        for(it=0;it<ltt;it++)
         {
           for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukfg[ir][it].real();
                 in3[ir][1]=ukfg[ir][it].imag();
              }
           fftwf_execute(p3);

           for(ir=0;ir<nr;ir++)
             {
               ufg[ir][it].real()=out3[ir+boul][0];
               ufg[ir][it].imag()=out3[ir+boul][1];
             }
         }

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
               ufm[ir][it].real()=out3[ir+boul][0];
               ufm[ir][it].imag()=out3[ir+boul][1];
             }
         }

        for(it=0;it<ltt;it++)
         {
           for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukfmg[ir][it].real();
                 in3[ir][1]=ukfmg[ir][it].imag();
              }
           fftwf_execute(p3);

           for(ir=0;ir<nr;ir++)
             {
               ufmg[ir][it].real()=out3[ir+boul][0];
               ufmg[ir][it].imag()=out3[ir+boul][1];
             }
         }

//        cout<<"2222"<<endl;

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

        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufg[ir][it].real();
                 in4[it][1]=ufg[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             ug[ir][it]=out4[it][0];
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
        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufmg[ir][it].real();
                 in4[it][1]=ufmg[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             umg[ir][it]=out4[it][0];
         }

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
           mg[ir][it]=um[ir][it]+umg[ir][it]+ug[ir][it];
//           mg[ir][it]=um[ir][it]+umg[ir][it];


    }//1111


    return 0;

}

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fn12[256],fn13[256],fn14[256],fn15[256],fn16[256],fn17[256],fn18[256],fn19[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda,n;
   float wd;
 
   ifstream swq;
   swq.open("csg_deghost.par");
   swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>lamda>>n>>wd;
   swq.close();

   cout<<"Fna of input shot gather is===="<<fn1<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting  is===="<<fn2<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting  is===="<<fn3<<endl; 
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

   cout<<"ifmax, ifmin===="<<ifmax<<" , "<<ifmin<<endl;

   cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

   for(is=0;is<nrb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dr*nrb);
   for(is=nrb/2+1;is<nrb;is++)
      kx[is]=-2*pai*1.0/float(2*dr)+2*pai*float(is-nrb/2)/float(dr*nrb);

//int deghosting_plus_extrapolation(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu,  complex <float> **ukfg, complex <float> **ukfm, complex<float> **ukfmg, complex <float> **ufu, complex <float> **ufg, complex <float> **ufm, complex<float>**ufmg,  float **uu, float **ug, float **um, float **umg, float **mg, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, float ifmin, float ifmax, float lamda, int flag)

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

   complex<float> **ukfg;
   ukfg=alloc2complex(ltt,nrb);

    complex<float> **ukfm;
   ukfm=alloc2complex(ltt,nrb);

   complex<float> **ukfmg;
   ukfmg=alloc2complex(ltt,nrb); 

    complex<float> **ufu;
   ufu=alloc2complex(ltt,nr);

   complex<float> **ufg;
   ufg=alloc2complex(ltt,nr);

    complex<float> **ufm;
   ufm=alloc2complex(ltt,nr);

   complex<float> **ufmg;
   ufmg=alloc2complex(ltt,nr);

   float **uu;
   uu=alloc2float(lt,nr);

   float **ug;
   ug=alloc2float(lt,nr);

   float **um;
   um=alloc2float(lt,nr);

   float **umg;
   umg=alloc2float(lt,nr);
 
   float **mg;
   mg=alloc2float(lt,nr);

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

  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
            swq1.read((char*)&uraw[ir][it],sizeof(uraw[ir][it]));

       flag=1;

       deghosting_plus_extrapolation(uraw,u, uf, ukf, ukfu, ukfg, ukfm, ukfmg, ufu, ufg, ufm, ufmg,uu,ug, um, umg, mg, omega, kx, sdep, rdep,  ds, dr, nr, lt,  ltt,  boul, bour,  nrb,  v, theta,ifmin, ifmax,  lamda,  flag,wd);

      for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
         {  
           swq2.write((char*)&uu[ir][it],sizeof(uu[ir][it]));
           swq3.write((char*)&ug[ir][it],sizeof(ug[ir][it]));
         }

       cout<<is+1<<" shot source-side  deghosting and wlrm plus ghosts prediction done..."<<endl;

/*
       flag=2;

        deghosting_plus_extrapolation(uraw, u, uf, ukf, ukfu, ukfd, ukfdfe, ukfube, ufu, ufd, ufdfe, ufube, uu, ud, udfe, uube, omega, kx,  sdep, rdep, ds, dr, nr, lt, ltt, boul, bour, nrb, v, theta, ifmin, ifmax, lamda, flag);

       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
            {
               swq14.write((char*)&ufu[ir][it].real(),sizeof(ufu[ir][it].real()));
               swq15.write((char*)&ufu[ir][it].imag(),sizeof(ufu[ir][it].imag()));
               swq16.write((char*)&ufd[ir][it].real(),sizeof(ufd[ir][it].real()));
               swq17.write((char*)&ufd[ir][it].imag(),sizeof(ufd[ir][it].imag()));
            } 

       for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
         {   
           swq18.write((char*)&uu[ir][it],sizeof(uu[ir][it]));
           swq19.write((char*)&ud[ir][it],sizeof(ud[ir][it]));
         }

        cout<<is+1<<" shot receiver deghosting done..."<<endl;
*/
        cout<<is+1<<" shot done!"<<endl;
    }
  
        return 0; 

}
































