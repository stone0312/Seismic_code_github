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

//zeroing the arrays
    for(ir=0;ir<nrb;ir++) 
       for(it=0;it<ltt;it++)
         {
            u[ir][it]=0.0; 
            uf[ir][it]=(0.0,0.0);
            ukf[ir][it]=(0.0,0.0);
            ukfu[ir][it]=(0.0,0.0);
            ukfg[ir][it]=(0.0,0.0);
            ukfm[ir][it]=(0.0,0.0);
            ukfmg[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
         {
            ufu[ir][it]=(0.0,0.0);
            ufg[ir][it]=(0.0,0.0);
            ufm[ir][it]=(0.0,0.0);
            ufmg[ir][it]=(0.0,0.0);
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



    for(ir=0;ir<nrb;ir++)
       for(it=0;it<ltt;it++)
          uf[ir][it]=(0.0,0.0); 

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

//    cout<<"22222"<<endl;  
 
    for(it=0;it<ltt;it++)
       {
           for(ir=0;ir<nrb;ir++)
              {
                 in2[ir][0]=uf[ir][it].real();
                 in2[ir][1]=uf[ir][it].imag();
              } 
           fftwf_execute(p2); 
          
//           cout<<it<<"   22222"<<endl;  

           for(ir=0;ir<nrb;ir++)
             {
               ukf[ir][it].real()=out2[ir][0]/nrb;
               ukf[ir][it].imag()=out2[ir][1]/nrb;
             }
       }

   if(flag==1)
    {//1111
     cout<<"Deghosting Begins..."<<endl;

     for(ir=0;ir<nrb;ir++)
        {
           for(it=0;it<ltt/2+1;it++)
              {
                 if(it<ifmin)
                   {
                      ukfu[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfg[ir][it].real()=0.0;//downward wavefield
                      ukfg[ir][it].imag()=0.0;
                      ukfm[ir][it].real()=0.0;//wlrm
                      ukfm[ir][it].imag()=0.0;
                      ukfmg[ir][it].real()=0.0;//ghosts of wlrm
                      ukfmg[ir][it].imag()=0.0;
                   }
                 else if(it>ifmax&&it<ltt/2+1)
                    {
                      ukfu[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfg[ir][it].real()=0.0;//downward wavefield
                      ukfg[ir][it].imag()=0.0;
                      ukfm[ir][it].real()=0.0;//wlrm
                      ukfm[ir][it].imag()=0.0;
                      ukfmg[ir][it].real()=0.0;//ghosts of wlrm
                      ukfmg[ir][it].imag()=0.0;
                   }
                 else
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
                      else if(pow(omega[it]/v,2)-pow(kx[ir],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                        {
                             kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));
                             ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                             ps1.imag()=-2*wd*kz;
                             ps1=exp(ps1);

                             ukfmg[ir][it]=ukf[ir][it]*ps1;   //upward field

                             ukfmg[ir][it].real()*=-1.0;
                             ukfmg[ir][it].imag()*=-1.0;


                             ukfmg[ir][ltt-it].real()=ukfmg[ir][it].real();
                             ukfmg[ir][ltt-it].imag()=-ukfmg[ir][it].imag();
 
                        }
                       else
                        {
                           ukfu[ir][it].real()=0.0;
                           ukfu[ir][it].imag()=0.0;
                           ukfg[ir][it].real()=0.0;
                           ukfg[ir][it].imag()=0.0;
                           ukfm[ir][it].real()=0.0;
                           ukfm[ir][it].imag()=0.0;
                           ukfmg[ir][it].real()=0.0;
                           ukfmg[ir][it].imag()=0.0;
                        } 
                   }
//                 if(it%1000==0)
//                    cout<<it<<" Frequency Slices Deghosting Done..."<<endl;
                }
           

            for(it=ltt/2+1;it<ltt;it++)
              {
                 ukfu[ir][it].real()=ukfu[ir][ltt-it].real();
                 ukfu[ir][it].imag()=-ukfu[ir][ltt-it].imag();
                 ukfg[ir][it].real()=ukfg[ir][ltt-it].real();
                 ukfg[ir][it].imag()=-ukfg[ir][ltt-it].imag();
                 ukfm[ir][it].real()=ukfm[ir][ltt-it].real();
                 ukfm[ir][it].imag()=-ukfm[ir][ltt-it].imag();
                 ukfmg[ir][it].real()=ukfmg[ir][ltt-it].real();
                 ukfmg[ir][it].imag()=-ukfmg[ir][ltt-it].imag();
              }
        }

//ifft

     cout<<"IFFT Begins..."<<endl;
 
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
           mg[ir][it]=umg[ir][it];

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
   swq.open("receiver_ghosts_obc.par");
   swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>lamda>>n>>wd;
   swq.close();

   cout<<"Fna of upward wavefields after source-side deghosting  is===="<<fn2<<endl; 
   cout<<"Fna of source-side ghosts and wlrm is===="<<fn4<<endl; 
   cout<<"Fna of all ghosts and wlrm  is===="<<fn4<<endl; 

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

   float **uraw1;
   uraw1=alloc2float(lt,nr);
   float **uraw2;
   uraw2=alloc2float(lt,nr);

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

   float **mg_tot;
   mg_tot=alloc2float(lt,nr);

  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl;


  ifstream swq2;
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
            swq1.read((char*)&uraw1[ir][it],sizeof(uraw1[ir][it]));
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
            swq2.read((char*)&uraw2[ir][it],sizeof(uraw2[ir][it]));

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
             uraw[ir][it]=uraw1[ir][it]+uraw2[ir][it];
 
       flag=1;

       deghosting_plus_extrapolation(uraw,u, uf, ukf, ukfu, ukfg, ukfm, ukfmg, ufu, ufg, ufm, ufmg,uu,ug, um, umg, mg, omega, kx, sdep, rdep,  ds, dr, nr, lt,  ltt,  boul, bour,  nrb,  v, theta,ifmin, ifmax,  lamda,  flag,wd);

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
           //  mg_tot[ir][it]=mg[ir][it]+uraw2[ir][it];
             mg_tot[ir][it]=mg[ir][it];

      for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
           swq3.write((char*)&mg_tot[ir][it],sizeof(mg_tot[ir][it]));

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
































