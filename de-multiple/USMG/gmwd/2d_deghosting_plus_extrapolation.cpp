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

int deghosting_plus_extrapolation(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu, complex <float> ** ukfd, complex <float> **ukfdfe, complex <float> **ukfube, complex <float> **ufu, complex <float> **ufd, complex <float> **ufdfe, complex <float> **ufube, float **uu, float **ud, float **udfe, float **uube, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, float ifmin, float ifmax, float lamda, int flag)
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
            ukfd[ir][it]=(0.0,0.0);
            ukfdfe[ir][it]=(0.0,0.0);
            ukfube[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
         {
            ufu[ir][it]=(0.0,0.0);
            ufd[ir][it]=(0.0,0.0);
            ufdfe[ir][it]=(0.0,0.0);
            ufube[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
         {
            uu[ir][it]=0.0;
            ud[ir][it]=0.0;
            udfe[ir][it]=0.0;
            uube[ir][it]=0.0;
         }

    for(ir=boul;ir<nr+boul;ir++)
      for(it=0;it<lt;it++)
           u[ir][it]=uraw[ir-boul][it];


    float dep1,dep2;
 
    if(flag==1)
      dep1=sdep;
    else
      dep1=rdep;
    dep2=fabs(rdep-sdep); 

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

     cout<<"Deghosting Begins..."<<endl;

     for(ir=0;ir<nrb;ir++)
        {
           for(it=0;it<ltt/2+1;it++)
              {
                 if(it<ifmin)
                   {
                      ukfu[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfd[ir][it].real()=0.0;//downward wavefield
                      ukfd[ir][it].imag()=0.0;
                   }
                 else if(it>ifmax&&it<ltt/2+1)
                    {
                      ukfu[ir][it].real()=0.0;//upward wavefield
                      ukfu[ir][it].imag()=0.0;
                      ukfd[ir][it].real()=0.0;//downward wavefield
                      ukfd[ir][it].imag()=0.0;
                   }
                 else
                   {
                      if(pow(omega[it]/v,2)-pow(kx[ir],2)<0.0)
                        {
                           ukfu[ir][it].real()=0.0;
                           ukfu[ir][it].imag()=0.0;
                           ukfd[ir][it].real()=0.0;
                           ukfd[ir][it].imag()=0.0;
                        }
                      else if(pow(omega[it]/v,2)-pow(kx[ir],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                        {
                            kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));
                             ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                             ps1.imag()=-2*dep1*kz;
                             ps1=exp(ps1);
                             ps2.real()=1.0-ps1.real()+lamda;
                             ps2.imag()=-ps1.imag();

                             ps3.real()=0.0;               //phase-shift extrapolating operator to achieve downward field
                             ps3.imag()=2*dep1*kz;
                             ps3=exp(ps3);
                             ps4.real()=1.0-ps3.real()+lamda;
                             ps4.imag()=-ps3.imag();

                             ukfu[ir][it]=ukf[ir][it]/ps2;   //upward field
                             ukfd[ir][it]=ukf[ir][it]/ps4;  //downward field

                                  ps5.real()=0.0;
                                  ps5.imag()=-dep2*kz;
                                  ps5=exp(ps5);
                                  ukfdfe[ir][it]=ukfd[ir][it]*ps5;
                                  
//                                  ukfdfe[ir][ltt-it].real()=ukfdfe[ir][it].real();
//                                  ukfdfe[ir][ltt-it].imag()=-ukfdfe[ir][it].imag();

                                  ps6.real()=0.0;
                                  ps6.imag()=dep2*kz;
                                  ps6=exp(ps6);
                                  ukfube[ir][it]=ukfd[ir][it]*ps6;

//                                  ukfube[ir][ltt-it].real()=ukfube[ir][it].real();
//                                  ukfube[ir][ltt-it].imag()=-ukfube[ir][it].imag();
 
                        }
                       else
                        {
                           ukfu[ir][it].real()=0.0;
                           ukfu[ir][it].imag()=0.0;
                           ukfd[ir][it].real()=0.0;
                           ukfd[ir][it].imag()=0.0;
                        } 
                   }
                }
           

            for(it=ltt/2+1;it<ltt;it++)
              {
                 ukfu[ir][it].real()=ukfu[ir][ltt-it].real();
                 ukfu[ir][it].imag()=-ukfu[ir][ltt-it].imag();
                 ukfd[ir][it].real()=ukfd[ir][ltt-it].real();
                 ukfd[ir][it].imag()=-ukfd[ir][ltt-it].imag();
                 ukfdfe[ir][it].real()=ukfdfe[ir][ltt-it].real();
                 ukfdfe[ir][it].imag()=-ukfdfe[ir][ltt-it].imag();
                 ukfube[ir][it].real()=ukfube[ir][ltt-it].real();
                 ukfube[ir][it].imag()=-ukfube[ir][ltt-it].imag();
              }
        }

//ifft

     cout<<"IFFT Begins..."<<endl;
 
     if(flag!=1)
      {//1111
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
                 in4[it][0]=ufd[ir][it].real();
                 in4[it][1]=ufd[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             ud[ir][it]=-out4[it][0];
         }
       }//
    
      else
        {//2222
          
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
                 in3[ir][0]=ukfdfe[ir][it].real();
                 in3[ir][1]=ukfdfe[ir][it].imag();
              }
            fftwf_execute(p3);

           for(ir=0;ir<nr;ir++)
             {
               ufdfe[ir][it].real()=out3[ir+boul][0];
               ufdfe[ir][it].imag()=out3[ir+boul][1];
             }
           }

//         cout<<"22222"<<endl;

          for(it=0;it<ltt;it++)
           {
             for(ir=0;ir<nrb;ir++)
              {
                 in3[ir][0]=ukfube[ir][it].real();
                 in3[ir][1]=ukfube[ir][it].imag();
              }
           fftwf_execute(p3);

           for(ir=0;ir<nr;ir++)
             {
               ufube[ir][it].real()=out3[ir+boul][0];
               ufube[ir][it].imag()=out3[ir+boul][1];
             }
           }

          for(ir=0;ir<nr;ir++)
           {
             for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufdfe[ir][it].real();
                 in4[it][1]=ufdfe[ir][it].imag();
              }
             fftwf_execute(p4);

             for(it=0;it<lt;it++)
               udfe[ir][it]=out4[it][0];
           }

        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=ufube[ir][it].real();
                 in4[it][1]=ufube[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             uube[ir][it]=out4[it][0];
         }
       }//2222

    return 0;

}

//int deghosting_plus_extrapolation(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu, complex <float> ** ukfd, complex <float> **ukfdfe, complex <float> **ukfube, complex <float> **ufu, complex <float> **ufd, complex <float> **ufdfe, complex <float> **ufube, float **uu, float **ud, float **udfe, float **uube, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour, int nrb,  float v, float theta, float ifmin, float ifmax, float lamda, int flag)
int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fn12[256],fn13[256],fn14[256],fn15[256],fn16[256],fn17[256],fn18[256],fn19[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda;
 
   ifstream swq;
   swq.open("2d_deghosting_plus_extrapolation.par");
   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>fn11>>fn12>>fn13>>fn14>>fn15>>fn16>>fn17>>fn18>>fn19>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>lamda;
   swq.close();

   cout<<"Fna of input shot gather is===="<<fn1<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting real parts is===="<<fn2<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting imaginary parts is===="<<fn3<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting real parts is===="<<fn4<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting imaginary parts is===="<<fn5<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting and redatumming real parts is===="<<fn6<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting and redatumming imaginary parts is===="<<fn7<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting and redatumming real parts is===="<<fn8<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting and redatumming imaginary parts is===="<<fn9<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting is===="<<fn10<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting is===="<<fn11<<endl; 
   cout<<"Fna of upward wavefields after source-side deghosting and redatumming is===="<<fn12<<endl; 
   cout<<"Fna of downward wavefields after source-side deghosting and redatumming is===="<<fn13<<endl; 
   cout<<"Fna of upward wavefields after receiver-side deghosting real parts is===="<<fn14<<endl;
   cout<<"Fna of upward wavefields after receiver-side deghosting imaginary parts is===="<<fn15<<endl;
   cout<<"Fna of downward wavefields after receiver-side deghosting real parts is===="<<fn16<<endl;
   cout<<"Fna of downward wavefields after receiver-side deghosting imaginary parts is===="<<fn17<<endl;
   cout<<"Fna of upward wavefields after receiver-side deghosting is===="<<fn18<<endl; 
   cout<<"Fna of downward wavefields after receiver-side deghosting is===="<<fn19<<endl; 

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

   cout<<"ifmax, ifmin"<<ifmax<<" , "<<ifmin<<endl;

   cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

   for(is=0;is<nrb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dr*nrb);
   for(is=nrb/2+1;is<nrb;is++)
      kx[is]=-2*pai*1.0/float(2*dr)+2*pai*float(is-nrb/2)/float(dr*nrb);

//int deghosting_plus_extrapolation(float **uraw, float **u, complex <float> **uf, complex <float> **ukf, complex <float> **ukfu, complex <float> ** ukfd, complex <float> **ukfdfe, complex <float> **ukfube, complex <float> **ufu, complex <float> **ufd, complex <float> **ufdfe, complex <float> **ufube, float **uu, float **ud, float **udfe, float **uube, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,int nrb,  float v, float theta, float ifmin, float ifmax, float lamda, int flag)

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

    complex<float> **ukfube;
   ukfube=alloc2complex(ltt,nrb);

   complex<float> **ukfdfe;
   ukfdfe=alloc2complex(ltt,nrb); 

    complex<float> **ufu;
   ufu=alloc2complex(ltt,nr);

   complex<float> **ufd;
   ufd=alloc2complex(ltt,nr);

    complex<float> **ufube;
   ufube=alloc2complex(ltt,nr);

   complex<float> **ufdfe;
   ufdfe=alloc2complex(ltt,nr);

   float **uu;
   uu=alloc2float(lt,nr);

   float **ud;
   ud=alloc2float(lt,nr);

   float **uube;
   uube=alloc2float(lt,nr);

   float **udfe;
   udfe=alloc2float(lt,nr);
 
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

  ofstream swq9;
  swq9.open(fn9,ios::binary);
  if(!swq9)
       cout<<"cannot open "<<fn9<<endl;

  ofstream swq10;
  swq10.open(fn10,ios::binary);
  if(!swq10)
       cout<<"cannot open "<<fn10<<endl;

  ofstream swq11;
  swq11.open(fn11,ios::binary);
  if(!swq11)
       cout<<"cannot open "<<fn11<<endl;

  ofstream swq12;
  swq12.open(fn12,ios::binary);
  if(!swq12)
       cout<<"cannot open "<<fn12<<endl;

  ofstream swq13;
  swq13.open(fn13,ios::binary);
  if(!swq13)
       cout<<"cannot open "<<fn13<<endl;

  ofstream swq14;
  swq14.open(fn14,ios::binary);
  if(!swq14)
       cout<<"cannot open "<<fn14<<endl;

  ofstream swq15;
  swq15.open(fn15,ios::binary);
  if(!swq15)
       cout<<"cannot open "<<fn15<<endl;

  ofstream swq16;
  swq16.open(fn16,ios::binary);
  if(!swq16)
       cout<<"cannot open "<<fn16<<endl;

  ofstream swq17;
  swq17.open(fn17,ios::binary);
  if(!swq17)
       cout<<"cannot open "<<fn17<<endl;

  ofstream swq18;
  swq18.open(fn18,ios::binary);
  if(!swq18)
       cout<<"cannot open "<<fn18<<endl;

  ofstream swq19;
  swq19.open(fn19,ios::binary);
  if(!swq19)
       cout<<"cannot open "<<fn19<<endl;


  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++) 
            swq1.read((char*)&uraw[ir][it],sizeof(uraw[ir][it]));

       flag=1;

       deghosting_plus_extrapolation(uraw, u, uf, ukf, ukfu, ukfd, ukfdfe, ukfube, ufu, ufd, ufdfe, ufube, uu, ud, udfe, uube, omega, kx,  sdep, rdep, ds, dr, nr, lt, ltt, boul, bour, nrb, v, theta, ifmin, ifmax, lamda, flag);

       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
            {
               swq2.write((char*)&ufu[ir][it].real(),sizeof(ufu[ir][it].real()));
               swq3.write((char*)&ufu[ir][it].imag(),sizeof(ufu[ir][it].imag()));
               swq4.write((char*)&ufd[ir][it].real(),sizeof(ufd[ir][it].real()));
               swq5.write((char*)&ufd[ir][it].imag(),sizeof(ufd[ir][it].imag()));
               swq6.write((char*)&ufube[ir][it].real(),sizeof(ufube[ir][it].real()));
               swq7.write((char*)&ufube[ir][it].imag(),sizeof(ufube[ir][it].imag()));
               swq8.write((char*)&ufdfe[ir][it].real(),sizeof(ufdfe[ir][it].real()));
               swq9.write((char*)&ufdfe[ir][it].imag(),sizeof(ufdfe[ir][it].imag()));
            }  

       for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
         {  
           swq10.write((char*)&uu[ir][it],sizeof(uu[ir][it]));
           swq11.write((char*)&ud[ir][it],sizeof(ud[ir][it]));
           swq12.write((char*)&uube[ir][it],sizeof(uube[ir][it]));
           swq13.write((char*)&udfe[ir][it],sizeof(udfe[ir][it]));
         }

       cout<<is+1<<" shot source deghosting and redatumming done..."<<endl;

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

        cout<<is+1<<" shot done!"<<endl;
    }
  
        return 0; 

}
































