#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265
#define fm 40.0
#define t0 0.20


int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256];
   int nx,ns,nr,lt,ltt,nsr,zmax,tmp;
   float dt,dsx,drx,fmin,fmax,dx,dz,v,dip_max;
   int itmin,itmax,bou;

   int is,ir,it,iz,is1,ir1;
   int wid,len;

    ifstream swq;
    swq.open("green_for_mwd.par");
    if(!swq)
      {
          cout<<"cannot open green_for_mwd.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou>>dip_max>>wid>>len;
    swq.close();

    cout<<"fna of original data is==== "<<fn1<<endl;
    cout<<"fna of predicted model is==== "<<fn2<<endl;
    cout<<"fna of water layer depth is==== "<<fn3<<endl;
    cout<<"fna of real part of fft of original data is==== "<<fn4<<endl;
    cout<<"fna of imaginary part of fft of original data is==== "<<fn5<<endl;
    cout<<"fna of imaginary part of backward Green's Function in kx-f domain is===="<<fn6<<endl;
    cout<<"fna of imaginary part of backward Green's Function in x-f domain is===="<<fn7<<endl;
    cout<<"No. of vel model width is==== "<<nx<<endl;
    cout<<"lateral interval of vel model is==== "<<dx<<endl;
    cout<<"vertical interval of vel model is==== "<<dz<<endl;
    cout<<"No. of shots is==== "<<ns<<endl;
    cout<<"No. of traces per shot is==== "<<nr<<endl;
    cout<<"No. of samples is==== "<<lt<<endl;
    cout<<"temperal interal (ms) is==== "<<dt<<endl;
    cout<<"spatial interval of sourses is==== "<<dsx<<endl;
    cout<<"spatial interval of receivers is==== "<<dsx<<endl;
    cout<<"minimum frequency is==== "<<fmin<<endl;
    cout<<"maximum frequency is==== "<<fmax<<endl;
    cout<<"velocity of water is==== "<<v<<endl;
    cout<<"maximum dip is===="<<dip_max<<endl;
    cout<<"width of taper for absorbing  is===="<<wid<<endl;
    cout<<"length of taper for frequency is===="<<len<<endl;

    ltt=2*lt; 
    nsr=ns*nr;

    nx=nx+2*bou;


//ordering the depth of water layer to find the maximum depth
    int *dep;
    dep=alloc1int(nx-2*bou);

    ifstream swq3;
    swq3.open(fn3,ios::binary);
    if(!swq3)
       {
          cout<<"cannot open "<<fn3<<endl;
          abort();
       } 
    for(is=0;is<nx-2*bou;is++)
       swq3.read((char*)&dep[is],sizeof(dep[is]));
    swq3.close();

    zmax=0;
    for(is=0;is<nx-2*bou;is++)
       {
          if(dep[is]>zmax)
              zmax=dep[is];
          else
              zmax=zmax;
       }
    cout<<"No. of grids for maximum depth is===="<<zmax<<endl;
    zmax=zmax+1;

    complex <float> *r;
    r=alloc1complex(zmax);
  
//reflectivity for horizonial layer
    for(iz=0;iz<zmax;iz++)
      r[iz]=(0.0,0.0);
    r[zmax-1].real()=1.0;
    r[zmax-1].imag()=0.0; 
 
    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nx);

    complex <float> *g_k;
    g_k=alloc1complex(ns);  

    complex <float> *g;
    g=alloc1complex(ns);

    complex <float> **g_for;
    g_for=alloc2complex(zmax,nx);

    complex <float> **g_for1;
    g_for1=alloc2complex(zmax,nx);

    complex <float> **g_back;
    g_back=alloc2complex(zmax,nx);

    complex <float> *gs;
    gs=alloc1complex(nx);

    complex <float> *gr;
    gr=alloc1complex(nx);
 
    complex <float> **g_all;
    g_all=alloc2complex(ltt,nx);

    complex <float> **g_all_xf;
    g_all_xf=alloc2complex(ltt,nx);

    float **g_final;
    g_final=alloc2float(ltt,nx);

    complex<float> **gz_x;
    gz_x=alloc2complex(zmax,nx);
    
    complex<float> **gback_z_x;
    gback_z_x=alloc2complex(zmax,nx);
 
    complex <float> gtmp1;
    gtmp1.real()=0.0;
    gtmp1.imag()=0.0;

   float * rwt1;
   rwt1=alloc1float(lt);
   float * rwt;
   rwt=alloc1float(ltt);
   complex<float> *rwf;
   rwf=alloc1complex(ltt);
   for(it=0;it<lt;it++)
     rwt1[it]=(1-2*pai*pai*fm*fm*(it*dt/1000-t0)*(it*dt/1000-t0))*exp(-pai*pai*fm*fm*(it*dt/1000-t0)*(it*dt/1000-t0));
   for(it=0;it<ltt;it++)
     rwt[it]=0.0;
   for(it=0;it<lt;it++)
     rwt[it]=rwt1[it];

/*
   for(it=0;it<lt;it++)
     cout<<it<<"  "<<rwt1[it]<<endl;
   return 0;
*/

    complex<float> *s;
    s=alloc1complex(nx);

    complex<float> *uxf;
    uxf=alloc1complex(nx);
    complex<float> *uxf1;
    uxf1=alloc1complex(nx);

    float *abso;
    abso=alloc1float(nx);

    for(is=0;is<nx;is++)
       abso[is]=1.0;

    for(is=0;is<wid;is++)
       abso[is]=sqrt(sin(pai/2*is/(wid-1)));

    for(is=nx-wid;is<nx;is++)
        abso[is]=abso[nx-is];    
      
/*
    for(is=0;is<nx;is++)
      cout<<is<<"   "<<abso[is]<<endl;
    return 0;
*/

    complex<float> *s_k;
    s_k=alloc1complex(nx);

//define the plans and arrays when fftw is complemented
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in5,*out5,*in6,*out6,*in7,*out7,*in8,*out8;
    fftwf_plan  p1,p2,p3,p4,p5,p6,p7,p8;
    
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in7=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out7=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in8=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out8=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);

    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(nx,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(nx,in4,out4,FFTW_FORWARD,FFTW_MEASURE);
    p5=fftwf_plan_dft_1d(nx,in5,out5,FFTW_BACKWARD,FFTW_MEASURE);
    p6=fftwf_plan_dft_1d(nx,in6,out6,FFTW_FORWARD,FFTW_MEASURE);
    p7=fftwf_plan_dft_1d(nx,in7,out7,FFTW_BACKWARD,FFTW_MEASURE);
    p8=fftwf_plan_dft_1d(nx,in8,out8,FFTW_FORWARD,FFTW_MEASURE);

    
    if((in1==NULL)||(out1==NULL))
       cout<<"memory insufficient"<<endl;
    else
       {
          for(it=0;it<ltt;it++)
             {
                in1[it][0]=rwt[it];
                in1[it][1]=0.0;
              }
       }

     fftwf_execute(p1);

     for(it=0;it<ltt;it++)
       {
          rwf[it].real()=out1[it][0]/ltt;
          rwf[it].imag()=out1[it][1]/ltt;
       }
/*  
     for(it=0;it<ltt;it++)
        cout<<it<<"  "<<rwf[it]<<endl;
     return 0;
*/ 
//zeroing all arrays;
   for(is=0;is<ns;is++)
       {
          g_k[is]=(0.0,0.0); 
          g[is]=(0.0,0.0);
       } 
   for(is=0;is<nx;is++)
     for(iz=0;iz<zmax;iz++)
        {
            g_for[is][iz]=(0.0,0.0);
            g_back[is][iz]=(0.0,0.0);
        } 

//calculate the omega and kx
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=-2*pai*1000/(2*dt)+2*pai*(it-ltt/2)*1000/(dt*ltt);

   itmin=int(fmin*dt*ltt/1000);
   itmax=int(fmax*dt*ltt/1000)+1;

   cout<<"totally "<<2*(itmax-itmin)<<" frequency slices need to be calculated"<<endl;

   float *ham;
   ham=alloc1float(itmax);

   for(it=0;it<itmax;it++)
      ham[it]=1.0;  

   for(it=itmin;it<(itmin+len);it++)
      ham[it]=sqrt(sin(pai*(it-itmin+1)/(2*(len))));
   for(it=itmax-len;it<itmax;it++)
      ham[it]=ham[itmin+itmax-it-1];

/*
   for(it=0;it<itmax;it++)
      cout<<it<<"  "<<ham[it]<<endl;
   return 0;
*/

   for(is=0;is<nx/2+1;is++)
      kx[is]=2*pai*float(is)/float(dx*nx); 
   for(is=nx/2+1;is<nx;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nx/2)/float(dx*nx);

/*
    for(is=0;is<nx;is++)
       cout<<is<<"   "<<kx[is]<<endl;
    return 0;
*/

//calculating the Green's functions in kx-w domain;

    cout<<"Calculation of Green's Functions for the positive frequency starts===="<<endl;

    for(it=0;it<ltt/2+1;it++)    
      {
       if((it+1)%100==0)
         cout<<it+1<<"th frequency slices done!"<<endl;
 
       if(it<itmin)     
         {
           for(is=0;is<nx;is++)
              g_all[is][it]=(0.0,0.0);
         }        
 
       else if(it>itmax-1&&it<ltt/2+1)
         {
           for(is=0;is<nx;is++)
              g_all[is][it]=(0.0,0.0);
         }
       else
        { 
          for(is=0;is<nx;is++)
             s[is]=(0.0,0.0);

//          gtmp1.real()=0.0;
//          gtmp1.imag()=omega[it]*t0;
//          s[200+bou]=rwf[it]*exp(gtmp1);

          s[200+bou].real()=1.0*ham[it];
          s[200+bou].imag()=0.0;
/*
          for(is=0;is<nx;is++)
             cout<<is<<"  "<<s[is]<<endl;
          return 0;
*/
          if((in4==NULL)||(out4==NULL))
             cout<<"memory insufficient"<<endl;
          else
            {
               for(is=0;is<nx;is++)
                 {
                    in4[is][0]=s[is].real();
                    in4[is][1]=s[is].imag();
                 }
            }

          fftwf_execute(p4);

          for(is=0;is<nx;is++)
            {
               s_k[is].real()=out4[is][0]/sqrt(nx);
               s_k[is].imag()=out4[is][1]/sqrt(nx);
            }
  
/* 
          cout<<"frequency slice is===="<<it<<endl; 
            for(is=0;is<nx;is++)
              cout<<is<<"   "<<s_k[is]<<endl;
            return 0;
*/         

          for(is=0;is<nx;is++)
            for(iz=0;iz<zmax;iz++) 
                g_for[is][iz]=(0.0,0.0);

          for(is=0;is<nx;is++)
             g_for[is][0]=s_k[is];

//forward extrapolation
          for(iz=1;iz<zmax;iz++)
            {   
              for(is=0;is<nx/2+1;is++)
                 {
                   if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                     {
                          g_for[is][iz].real()=0.0;
                          g_for[is][iz].imag()=0.0;
                     }
                
//                   else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2)||pow(omega[it]/v,2)-pow(kx[is],2)==pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                   else
                    {
                      gtmp1.real()=0.0;               //phase-shift extrapolating operator
                      gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                      g_for[is][iz]=g_for[is][iz-1]*exp(gtmp1);
                    }
/*
                   else
                    {
                       g_for[is][iz].real()=0.0;
                       g_for[is][iz].imag()=0.0;
                    }
*/
                 }
             for(is=nx/2+1;is<nx;is++)
                {
                   g_for[is][iz].real()=g_for[nx-is][iz].real();
                   g_for[is][iz].imag()=-g_for[nx-is][iz].imag();
                }

              if((in7==NULL)||(out7==NULL))
                cout<<"memory insufficient"<<endl;
              else
               {
                  for(is=0;is<nx;is++)
                    {
                       in7[is][0]=g_for[is][iz].real();
                       in7[is][1]=g_for[is][iz].imag();
                    }
               }

              fftwf_execute(p7);

              for(is=0;is<nx;is++)
                {
                  uxf[is].real()=out7[is][0]/sqrt(nx);
                  uxf[is].imag()=out7[is][1]/sqrt(nx);
                }

              for(is=0;is<nx;is++)
                uxf1[is]=uxf[is]*abso[is];

/*
              for(is=0;is<nx;is++)
                cout<<is<<"  "<<uxf[is]<<endl;
              return 0;
*/

              if((in8==NULL)||(out8==NULL))
                cout<<"memory insufficient"<<endl;
              else
               {
                  for(is=0;is<nx;is++)
                    {
                       in8[is][0]=uxf1[is].real();
                       in8[is][1]=uxf1[is].imag();
                    }
               }

              fftwf_execute(p8);

              for(is=0;is<nx;is++)
                {
                  g_for[is][iz].real()=out8[is][0]/sqrt(nx);
                  g_for[is][iz].imag()=out8[is][1]/sqrt(nx);
                }
           }//end of iz for all depth  forward extraplotation
           
//backward extrapolaion
         for(is=0;is<nx;is++)
           for(iz=0;iz<zmax;iz++)
             g_for1[is][iz]=g_for[is][iz]*r[iz];

         for(is=0;is<nx;is++)
            g_back[is][0]=g_for1[is][zmax-1];   
         
        for(iz=1;iz<zmax-1;iz++)
         { 
            for(is=0;is<nx/2+1;is++)
               {
                 if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                   {
                       g_back[is][iz].real()=0.0;
                       g_back[is][iz].imag()=0.0;
                   }
               
//                 else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2)||pow(omega[it]/v,2)-pow(kx[is],2)==pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))

                 else
                   {
                      gtmp1.real()=0.0;               //phase-shift extrapolating operator
                      gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                      g_back[is][iz]=g_back[is][iz-1]*exp(gtmp1)+g_for1[is][zmax-1-iz];  

                   }
/*
                 else
                    {
                      g_back[is][iz].real()=0.0;
                      g_back[is][iz].imag()=0.0;
                    }
*/
                }//end of is

            for(is=nx/2+1;is<nx;is++)
               {
                  g_back[is][iz].real()=g_back[nx-is][iz].real();
                  g_back[is][iz].imag()=-g_back[nx-is][iz].imag();
               }
              
              if((in5==NULL)||(out5==NULL))
                cout<<"memory insufficient"<<endl;
              else
               {
                  for(is=0;is<nx;is++)
                    {
                       in5[is][0]=g_back[is][iz].real();
                       in5[is][1]=g_back[is][iz].imag();
                    }
               }

              fftwf_execute(p5);

              for(is=0;is<nx;is++)
                {
                  uxf[is].real()=out5[is][0]/sqrt(nx);
                  uxf[is].imag()=out5[is][1]/sqrt(nx);
                }

              for(is=0;is<nx;is++)
                uxf[is]=uxf[is]*abso[is];

              if((in6==NULL)||(out6==NULL))
                cout<<"memory insufficient"<<endl;
              else
               {
                  for(is=0;is<nx;is++)
                    {
                       in6[is][0]=uxf[is].real();
                       in6[is][1]=uxf[is].imag();
                    }
               }

              fftwf_execute(p6);

              for(is=0;is<nx;is++)
                {
                  g_back[is][iz].real()=out6[is][0]/sqrt(nx);
                  g_back[is][iz].imag()=out6[is][1]/sqrt(nx);
                }
          }//end of iz for back extrapolation
 
         for(is=0;is<nx;is++)
             g_back[is][zmax-1]=g_back[is][zmax-2]*exp(gtmp1);
          
//final Green's functions in kx-w domain for a certain (kx,w) pair;
         for(is=0;is<nx;is++)
            g_all[is][it]=g_back[is][zmax-1];

       } //end of else for frequency from itmin to itmax                                                                   

      }//end of all positive frequency
 
    cout<<"Calculation of Green's Function for Positive Frequency in Kx domain done!"<<endl;      
    cout<<"Calculation of Green's Function for Necgtive Frequency in Kx domain starts:"<<endl;              

    for(is=0;is<nx;is++)
      {
         for(it=ltt/2+1;it<ltt;it++)
           {
              g_all[is][it].real()=g_all[is][ltt-it].real();
              g_all[is][it].imag()=-g_all[is][ltt-it].imag();
           }
       }

    cout<<"Calculation of Green's Function for All Frequency in Kx domain done!"<<endl;


    ofstream swq8;
    swq8.open("/data2/swq/kx_omega_domain_green_function_real.dat",ios::binary);

    ofstream swq9;
    swq9.open("/data2/swq/kx_omega_domain_green_function_imag.dat",ios::binary);

    for(is=0;is<nx;is++)
        for(it=0;it<ltt;it++)
           {
               swq8.write((char*)&g_all[is][it].real(),sizeof(g_all[is][it].real()));
               swq9.write((char*)&g_all[is][it].imag(),sizeof(g_all[is][it].imag()));
           }
    swq8.close();
    swq9.close();


//ifft to transform the Green's function from kx domain to x domain 
         for(it=0;it<ltt;it++)
            {      
              if((in3==NULL)||(out3==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                    for(is=0;is<nx;is++)
                      {
                         in3[is][0]=g_all[is][it].real();
                         in3[is][1]=g_all[is][it].imag();
                      }
                 }

              fftwf_execute(p3);

              for(is=0;is<nx;is++)
                 {
                    g_all_xf[is][it].real()=out3[is][0]/sqrt(nx);
                    g_all_xf[is][it].imag()=out3[is][1]/sqrt(nx);
                 }
           }

         ofstream swq6;
         swq6.open(fn6,ios::binary);
         if(!swq6)
            {
               cout<<"cannot open "<<fn6<<endl;
               abort();
            }
          for(is=0;is<nx;is++)
            for(it=0;it<ltt;it++)
               swq6.write((char*)&g_all_xf[is][it].real(),sizeof(g_all_xf[is][it].real()));
          swq6.close();
 
         ofstream swq7;
         swq7.open(fn7,ios::binary);
         if(!swq7)
            {
               cout<<"cannot open "<<fn7<<endl;
               abort();
            }

          for(is=0;is<nx;is++)
            for(it=0;it<ltt;it++)
               swq7.write((char*)&g_all_xf[is][it].imag(),sizeof(g_all_xf[is][it].imag()));

          swq7.close();

          for(is=0;is<nx;is++)
            {
               for(it=0;it<ltt;it++)
                  {          
                    in2[it][0]=g_all_xf[is][it].real();
                    in2[it][1]=g_all_xf[is][it].imag();
                  }
           
               fftwf_execute(p2);

               for(it=0;it<ltt;it++)
                  g_final[is][it]=out2[it][0]/sqrt(ltt);
            }
       
        cout<<"wirting disk starts... "<<endl;
 
        ofstream swq2;   
        swq2.open(fn2,ios::binary);
        if(!swq2)
          {
             cout<<"cannot open "<<fn2<<endl;
             abort();
          }
 
        for(is=0;is<nx;is++)
           for(it=0;it<ltt;it++)        
             swq2.write((char*)&g_final[is][it],sizeof(g_final[is][it]));  
        swq2.close();

        cout<<"disk written."<<endl;

//free the memory
        free1float(omega);
        free1float(kx);
        free1complex(g_k);
        free1complex(g);
        free2complex(g_for);
        free2complex(g_back);
        free1complex(gs);
        free1complex(gr);

        cout<<"all done!"<<endl;

        return 0;
   
}








