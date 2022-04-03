#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256];
   int nx,ns,nr,lt,ltt,nsr,zmax;
   float dt,dsx,drx,fmin,fmax,dx,dz,v;
   int itmin,itmax;

   int is,ir,it,iz,is1,ir1;
   int wid,wid1;
   int k;

    ifstream swq;
    swq.open("modified_mwd_v3.par");
    if(!swq)
      {
          cout<<"cannot open modified_mwd_v3.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>fn5>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>wid>>wid1;
    swq.close();

    cout<<"fna of original data is==== "<<fn1<<endl;
    cout<<"fna of predicted model is==== "<<fn2<<endl;
    cout<<"fna of water layer depth is==== "<<fn3<<endl;
    cout<<"fna of real part of fft of original data is==== "<<fn4<<endl;
    cout<<"fna of imaginary part of fft of original data is==== "<<fn5<<endl;
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
    cout<<"width of Hamming window is===="<<wid<<endl;
    cout<<"width of Gaussian window is===="<<wid1<<endl;

    ltt=2*lt; 
    nsr=ns*nr;

//ordering the depth of water layer to find the maximum depth
    int *dep;
    dep=alloc1int(nx);

    ifstream swq3;
    swq3.open(fn3,ios::binary);
    if(!swq3)
       {
          cout<<"cannot open "<<fn3<<endl;
          abort();
       } 
    for(is=0;is<nx;is++)
       swq3.read((char*)&dep[is],sizeof(dep[is]));
    swq3.close();

    zmax=0;
    for(is=0;is<nx;is++)
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

    complex<float> *s;
    s=alloc1complex(nx);

    for(is=0;is<nx;is++)
      s[is]=(0.0,0.0);
    s[1000].real()=1.0;
    s[1000].imag()=0.0;

    float *gauss;
    gauss=alloc1float(nx);

    for(is=0;is<nx;is++)
       gauss[is]=1.0;
    
    for(is=nx-wid1-1;is<nx;is++)
       gauss[is]=0.5+0.5*cos(pai*(is-(nx-wid1-1))/(wid1+1));

    for(is=0;is<wid1+1;is++)
       {
          k=nx-is-1;
          gauss[is]=gauss[k];
       }
     
    complex<float> *s_k;
    s_k=alloc1complex(nx);

//define the plans and arrays when fftw is complemented
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in33,*out33;
    fftwf_plan p1, p2, p3,p4,p33;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in33=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out33=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);

    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(ns,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(nx,in4,out4,FFTW_FORWARD,FFTW_MEASURE);
    p33=fftwf_plan_dft_1d(ns,in33,out33,FFTW_BACKWARD,FFTW_MEASURE);

    float mem;
    mem=0.0;
    mem+=(nsr*ltt+ltt+nx+ns*nr*2+nsr*lt+nx*12+zmax*2)*4.0/1024.0/1024.0;
    cout<<"Memory needed to be allocated is==== "<<mem<<"MB"<<endl;


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
   for(it=0;it<lt;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=lt;it<ltt;it++)
      omega[it]=-2*pai*1000/(2*dt)+2*pai*(it-lt)*1000/(dt*ltt);

   itmin=int(fmin*dt*ltt/1000);
   itmax=int(fmax*dt*ltt/1000)+1;

   cout<<"Mininum and Maximum frequency slice is==== "<<itmin<<" , "<<itmax<<endl;
   cout<<"totally "<<2*(itmax-itmin)<<" frequency slices need to be calculated"<<endl;
  
   float *ham;
   ham=alloc1float(ltt);

   for(it=0;it<ltt;it++)
      ham[it]=1.0;  
   for(it=lt-wid;it<lt+1;it++)
      ham[it]=0.5+0.5*cos(pai*(it-(lt-wid))/wid);
   for(it=lt+1;it<lt+wid+2;it++)
     {
       k=ltt+1-it;
       ham[it]=ham[k];
     }

   for(is=0;is<nx/2;is++)
      kx[is]=2*pai*float(is)/float(dx*nx); 
   for(is=nx/2;is<nx;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nx/2)/float(dx*nx);

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
          s_k[is].real()=out4[is][0];
          s_k[is].imag()=out4[is][1];
      }

    for(it=100;it<101;it++)
//    for(it=0;it<ltt;it++)    
        {
//calculating the Green's functions in kx-w domain;
          cout<<it<<" Green's function of frequency slices====="<<endl;

              for(is=0;is<nx;is++)
                  s_k[is]=s_k[is];
              for(is=0;is<nx;is++)
                 {
                   if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                     {
                       for(iz=0;iz<zmax;iz++)
                        {
                          g_for[is][iz].real()=0.0;
                          g_for[is][iz].imag()=0.0;
                          g_back[is][iz].real()=0.0;
                          g_back[is][iz].imag()=0.0;
                       }
                      g_all[is][it]=(0.0,0.0);
                     }

                   else if(pow(omega[it]/v,2)-pow(kx[is],2)==0.0)
                     {
                       for(iz=0;iz<zmax;iz++)
                        {
                          g_for[is][iz].real()=0.0;
                          g_for[is][iz].imag()=0.0;
                          g_back[is][iz].real()=0.0;
                          g_back[is][iz].imag()=0.0;
                       }
                      g_all[is][it]=(0.0,0.0);
                     }

                   else
                    {
                     gtmp1.real()=0.0;               //phase-shift extrapolating operator
                     if(omega[it]>0.0) 
                        gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                     else
                        gtmp1.imag()=dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));

//forward extrapolation
                     g_for[is][0].real()=s_k[is].real();
                     g_for[is][0].imag()=s_k[is].imag();

                     for(iz=1;iz<zmax;iz++)
                       g_for[is][iz]=g_for[is][iz-1]*exp(gtmp1);
//backward extrapolaion
                     for(iz=0;iz<zmax;iz++)
                        g_for1[is][iz]=g_for[is][iz]*r[iz];
  
                     g_back[is][0]=g_for1[is][zmax-1];
                     for(iz=1;iz<zmax-1;iz++) 
                      {   
                        g_back[is][iz]=g_back[is][iz-1]*exp(gtmp1)+g_for1[is][zmax-1-iz];     
                        g_back[is][iz]=g_back[is][iz];
                      }
                     g_back[is][zmax-1]=g_back[is][zmax-2]*exp(gtmp1);
                     g_back[is][zmax-1]=g_back[is][zmax-1];
//final Green's functions in kx-w domain for a certain (kx,w) pair;
                       g_all[is][it]=g_back[is][zmax-1];
                     }                                                 
                 }


     
           for(iz=0;iz<zmax;iz++)
            {      
              if((in3==NULL)||(out3==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                    for(is=0;is<nx;is++)
                      {
                         in3[is][0]=g_for[is][iz].real();
                         in3[is][1]=g_for[is][iz].imag();
                      }
                 }

              fftwf_execute(p3);

              for(is=0;is<nx;is++)
                 {
                    gz_x[is][iz].real()=out3[is][0];
                    gz_x[is][iz].imag()=out3[is][1];
                 }
             }  
            
             ofstream swq22;
             swq22.open("/data1/swq/srme_data/check_forward_Green_real_shot100.dat",ios::binary);
             if(!swq22)
          	{
           	  cout<<"cannot open "<<fn2<<endl;
             	  abort();
          	}

             for(is=0;is<nx;is++)
                for(iz=0;iz<zmax;iz++)
                  swq22.write((char*)&gz_x[is][iz].real(),sizeof(gz_x[is][iz].real()));
             swq22.close();

             for(iz=0;iz<zmax;iz++)
               {
                 if((in33==NULL)||(out33==NULL))
                   cout<<"memory insufficient"<<endl;
                 else
                   {
                     for(is=0;is<nx;is++)
                       {
                         in33[is][0]=g_back[is][iz].real();
                         in33[is][1]=g_back[is][iz].imag();
                       }
                    }

                 fftwf_execute(p33);

                 for(is=0;is<nx;is++)
                   {
                      gback_z_x[is][iz].real()=out33[is][0];
                      gback_z_x[is][iz].imag()=out33[is][1];
                   }
                }

             ofstream swq222;
             swq222.open("/data1/swq/srme_data/check_backward_Green_real_shot100.dat",ios::binary);
             if(!swq222)
                {
                  cout<<"cannot open "<<endl;
                  abort();
                }

             for(is=0;is<nx;is++)
                for(iz=0;iz<zmax;iz++)
                  swq222.write((char*)&gback_z_x[is][iz].real(),sizeof(gback_z_x[is][iz].real()));
             swq222.close();

             cout<<"disk written."<<endl;
 
//             return 0;



         
          }//end of all frequency
                
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
                    g_all_xf[is][it].real()=out3[is][0];
                    g_all_xf[is][it].imag()=out3[is][1];
                 }
           }

         ofstream swq22222;
         swq22222.open("/data2/swq/check_kf_backward_Green_imaginary_shot100_taper100.dat",ios::binary);
         if(!swq22222)
            {
               cout<<"cannot open "<<endl;
               abort();
            }

          for(is=0;is<nx;is++)
            for(it=0;it<ltt;it++)
               swq22222.write((char*)&g_all[is][it].imag(),sizeof(g_all[is][it].imag()));
          swq22222.close();

 
         ofstream swq2222;
         swq2222.open("/data2/swq/check_xf_backward_Green_imaginary_shot100_taper100.dat",ios::binary);
         if(!swq2222)
            {
               cout<<"cannot open "<<endl;
               abort();
            }

          for(is=0;is<nx;is++)
            for(it=0;it<ltt;it++)
               swq2222.write((char*)&g_all_xf[is][it].imag(),sizeof(g_all_xf[is][it].imag()));
          swq2222.close();
 
//ifft the (x,w) domain predicted multiple to (x,t) domain
        for(is=0;is<nx;is++)
           {
              if((in2==NULL)||(out2==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                   for(it=0;it<ltt;it++)
                     {
                         in2[it][0]=g_all_xf[is][it].real();
                         in2[it][1]=g_all_xf[is][it].imag();
                     }
                 }

              fftwf_execute(p2);

              for(it=0;it<ltt;it++)
                 g_final[is][it]=out2[it][0];

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




































































