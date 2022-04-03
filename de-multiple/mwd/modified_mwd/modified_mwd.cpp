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
   char fn1[256],fn2[256],fn3[256];
   int nx,ns,nr,lt,ltt,nsr,zmax;
   float dt,dsx,drx,fmin,fmax,dx,dz,v;
   int itmin,itmax;

   int is,ir,it,iz,is1,ir1;

    ifstream swq;
    swq.open("modified_mwd.par");
    if(!swq)
      {
          cout<<"cannot open modified_mwd.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v;
    swq.close();

    cout<<"fna of original data is==== "<<fn1<<endl;
    cout<<"fna of predicted model is==== "<<fn2<<endl;
    cout<<"fna of water layer depth is==== "<<fn3<<endl;
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
 
    float **u;
    u=alloc2float(ltt,nsr);

    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nx);

    complex <float> **uf;
    uf=alloc2complex(ltt,nsr);

    complex <float> **d;
    d=alloc2complex(nr,ns);

    complex <float> **mf;
    mf=alloc2complex(ltt,nsr);
  
    float **m;
    m=alloc2float(ltt,nsr);

    complex <float> *g_k;
    g_k=alloc1complex(ns);  

    complex <float> *g;
    g=alloc1complex(ns);

    complex <float> *g_for;
    g_for=alloc1complex(zmax);

    complex <float> *g_back;
    g_back=alloc1complex(zmax);

    complex <float> *gs;
    gs=alloc1complex(nx);

    complex <float> *gr;
    gr=alloc1complex(nx);
 
    complex <float> *A;
    A=alloc1complex(nx);

    complex <float> *B;
    B=alloc1complex(nx);

    complex <float> gtmp1;
    complex <float> gtmp2;
    complex <float> gtmp3;
    
    complex <float> *gsd;
    gsd=alloc1complex(nx);

    gtmp1=(0.0,0.0);
    gtmp2=(0.0,0.0);
    gtmp3=(0.0,0.0);

    for(is=0;is<nx;is++)
      gsd[is]=(0.0,0.0);

/*
    complex <float> **gf_all;
    gf_all=alloc2complex(ltt,ns);
       
    float **gt1;
    gt1=alloc2float(ltt,ns);

    float **gt;
    gt=alloc2float(lt,ns);
*/

//define the plans and arrays when fftw is complemented
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3;
    fftwf_plan p1, p2, p3;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);

    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(ns,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);

    float mem;
    mem=0.0;
    mem+=(nsr*ltt*4+ltt+nx+ns*nr*2+nsr*lt+nx*12+zmax*2)*4.0/1024.0/1024.0;
    cout<<"Memory needed to be allocated is==== "<<mem<<"MB"<<endl;

//zeroing all arrays;
    for(is=0;is<nsr;is++)
      for(it=0;it<ltt;it++)
        {  
          u[is][it]=0.0;
          uf[is][it]=(0.0,0.0);
          mf[is][it]=(0.0,0.0);
          m[is][it]=0.0;
        }

   for(is=0;is<ns;is++)
       {
          g_k[is]=(0.0,0.0); 
          g[is]=(0.0,0.0);
       } 

   for(iz=0;iz<zmax;iz++)
      {
          g_for[iz]=(0.0,0.0);
          g_back[iz]=(0.0,0.0);
      } 

//calculate the omega and kx
   for(it=0;it<lt;it++)
      omega[it]=2*pai*it*1000/(2*dt*ltt);
   for(it=lt;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-lt)*1000/(2*dt*ltt));

   itmin=int(fmin*2*dt*ltt/1000);
   itmax=int(fmax*2*dt*ltt/1000)+1;

   cout<<"totally "<<2*(itmax-itmin)<<" frequency slices need to be calculated"<<endl;
   
   for(is=0;is<nx/2;is++)
      kx[is]=2*pai*is/(2*dx*nx); 
   for(is=nx/2;is<nx;is++)
      kx[is]=2*pai*(-1/(2*dx)+(is-nx/2)/(2*dx*nx));

//read the original data
    ifstream swq1;
    swq1.open(fn1,ios::binary);
    if(!swq1)
       {
          cout<<"cannot open "<<fn1<<endl;
          abort();
       }
    for(is=0;is<nsr;is++)
      for(it=0;it<lt;it++)
        swq1.read((char*)&u[is][it],sizeof(u[is][it]));
    swq1.close();

//fft of original data     
    for(is=0;is<nsr;is++)
      {
         if((in1==NULL)||(out1==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
              for(it=0;it<ltt;it++)
                 {
                    in1[it][0]=u[is][it];
                    in1[it][1]=0.0;
                 }
           }

        fftwf_execute(p1);

        for(it=0;it<ltt;it++)
          {
             uf[is][it].real()=out1[it][0];
             uf[is][it].imag()=out1[it][1];
          }
      }
 
    cout<<"FFT of original data finished!"<<endl;

    for(it=0;it<ltt;it++)    
        {
          if(it<itmin)
            {
               for(is=0;is<nsr;is++)
                 mf[is][it]=(0.0,0.0);       
            }      
          else if(it>itmax&&it<ltt-itmax+itmin)
            {
               for(is=0;is<nsr;is++)
                 mf[is][it]=(0.0,0.0);
            } 
          else
           { 
              if((it+1)%5==0)
                 cout<<it<<"frequency slices have been finished!"<<endl;

//forming the data matrix for a certain frequcncy
              for(is=0;is<ns;is++)
                for(ir=0;ir<nr;ir++)
                   d[is][ir]=uf[is*nr+ir][it];                    
      
//calculating the Green's functions in kx-w domain;
              for(is=0;is<nx;is++)
                 {
                     gtmp1=(0.0,dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2)));//phase-shift extrapolating operator

//forward extrapolation
                     g_for[0]=(1.0,0.0);
                     for(iz=1;iz<zmax;iz++)
                        g_for[iz]=g_for[iz-1]*exp(gtmp1);

//backward extrapolaion
                     g_back[0]=g_for[zmax-1];
                     for(iz=1;iz<zmax-1;iz++)   
                        g_back[iz]=g_back[iz-1]*gtmp1+g_for[zmax-1-iz];     
                     g_back[zmax-1]=g_back[zmax-2]*gtmp1;

//final Green's functions in kx-w domain for a certain (kx,w) pair;
                     g_k[is]=g_back[zmax-1];
              
                     for(iz=0;iz<zmax;iz++)
                        {
                           g_for[iz]=(0.0,0.0);
                           g_back[iz]=(0.0,0.0);
                        }                                                 
                 }
//ifft to transform the Green's function from kx domain to x domain 
              if((in3==NULL)||(out3==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                    for(is=0;is<nx;is++)
                      {
                         in3[is][0]=g_k[is].real();
                         in3[is][1]=g_k[is].imag();
                      }
                 }

              fftwf_execute(p3);

              for(is=0;is<nx;is++)
                 {
                    g[is].real()=out3[is][0];
                    g[is].imag()=out3[is][1];
                 }
//for a certain frequency, the Green's function in (x,w) domain calculation done.
               for(is=0;is<nx;is++)
                 {
                    for(is1=0;is1<nx;is1++)
                       A[is1]=(0.0,0.0);
                    A[is]=(1.0,0.0);
                    for(is1=0;is1<nx;is1++)                    
                      gs[is1]=g[is1]+A[is1]; 

                    for(ir=0;ir<nx;ir++)
                       {
                          for(is1=0;is1<nx;is1++)
                             gsd[is1]=(0.0,0.0);
                          gtmp2=(0.0,0.0);

                          for(ir1=0;ir1<nx;ir1++)
                            B[ir1]=(0.0,0.0);
                          B[ir]=(1.0,0.0);
                          for(ir1=0;ir1<nx;ir1++)
                             gr[ir1]=g[ir1]+B[ir1];
                       
                          for(ir1=0;ir1<nx;ir1++)
                             {
                                for(is1=0;is1<nx;is1++)
                                   gsd[ir1]+=gs[is1]*d[is1][ir1];
                             } 
                           for(ir1=0;ir1<nx;ir1++)
                             gtmp2+=gsd[ir1]*gr[ir1];
                           mf[is*ns+ir][it]=gtmp2-d[is][ir];
                           gtmp2=(0.0,0.0);

                       } 
                 }        
            }//end of else,.ie. a certain frequency
         }//end of all frequency
//ifft the (x,w) domain predicted multiple to (x,t) domain
        for(is=0;is<nsr;is++)
           {
              if((in2==NULL)||(out2==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                   for(it=0;it<ltt;it++)
                     {
                        in2[it][0]=mf[is][it].real();
                        in2[it][1]=mf[is][it].imag();
                     }
                 }

              fftwf_execute(p2);

              for(it=0;it<ltt;it++)
                 m[is][it]=out2[it][0];

              for(it=0;it<ltt;it++)
                  {
                      mf[is][it].real()=0.0;
                      mf[is][it].imag()=0.0;
                  }
            }
       
        cout<<"wirting disk starts... "<<endl;
 
        ofstream swq2;   
        swq2.open(fn2,ios::binary);
        if(!swq2)
          {
             cout<<"cannot open "<<fn2<<endl;
             abort();
          }
 
        for(is=0;is<nsr;is++)
           for(it=0;it<lt;it++)        
               swq2.write((char*)&m[is][it],sizeof(m[is][it]));  
        swq2.close();

        cout<<"disk written."<<endl;

//free the memory
        free2float(u);
        free1float(omega);
        free1float(kx);
        free2complex(uf);
        free2complex(d);
        free2complex(mf);
        free2float(m);
        free1complex(g_k);
        free1complex(g);
        free1complex(g_for);
        free1complex(g_back);
        free1complex(gs);
        free1complex(gr);
        free1complex(A);
        free1complex(B);
        free1complex(gsd);  

        cout<<"all done!"<<endl;

        return 0;
 
   
}




































































