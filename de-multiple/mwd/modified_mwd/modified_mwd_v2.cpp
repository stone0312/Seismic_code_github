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
   int nx,ns,nr,lt,ltt,nsr,zmax,zmin;
   float dt,dsx,drx,fmin,fmax,dx,dz,v;
   int itmin,itmax;

   int is,ir,it,iz,is1,ir1,is2,ir2;

    ifstream swq;
    swq.open("modified_mwd_v2.par");
    if(!swq)
      {
          cout<<"cannot open modified_mwd_v2.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>fn5>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v;
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

    zmin=9999999;
    for(is=0;is<nx;is++)
       {
          if(dep[is]<zmin)
              zmin=dep[is];
          else
              zmin=zmin;
       }
    cout<<"No. of grids for minimum depth is===="<<zmin<<endl;

    float *s;
    s=alloc1float(nx);
    for(is=0;is<nx;is++)
      s[is]=0.0;
 
    complex<float> *r1;
    r1=alloc1complex(zmax);
    for(is=0;is<zmax;is++)
      r1[is]=(0.0,0.0);

    if(zmin==zmax)
      {
         r1[zmax-1].real()=1.0;
         r1[zmax-1].imag()=0.0;
      }
    else
      {
         for(iz=zmin;iz<zmax;iz++)
           {
             r1[iz].real()=1.0;
             r1[iz].imag()=0.0;
           }
      }

    complex<float> *s;
    s=alloc1complex(nx);

    complex<float> *r;
    r=alloc1complex(nx);

    complex<float> *s_k;
    s_k=alloc1complex(nx);

    complex<float> *r_k;
    r_k=alloc1complex(nx);
 
    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nx);

    complex <float> **d;
    d=alloc2complex(nr,ns);

    complex <float> **mf;
    mf=alloc2complex(ltt,nsr);
  
    float **m;
    m=alloc2float(ltt,nsr);

    complex <float> *gs_k;
    gs_k=alloc1complex(ns);  
    complex <float> *gr_k;
    gr_k=alloc1complex(ns);

    complex <float> *gs;
    gs=alloc1complex(ns);
    complex <float> *gr;
    gr=alloc1complex(ns);

    complex <float> *g;
    g=alloc1complex(ns);

    complex <float> **gs_for;
    gs_for=alloc2complex(zmax,nx);
    complex <float> **gr_for;
    gr_for=alloc2complex(zmax,nx);

    complex <float> **gs_for1;
    gs_for1=alloc2complex(zmax,nx);
    complex <float> **gr_for1;
    gr_for1=alloc2complex(zmax,nx);
 
    complex <float> **gs_back;
    gs_back=alloc2complex(zmax,nx); 
    complex <float> **gr_back;
    gr_back=alloc2complex(zmax,nx);

    complex <float> *A;
    A=alloc1complex(nx);

    complex <float> *B;
    B=alloc1complex(nx);

    complex <float> gtmp1;   
    complex <float> gtmp2;

    complex <float> *gsd;
    gsd=alloc1complex(nx);

    gtmp1=(0.0,0.0);
    gtmp2=(0.0,0.0);

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
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in5,*out5,*in6,*out6;
    fftwf_plan p1,p2,p3,p4,p5,p6;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nx);
    in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
 
    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(ns,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(nx,in4,out4,FFTW_FORWARD,FFTW_MEASURE);
    p5=fftwf_plan_dft_1d(nx,in5,out5,FFTW_FORWARD,FFTW_MEASURE);
    p6=fftwf_plan_dft_1d(ns,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);

    float mem;
    mem=0.0;
    mem+=(nsr*ltt*3+ltt+nx+ns*nr*2+nsr*lt+nx*12+zmax*2)*4.0/1024.0/1024.0;
    cout<<"Memory needed to be allocated is==== "<<mem<<"MB"<<endl;

//zeroing all arrays;
    for(is=0;is<nsr;is++)
      for(it=0;it<ltt;it++)
        {  
          mf[is][it]=(0.0,0.0);
          m[is][it]=0.0;
        }

   for(is=0;is<ns;is++)
       {
          gs_k[is]=(0.0,0.0); 
          gr_k[is]=(0.0,0.0); 
       } 

   for(is=0;is<nx;is++)
    for(iz=0;iz<zmax;iz++)
      {
          gs_for[is][iz]=(0.0,0.0);
          gs_back[is][iz]=(0.0,0.0);
          gr_for[is][iz]=(0.0,0.0);
          gr_back[is][iz]=(0.0,0.0);
          gs_for1[is][iz]=(0.0,0.0);
          gr_for1[is][iz]=(0.0,0.0);
      } 

//calculate the omega and kx
   for(it=0;it<lt;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=lt;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-lt)*1000/(dt*ltt));

   itmin=int(fmin*dt*ltt/1000);
   itmax=int(fmax*dt*ltt/1000)+1;

   cout<<"totally "<<2*(itmax-itmin)<<" frequency slices need to be calculated"<<endl;
   
   for(is=0;is<nx/2;is++)
      kx[is]=2*pai*float(is)/float(dx*nx); 
   for(is=nx/2;is<nx;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nx/2)/float(dx*nx);
 
//read the original data
    ifstream swq1;
    swq1.open(fn4,ios::binary);
    if(!swq1)
       {
          cout<<"cannot open "<<fn4<<endl;
          abort();
       }
    ifstream swq11;
    swq11.open(fn5,ios::binary);
    if(!swq11)
       {
          cout<<"cannot open "<<fn5<<endl;
          abort();
       }

    for(it=0;it<ltt;it++)    
        {
          if(it<itmin)
            {
               for(is=0;is<nsr;is++)
                 mf[is][it]=(0.0,0.0);       
            }      
          else if(it>itmax&&it<ltt-itmax)
            {
               for(is=0;is<nsr;is++)
                 mf[is][it]=(0.0,0.0);
            } 
          else if(it>ltt-itmin&&it<ltt)
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
                {
                   for(ir=0;ir<nr;ir++)
                    {
                      swq1.seekg(((is*nr+ir)*ltt+it)*4,ios::beg);  
                      swq11.seekg(((is*nr+ir)*ltt+it)*4,ios::beg);  
                      swq1.read((char*)&(d[is][ir].real()),sizeof(d[is][ir].real()));
                      swq11.read((char*)&(d[is][ir].imag()),sizeof(d[is][ir].imag()));    
                    }                
                }

//calculating the Green's functions in kx-w domain;
          for(is2=0;is2<nx;is2++)
            {
              for(is=0;is<nx;is++)
                 s[is]=(0.0,0.0);
              s[is2].real()=1.0;
              s[is2].imag()=0.0;
              
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
              
              for(ir2=0;ir2<nx;ir2++)
               { 
                for(is=0;is<nx;is++)
                  r[is]=(0.0,0.0); 
                r[ir2].real()=1.0;
                r[ir2].imag()=0.0;  
             
                if((in5==NULL)||(out5==NULL))
                  cout<<"memory insufficient"<<endl;
                else
                  {
                    for(ir=0;ir<nx;ir++)
                      {
                        in5[ir][0]=r[ir].real();
                        in5[ir][1]=r[ir].imag();
                      }
                  }

                fftwf_execute(p5);

                for(ir=0;ir<nx;ir++)
                   {
                      r_k[ir].real()=out5[ir][0];
                      r_k[ir].imag()=out5[ir][1];
                   }

                for(is=0;is<nx;is++)//wavenumber loop 
                 {
                   if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                     {
                       gs_k[is].real()=0.0;
                       gs_k[is].imag()=0.0;
                       gr_k[is].real()=0.0;
                       gr_k[is].imag()=0.0;
                     }
                   else
                    {
                     gtmp1.real()=0.0;//phase-shift extrapolating operator
                     gtmp1.imag()=dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
//forward extrapolation
                     gs_for[is][0].real()=s_k[is].real();
                     gs_for[is][0].imag()=s_k[is].imag();
                     gr_for[is][0].real()=r_k[is].real();
                     gr_for[is][0].imag()=r_k[is].imag();

                     for(iz=1;iz<zmax;iz++)
                       {
                         gs_for[is][iz]=gs_for[is][iz-1]*exp(gtmp1);
                         gr_for[is][iz]=gr_for[is][iz-1]*exp(gtmp1);
                       }
//backward extrapolaion
                     for(iz=0;iz<zmax;iz++)
                       {
                         gs_for1[is][iz]=gs_for[is][iz-1]*r[iz];
                         gr_for1[is][iz]=gr_for[is][iz-1]*r[iz];
                       }

                     gs_back[is][0]=gs_for1[is][zmax-1];
                     gr_back[is][0]=gr_for1[is][zmax-1];

                     for(iz=1;iz<zmax-1;iz++)   
                      {
                         gs_back[is][iz]=gs_back[is][iz-1]*exp(gtmp1)+gs_for1[is][zmax-1-iz];     
                         gr_back[is][iz]=gr_back[is][iz-1]*exp(gtmp1)+gr_for1[is][zmax-1-iz];
                      }
                     gs_back[is][zmax-1]=gs_back[is][zmax-2]*exp(gtmp1);
                     gr_back[is][zmax-1]=gr_back[is][zmax-2]*exp(gtmp1);

//final Green's functions in kx-w domain for a certain (kx,w) pair;
                     gs_k[is]=gs_back[is][zmax-1];
                     gr_k[is]=gr_back[is][zmax-1];
                    }
                 }

//ifft to transform the Green's function from kx domain to x domain 
              if((in3==NULL)||(out3==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                    for(is=0;is<nx;is++)
                      {
                         in3[is][0]=gs_k[is].real();
                         in3[is][1]=gs_k[is].imag();
                      }
                 }
              fftwf_execute(p3);

              for(is=0;is<nx;is++)
                 {
                    gs[is].real()=out3[is][0];
                    gs[is].imag()=out3[is][1];
                 }

              if((in6==NULL)||(out6==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                    for(is=0;is<nx;is++)
                      {
                         in6[is][0]=gr_k[is].real();
                         in6[is][1]=gr_k[is].imag();
                      }
                 }
              fftwf_execute(p6);

              for(is=0;is<nx;is++)
                 {
                    gr[is].real()=out6[is][0];
                    gr[is].imag()=out6[is][1];
                 }

//for a certain frequency, the Green's function in (x,w) domain calculation done.
                    for(is1=0;is1<nx;is1++)
                       A[is1]=(0.0,0.0);
                    A[is2].real()=1.0;
                    A[is2].imag()=0.0;

                    for(is1=0;is1<nx;is1++)                    
                      gs[is1]=gs[is1]+A[is1]; 

                    for(ir=0;ir<nx;ir++)
                       {
                          for(is1=0;is1<nx;is1++)
                             gsd[is1]=(0.0,0.0);
                          gtmp2=(0.0,0.0);

                          for(ir1=0;ir1<nx;ir1++)
                            B[ir1]=(0.0,0.0);
                          B[ir2].real()=1.0;
                          B[ir2].imag()=0.0;      

                          for(ir1=0;ir1<nx;ir1++)
                             gr[ir1]=gr[ir1]+B[ir1];
                       
                          for(ir1=0;ir1<nx;ir1++)
                             {
                                for(is1=0;is1<nx;is1++)
                                   gsd[ir1]+=gs[is1]*d[is1][ir1];
                             } 
                           for(ir1=0;ir1<nx;ir1++)
                             gtmp2+=gsd[ir1]*gr[ir1];
                           mf[is2*ns+ir2][it]=gtmp2-d[is2][ir2];
                           gtmp2=(0.0,0.0);
                       } 
                   }//end of receiver loop for a cetain frequency 
                 }//end of shot loop for a cetain frequency         
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
        free1float(omega);
        free1float(kx);
        free2complex(d);
        free2complex(mf);
        free2float(m);
        free1complex(gs_k);
        free1complex(gr_k); 
        free1complex(g);
        free2complex(gs_for);
        free2complex(gr_for);
        free2complex(gs_back);
        free2complex(gr_back);
        free1complex(gs);
        free1complex(gr);
        free1complex(A);
        free1complex(B);
        free1complex(gsd);  

        cout<<"all done!"<<endl;

        return 0;
 
   
}






