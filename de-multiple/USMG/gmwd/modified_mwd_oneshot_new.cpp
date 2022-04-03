#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265
#define t0 0.2
#define fm 40.0
int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256];
   int nx,ns,nr,lt,ltt,nsr,zmax,zmin,bou,nxb,wid;
   float dt,dsx,drx,fmin,fmax,dx,dz,v,dip_max;
   int itmin,itmax,shot,len;

   int is,ir,it,iz,is1,ir1,is2,ir2,itmp,ix,tmp,tmp1;

    ifstream swq;
    swq.open("modified_mwd_oneshot_new.par");
    if(!swq)
      {
          cout<<"cannot open modified_mwd_oneshot_new.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou>>shot>>dip_max>>wid>>len;
    swq.close();

    cout<<"fna of original data is==== "<<fn1<<endl;
    cout<<"fna of predicted model is==== "<<fn2<<endl;
    cout<<"fna of water layer depth is==== "<<fn3<<endl;
    cout<<"fna of real part of fft of original data is==== "<<fn4<<endl;
    cout<<"fna of imaginary part of fft of original data is==== "<<fn5<<endl;
    cout<<"fna of real part of fft for positive frequency of predicted model is==== "<<fn6<<endl;
    cout<<"fna of imaginary part of fft for positive frequency of predicted model is==== "<<fn7<<endl;
    cout<<"fna of real part of fft for all frequency of predicted model is==== "<<fn8<<endl;
    cout<<"fna of imaginary part of fft for all frequency of predicted model is==== "<<fn9<<endl;
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
    cout<<"half of width of boundary is==== "<<bou<<endl;
    cout<<"No. of shot to be predicted is===="<<shot<<endl;
    cout<<"maximum dip is===="<<dip_max<<endl;
    cout<<"width of taper for absorbing  is===="<<wid<<endl;
    cout<<"length of taper for gibbs effects of frequency is===="<<len<<endl;

    ltt=2*lt; 
    nsr=ns*nr;
//    nxb=nx+2*bou;
    nxb=nx*(dsx/dx)+2*bou;
/*
    cout<<"nxb==="<<nxb<<endl;
    return 0;
*/
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
/*
    for(is=0;is<nx;is++)
       cout<<is<<"  "<<dep[is]<<endl;
*/
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

    zmin=9999999;
    for(is=0;is<nx;is++)
       {
          if(dep[is]<zmin)
              zmin=dep[is];
          else
              zmin=zmin;
       }
    cout<<"No. of grids for minimum depth is===="<<zmin<<endl;

//    return 0;

    complex<float> *r1;
    r1=alloc1complex(zmax);
    for(is=0;is<zmax;is++)
      r1[is]=(0.0,0.0);

    if(zmin==zmax-1)
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
    s=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
       s[is]=(0.0,0.0);

    complex<float> *r;
    r=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
       r[is]=(0.0,0.0);

    complex<float> *s_k;
    s_k=alloc1complex(nxb);
 
    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nxb);

    complex <float> *dd;
    dd=alloc1complex(nsr);

    complex <float> **d;
    d=alloc2complex(nr,ns);

    complex <float> mf;
  
    float *m;
    m=alloc1float(ltt);

    complex <float> **g_k;
    g_k=alloc2complex(nxb,nx);

    complex <float> **g;
    g=alloc2complex(nxb,nx);

    complex<float> *uxf;
    uxf=alloc1complex(nxb);

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
   for(it=0;it<ltt;it++)
     cout<<it<<"  "<<rwt[it]<<endl;
   return 0;
*/
    float *abso;
    abso=alloc1float(nxb);
    
    for(is=0;is<nxb;is++)
       abso[is]=1.0;
    
    for(is=0;is<wid;is++)
       abso[is]=sqrt(sin(pai/2*is/(wid-1)));
    
    for(is=nxb-wid;is<nxb;is++)
        abso[is]=abso[nxb-is];
      
    complex <float> **gs_for;
    gs_for=alloc2complex(zmax,nxb);

    complex <float> **gs_for1;
    gs_for1=alloc2complex(zmax,nxb);
 
    complex <float> **gs_back;
    gs_back=alloc2complex(zmax,nxb); 

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
    gtmp1=(0.0,0.0);
    gtmp2=(0.0,0.0);
    gtmp3=(0.0,0.0);
    
    complex <float> *gsd;
    gsd=alloc1complex(nx);

    for(is=0;is<nx;is++)
      gsd[is]=(0.0,0.0);

    complex<float> *mf_tmp;
    mf_tmp=alloc1complex(nsr);
    for(is=0;is<nsr;is++)
       mf_tmp[is]=(0.0,0.0);

    complex<float> *f_tmp;
    f_tmp=alloc1complex(ltt);

//define the plans and arrays when fftw is complemented
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in5,*out5,*in6,*out6,*in7,*out7,*in8,*out8;
    fftwf_plan p1,p2,p3,p4,p5,p6,p7,p8;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in7=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out7=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in8=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out8=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
 
    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(nxb,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(nxb,in4,out4,FFTW_FORWARD,FFTW_MEASURE);
    p5=fftwf_plan_dft_1d(nxb,in5,out5,FFTW_BACKWARD,FFTW_MEASURE);
    p6=fftwf_plan_dft_1d(nxb,in6,out6,FFTW_FORWARD,FFTW_MEASURE);
    p7=fftwf_plan_dft_1d(nxb,in7,out7,FFTW_BACKWARD,FFTW_MEASURE);
    p8=fftwf_plan_dft_1d(nxb,in8,out8,FFTW_FORWARD,FFTW_MEASURE);

    float mem;
    mem=0.0;
    mem+=(nx*nxb*2+ltt+ltt+nx*4+ns*nr*2+nxb*6+zmax*nxb*4)*4.0/1024.0/1024.0;
    cout<<"Memory needed to be allocated is==== "<<mem<<"MB"<<endl;

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
          rwf[it].real()=out1[it][0];
          rwf[it].imag()=out1[it][1];
       }
/*
     for(it=0;it<ltt;it++)
       cout<<it<<"  "<<rwf[it]<<endl;
     return 0;
*/
//zeroing all arrays;
      for(it=0;it<ltt;it++)
        m[it]=0.0;

   for(is=0;is<nxb;is++)
    for(iz=0;iz<zmax;iz++)
      {
          gs_for[is][iz]=(0.0,0.0);
          gs_back[is][iz]=(0.0,0.0);
          gs_for1[is][iz]=(0.0,0.0);
      } 

//calculate the omega and kx
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   itmin=int(fmin*dt*ltt/1000);
   itmax=int(fmax*dt*ltt/1000)+1;

   cout<<"totally "<<(itmax-itmin)<<" frequency slices need to be calculated"<<endl;

   float *ham;
   ham=alloc1float(itmax);
   for(it=0;it<itmax;it++)
      ham[it]=1.0;
   for(it=itmin;it<(itmin+len);it++)
      ham[it]=sqrt(sin(pai*(it-itmin+1)/(2*(len))));
   for(it=itmax-len;it<itmax;it++)
      ham[it]=ham[itmin+itmax-it-1];

/*    
    for(ix=0;ix<nxb;ix++)
       cout<<ix<<"===="<<hann1[ix]<<endl;
    return 0;
*/
   for(is=0;is<nxb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dx*nxb); 
   for(is=nxb/2+1;is<nxb;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nxb/2)/float(dx*nxb);
 
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

    ofstream swq6;
    swq6.open(fn6,ios::binary);
    if(!swq6)
       {
          cout<<"cannot open "<<fn6<<endl;
          abort();
       }

    ofstream swq7;
    swq7.open(fn7,ios::binary);
    if(!swq7)
       {
          cout<<"cannot open "<<fn7<<endl;
          abort();
       }

    ofstream swq666;
    swq666.open(fn8,ios::binary);
    if(!swq666)
       {
          cout<<"cannot open "<<fn6<<endl;
          abort();
       }
    
    ofstream swq777;
    swq777.open(fn9,ios::binary);
    if(!swq777)
       {
          cout<<"cannot open "<<fn9<<endl;
          abort();
       }

    ofstream swq2; 
    swq2.open(fn2,ios::binary);
    if(!swq2)
       {
          cout<<"cannot open "<<fn2<<endl;
          abort();
       }

    cout<<"Prediction Starts..."<<endl;

    for(it=0;it<ltt/2+1;it++)    
        {

          if(it<itmin)
            {
               for(is=0;is<nr;is++)
                 {
                    mf=(0.0,0.0);  
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));  
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));                
                 }   
               cout<<it<<"frequency slices have been finished!"<<endl;
            }      
          else if(it>itmax&&it<ltt/2+1)
            {
                for(is=0;is<nr;is++)
                 {
                    mf=(0.0,0.0);
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                 } 
                cout<<it<<"frequency slices have been finished!"<<endl;
            } 
 
          else
           { 
             if((it+1)%5==0)
                cout<<it<<"frequency slices have been finished!"<<endl;

//forming the data matrix for a certain frequcncy
               swq1.seekg(it*4,ios::beg);
               swq1.read((char*)&(dd[0].real()),sizeof(dd[0].real()));
               swq11.seekg(it*4,ios::beg);
               swq11.read((char*)&(dd[0].imag()),sizeof(dd[0].imag()));

//               cout<<"222"<<endl;
//               return 0;
                            

               for(is=1;is<nsr;is++)
                {
                  swq1.seekg((ltt-1)*4,ios::cur);  
                  swq1.read((char*)&(dd[is].real()),sizeof(dd[is].real()));
                  swq11.seekg((ltt-1)*4,ios::cur);  
                  swq11.read((char*)&(dd[is].imag()),sizeof(dd[is].imag()));    
                }

               for(is=0;is<ns;is++)
                 {
                   for(ir=0;ir<nr;ir++)
                     {
                        itmp=is*nr+ir;
                        d[is][ir].real()=dd[itmp].real();
                        d[is][ir].imag()=dd[itmp].imag();
                     }
                 }
/*
               for(is=0;is<ns;is++)
                 {
                   for(ir=0;ir<nr;ir++)
                     {
                        if(fabs(d[is][ir].real())>10||fabs(d[is][ir].imag())>10)   
                        d[is][ir].real()=0.0;
                        d[is][ir].imag()=0.0;
                     }
                 }
*/
/*
               for(is=0;is<ns;is++)
                 for(ir=0;ir<nr;ir++)
                     cout<<"("<<is<<","<<ir<<")===="<<d[is][ir]<<endl;
               return 0;
*/
           
//calculating Green's functions of each source (or equally receiver) in kx-w domain;    
           for(is=0;is<nx;is++)
             for(ir=0;ir<nxb;ir++)
               {
                  g_k[is][ir]=(0.0,0.0);    
                  g[is][ir]=(0.0,0.0);  
               }

           for(is2=0;is2<nx;is2++)
            {
              for(is=0;is<nxb;is++)
                 s[is]=(0.0,0.0);
              s[is2*int(dsx/dx)+bou].real()=1.0*ham[it];
              s[is2*int(dsx/dx)+bou].imag()=0.0;

//              gtmp3.real()=0.0;
//              gtmp3.imag()=omega[it]*t0;
//              s[is2*int(dsx/dx)+bou]=rwf[it]*exp(gtmp3);

              if((in4==NULL)||(out4==NULL))
                cout<<"memory insufficient"<<endl;
              else
                {
                  for(is=0;is<nxb;is++)
                    {
                       in4[is][0]=s[is].real();
                       in4[is][1]=s[is].imag();
                    }
                }

              fftwf_execute(p4);

              for(is=0;is<nxb;is++)
                {
                   s_k[is].real()=out4[is][0];
                   s_k[is].imag()=out4[is][1];
                }

              for(is=0;is<nxb;is++)
                for(iz=0;iz<zmax;iz++)
                   {
                      gs_for[is][iz]=(0.0,0.0);
                      gs_for1[is][iz]=(0.0,0.0);
                      gs_back[is][iz]=(0.0,0.0);
                   } 

//forward extrapolation
                for(is=0;is<nxb;is++)
                   gs_for[is][0]=s_k[is]; 
                for(is=0;is<nxb;is++)
                  {
                     if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                        {
                           gs_for[is][1].real()=0.0;
                           gs_for[is][1].imag()=0.0;
                        }

                     else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                        {
                          gtmp1.real()=0.0;               //phase-shift extrapolating operator
                          gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                          gs_for[is][1]=gs_for[is][0]*exp(gtmp1);
                        }
                     else
                        {
                           gs_for[is][1].real()=0.0;
                           gs_for[is][1].imag()=0.0;
                        }
                   }
                if((in7==NULL)||(out7==NULL))
                   cout<<"memory insufficient"<<endl;
                else
                  {
                    for(is=0;is<nxb;is++)
                        {
                           in7[is][0]=gs_for[is][1].real();
                           in7[is][1]=gs_for[is][1].imag();
                        }
                   }

                   fftwf_execute(p7);

                 for(is=0;is<nxb;is++)
                   {
                      uxf[is].real()=out7[is][0]/sqrt(nxb);
                         uxf[is].imag()=out7[is][1]/sqrt(nxb);
                   }

                 for(is=0;is<nxb;is++)
                     uxf[is]=uxf[is]*abso[is];

                 if((in8==NULL)||(out8==NULL))
                    cout<<"memory insufficient"<<endl;
                 else
                   {
                     for(is=0;is<nxb;is++)
                      {
                        in8[is][0]=uxf[is].real();
                        in8[is][1]=uxf[is].imag();
                      }
                   }

                  fftwf_execute(p8);

                  for(is=0;is<nxb;is++)
                     {
                       gs_for[is][1].real()=out8[is][0]/sqrt(nxb);
                       gs_for[is][1].imag()=out8[is][1]/sqrt(nxb);
                     }

                for(iz=2;iz<zmax;iz++)
                   {
                     for(is=0;is<nxb;is++)
                       {
                         if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                           {
                              gs_for[is][iz].real()=0.0;
                              gs_for[is][iz].imag()=0.0;
                           }
            
                         else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                           {
                              gtmp1.real()=0.0;               //phase-shift extrapolating operator
                              gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                              gs_for[is][iz]=gs_for[is][iz-1]*exp(gtmp1);
                           }
                         else
                           {
                              gs_for[is][iz].real()=0.0;
                              gs_for[is][iz].imag()=0.0;
                           }
                       }
                 }//end of iz for all depth  forward extraplotation

/*
                for(is=0;is<nxb;is++)
                   for(iz=0;iz<zmax;iz++)
                        cout<<gs_for[is][iz]<<endl;
                return 0;
*/


//backward extrapolaion
                 for(is=0;is<nxb;is++)
                     for(iz=0;iz<zmax;iz++)
                       gs_for1[is][iz]=gs_for[is][iz]*r1[iz];

                 for(is=0;is<nxb;is++)
                   gs_back[is][0]=gs_for1[is][zmax-1];

                 for(iz=1;iz<zmax-1;iz++)
                   {
                      for(is=0;is<nxb;is++)
                        {
                           if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                              {
                                 gs_back[is][iz].real()=0.0;
                                 gs_back[is][iz].imag()=0.0;
                              }

                           else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                              {
                                   gtmp1.real()=0.0;               //phase-shift extrapolating operator
                                   gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                                   gs_back[is][iz]=gs_back[is][iz-1]*exp(gtmp1)+gs_for1[is][zmax-1-iz];
                              }
                          else
                              {
                    		   gs_back[is][iz].real()=0.0;
                                   gs_back[is][iz].imag()=0.0;
                              }
                          }//end of is
                }//end of iz for back extrapolation
             
/*
               for(is=0;is<nxb;is++)
                   for(iz=0;iz<zmax;iz++)
                        cout<<gs_back[is][iz]<<endl;
                return 0;
*/

                for(is=0;is<nxb;is++)
                   gs_back[is][zmax-1]=gs_back[is][zmax-2]*exp(gtmp1);
//final Green's functions in kx-w domain for a certain (kx,w) pair;
                for(is=0;is<nxb;is++)   
                   g_k[is2][is]=gs_back[is][zmax-1];

             }//end of is2 loop, .ie. Green's function of all shots for a certain frequency 
/*
             for(is2=0;is2<nx;is2++)
                for(is=0;is<nxb;is++)
                   cout<<g_k[is2][is]<<endl;
                return 0;
*/



//ifft the Green's functions for frequency it from kx to x domain 
             for(is=0;is<nx;is++)
               { 
                 if((in3==NULL)||(out3==NULL))
                   cout<<"memory insufficient"<<endl;
                 else
                   {
                      for(ir=0;ir<nxb;ir++)
                        {
                           in3[ir][0]=g_k[is][ir].real();
                           in3[ir][1]=g_k[is][ir].imag();
                        }
                   }
                 fftwf_execute(p3);

                 for(ir=0;ir<nxb;ir++)
                   {
                      g[is][ir].real()=out3[ir][0];
                      g[is][ir].imag()=out3[ir][1];
                   }
                }
//calculating the predicted results in x-w domain;

                for(is1=0;is1<nx;is1++)                    
                   gs[is1]=g[shot-1][is1*int(dsx/dx)+bou]; 

                for(is1=0;is1<nx;is1++)
                   gsd[is1]=(0.0,0.0);
          
                for(ir=0;ir<nx;ir++)
                   {
                      gtmp2.real()=0.0;
                      gtmp2.imag()=0.0;
                      gtmp3.real()=0.0;
                      gtmp3.imag()=0.0;

                      for(ir1=0;ir1<nx;ir1++)
                        {
                          gtmp2+=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou];
                          gtmp3+=d[ir1][ir]*gs[ir1];
                        }
                      mf=gtmp2+gtmp3;

                      mf.real()=mf.real()/nx/nx;
                      mf.imag()=mf.imag()/nx/nx;

                      swq6.write((char*)&mf.real(),sizeof(mf.real()));
                      swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                   } //end of receiver for a certain shot and a certain frequency
            }//end of else,.ie. a certain frequency
         }//end of all positive frequency

        swq6.close();
        swq7.close();

        cout<<"Positive frequcies prediction Done!"<<endl;

       ifstream swq66;
       swq66.open(fn6,ios::binary);
       if(!swq66)
         {
           cout<<"cannot open "<<fn6<<endl;
           abort();
         }

       ifstream swq77;
       swq77.open(fn7,ios::binary);
       if(!swq77)
         {
           cout<<"cannot open "<<fn7<<endl;
           abort();
         }

//start to compute the predicticted model of the negative frequency
        cout<<"Negative frequcies prediction Done!"<<endl;
        cout<<"IFFT Starts..."<<endl;
        
        for(is=0;is<nr;is++)
          {
           swq66.seekg(is*4,ios::beg); 
           swq66.read((char*)&f_tmp[0].real(),sizeof(f_tmp[0].real())); 
           swq77.seekg(is*4,ios::beg);  
           swq77.read((char*)&f_tmp[0].imag(),sizeof(f_tmp[0].imag()));

           for(it=1;it<ltt/2+1;it++)
             {
               swq66.seekg((nr-1)*4,ios::cur);
               swq77.seekg((nr-1)*4,ios::cur);
               swq66.read((char*)&f_tmp[it].real(),sizeof(f_tmp[it].real()));
               swq77.read((char*)&f_tmp[it].imag(),sizeof(f_tmp[it].imag()));
             } 
           for(it=ltt/2+1;it<ltt;it++)
               {
                 f_tmp[it].real()=f_tmp[ltt-it].real();
                 f_tmp[it].imag()=-f_tmp[ltt-it].imag();
               } 
           for(it=0;it<ltt;it++)
               {
                 swq666.write((char*)&f_tmp[it].real(),sizeof(f_tmp[it].real()));
                 swq777.write((char*)&f_tmp[it].imag(),sizeof(f_tmp[it].imag()));
               }
              
           if((in2==NULL)||(out2==NULL))
              cout<<"memory insufficient"<<endl;
           else
              {
                for(it=0;it<ltt;it++)
                  {
                     in2[it][0]=f_tmp[it].real();
                     in2[it][1]=f_tmp[it].imag();
                  }
              }

              fftwf_execute(p2);

              for(it=0;it<ltt;it++)
                 m[it]=-out2[it][0]/ltt;

              for(it=0;it<lt;it++)
                swq2.write((char*)&m[it],sizeof(m[it]));
   
          }

        cout<<"IFFT Ends..."<<endl;
                                       
       
        swq1.close(); 
        swq2.close();
        swq11.close();
        swq66.close();
        swq77.close();
        swq666.close();
        swq777.close();

        cout<<"disk written."<<endl;

//free the memory
        free1float(omega);
        free1float(kx);
        free2complex(d);
        free1float(m);
        free2complex(gs_for);
        free2complex(gs_back);
        free2complex(g_k);
        free2complex(g);
        free1complex(A);
        free1complex(B);
        free1complex(gsd);  

        cout<<"all done!"<<endl;

        return 0;
 
   
}






