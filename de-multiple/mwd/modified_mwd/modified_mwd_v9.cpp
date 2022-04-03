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
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256];
   int nx,ns,nr,lt,nsr,zmax,zmin,bou,nxb;
   float dt,dsx,drx,fmin,fmax,dx,dz,v;
   int itmin,itmax;

   int is,ir,it,iz,is1,ir1,is2,ir2;

    ifstream swq;
    swq.open("modified_mwd_v9.par");
    if(!swq)
      {
          cout<<"cannot open modified_mwd_v9.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou;
    swq.close();

    cout<<"fna of original data is==== "<<fn1<<endl;
    cout<<"fna of predicted model is==== "<<fn2<<endl;
    cout<<"fna of water layer depth is==== "<<fn3<<endl;
    cout<<"fna of real part of fft of original data is==== "<<fn4<<endl;
    cout<<"fna of imaginary part of fft of original data is==== "<<fn5<<endl;
    cout<<"fna of real part of fft of predicted model is==== "<<fn6<<endl;
    cout<<"fna of imaginary part of fft of predicted model is==== "<<fn7<<endl;
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

    nsr=ns*nr;
    nxb=nx+2*bou;

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

 //   return 0;

    complex<float> *r1; //define the reflectivity of each layer in the x-f domain 
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

    complex<float> *s; // define the wavelets of source-side in the x-f domain
    s=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
       s[is]=(0.0,0.0);

    complex<float> *r; // define the wavelets of receiver-side in the x-f domain
    r=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
       r[is]=(0.0,0.0);

    complex<float> *s_k; // define the wavelets of source-side in the kx-f domain
    s_k=alloc1complex(nxb);
 
    float *omega;
    omega=alloc1float(lt);
    float *kx;
    kx=alloc1float(nxb);

    complex <float> **d; // define the data matrix for a certain frequency 
    d=alloc2complex(nr,ns);

    complex <float> mf; //define the predicted model result for a certain source and receiver pair at a certain frequency
  
    float *m;  // define the final predicted model result for a certain trace
    m=alloc1float(lt);

    complex <float> **g_k;
    g_k=alloc2complex(nxb,nx);

    complex <float> **g;
    g=alloc2complex(nxb,nx);

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

    complex <float> *gsd;
    gsd=alloc1complex(nx);

    gtmp1=(0.0,0.0);
    gtmp2=(0.0,0.0);

    for(is=0;is<nx;is++)
      gsd[is]=(0.0,0.0);

//define the plans and arrays when fftw is complemented
    fftwf_complex *in2,*out2,*in3,*out3,*in4,*out4;
    fftwf_plan p2,p3,p4;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
 
    p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(nxb,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(nxb,in4,out4,FFTW_FORWARD,FFTW_MEASURE);

    float mem;
    mem=0.0;
    mem+=(nx*nxb*2+lt+lt+nx*4+ns*nr*2+nxb*6+zmax*nxb*4)*4.0/1024.0/1024.0;
    cout<<"Memory needed to be allocated is==== "<<mem<<"MB"<<endl;

//zeroing all arrays;
      for(it=0;it<lt;it++)
        m[it]=0.0;

   for(is=0;is<nxb;is++)
    for(iz=0;iz<zmax;iz++)
      {
          gs_for[is][iz]=(0.0,0.0);
          gs_back[is][iz]=(0.0,0.0);
          gs_for1[is][iz]=(0.0,0.0);
      } 

//calculate the omega and kx
   for(it=0;it<lt;it++)
      omega[it]=2*pai*it*1000/(dt*lt);
   for(it=lt;it<lt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-lt)*1000/(dt*lt));

   itmin=int(fmin*dt*lt/1000);
   itmax=int(fmax*dt*lt/1000)+1;

   cout<<"totally "<<2*(itmax-itmin)<<" frequency slices need to be calculated"<<endl;
   
   for(is=0;is<nxb/2;is++)
      kx[is]=2*pai*float(is)/float(dx*nxb); 
   for(is=nxb/2;is<nxb;is++)
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

    ofstream swq2; 
    swq2.open(fn2,ios::binary);
    if(!swq2)
       {
          cout<<"cannot open "<<fn2<<endl;
          abort();
       }

    for(it=0;it<lt;it++)    
        {
          if(it<itmin)
            {
               for(is=0;is<nsr;is++)
                 {
                    mf=(0.0,0.0);  
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));  
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));                
                 }   
            }      
          else if(it>itmax&&it<lt-itmax)
            {
                for(is=0;is<nsr;is++)
                 {
                    mf=(0.0,0.0);
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                 } 
            } 
          else if(it>lt-itmin&&it<lt)
            {
               for(is=0;is<nsr;is++)
                 {
                    mf=(0.0,0.0);
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                 }
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
                      swq1.seekg(((is*nr+ir)*lt+it)*4,ios::beg);  
                      swq11.seekg(((is*nr+ir)*lt+it)*4,ios::beg);  
                      swq1.read((char*)&(d[is][ir].real()),sizeof(d[is][ir].real()));
                      swq11.read((char*)&(d[is][ir].imag()),sizeof(d[is][ir].imag()));    
                    }                
                }
               
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
              s[is2+bou].real()=1.0;
              s[is2+bou].imag()=0.0;

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

              for(is=0;is<nxb;is++)//wavenumber loop 
                 {
                   if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                     {
                       g_k[is2][is].real()=0.0;
                       g_k[is2][is].imag()=0.0;
                     }
                   else
                    {
                     gtmp1.real()=0.0;//phase-shift extrapolating operator
                     gtmp1.imag()=dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
//forward extrapolation
                     gs_for[is][0].real()=s_k[is].real();
                     gs_for[is][0].imag()=s_k[is].imag();
                     
                     for(iz=1;iz<zmax;iz++)
                         gs_for[is][iz]=gs_for[is][iz-1]*exp(gtmp1);

//backward extrapolaion
                     for(iz=0;iz<zmax;iz++)
                        gs_for1[is][iz]=gs_for[is][iz-1]*r1[iz];

                     gs_back[is][0]=gs_for1[is][zmax-1];

                     for(iz=1;iz<zmax-1;iz++)
                         gs_back[is][iz]=gs_back[is][iz-1]*exp(gtmp1)+gs_for1[is][zmax-1-iz];
                     gs_back[is][zmax-1]=gs_back[is][zmax-2]*exp(gtmp1);

//final Green's functions in kx-w domain for a certain (kx,w) pair;
                     g_k[is2][is]=gs_back[is][zmax-1];
                    }
                 }//end of is, .ie. Green's function of a shot for a certain frequency
             }//end of is2 loop, .ie. Green's function of all shots for a certain frequency 
              
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
          for(is2=0;is2<nx;is2++)
            {
               for(is1=0;is1<nx;is1++)
                  A[is1]=(0.0,0.0);
                A[is2].real()=1.0;
                A[is2].imag()=0.0;

                for(is1=0;is1<nx;is1++)                    
                   gs[is1]=g[is2][is1+bou]+A[is1]; 

                for(ir=0;ir<nx;ir++)
                   {
                      mf=(0.0,0.0);

                      for(is1=0;is1<nx;is1++)
                        gsd[is1]=(0.0,0.0);
                      gtmp2=(0.0,0.0);

                      for(ir1=0;ir1<nx;ir1++)
                         B[ir1]=(0.0,0.0);
                      B[ir].real()=1.0;
                      B[ir].imag()=0.0;      

                      for(ir1=0;ir1<nx;ir1++)
                         gr[ir1]=g[ir][ir1+bou]+B[ir1];
                       
                      for(ir1=0;ir1<nx;ir1++)
                         {
                            for(is1=0;is1<nx;is1++)
                              gsd[ir1]+=gs[is1]*d[is1][ir1];
                         } 
                      for(ir1=0;ir1<nx;ir1++)
                         gtmp2+=gsd[ir1]*gr[ir1];
                      
                      mf=gtmp2-d[is2][ir];
                      mf.real()=mf.real()/nx/nx;
                      mf.imag()=mf.imag()/nx/nx;
 
                      swq6.write((char*)&mf.real(),sizeof(mf.real()));
                      swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                       
                      gtmp2=(0.0,0.0);

                    }//end of receiver for a certain shot and a certain frequency
               } //end of shot for  a certain frequency
            }//end of else,.ie. a certain frequency
         }//end of all frequency

//ifft the (x,w) domain predicted multiple to (x,t) domain
        for(is=0;is<nsr;is++)
           {
              if((in2==NULL)||(out2==NULL))
                 cout<<"memory insufficient"<<endl;
              else
                 {
                   for(it=0;it<lt;it++)
                     {
                        swq66.seekg((it*nsr+is)*4,ios::beg); 
                        swq77.seekg((it*nsr+is)*4,ios::beg);
                        swq66.read((char*)&in2[it][0],sizeof(in2[it][0]));    
                        swq77.read((char*)&in2[it][1],sizeof(in2[it][1]));                   
                     }
                 }

              fftwf_execute(p2);

              for(it=0;it<lt;it++)
                 m[it]=out2[it][0];

              for(it=0;it<lt;it++)
                swq2.write((char*)&m[it],sizeof(m[it])); 
                                       
            }
       
        swq1.close(); 
        swq2.close();
        swq6.close();
        swq7.close();
        swq11.close();
        swq66.close();
        swq77.close();

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






