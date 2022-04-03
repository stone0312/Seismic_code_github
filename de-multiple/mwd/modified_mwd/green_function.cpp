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
   char fn3[256];
   int nx,ns,nr,lt,ltt,nsr,zmax,zmin,bou,nxb;
   float dt,dsx,drx,fmin,fmax,dx,dz,v;
   int itmin,itmax,shot;

   int is,ir,it,iz,is1,ir1,is2,ir2;

    ifstream swq;
    swq.open("green_function.par");
    if(!swq)
      {
          cout<<"cannot open modified_mwd_v11.par"<<endl;
          abort();
      }
    swq>>fn3>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou>>shot;
    swq.close();

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
    cout<<"half of width of boundary is==== "<<bou<<endl;
    cout<<"No. of shot to be predicted is===="<<shot<<endl;

    ltt=2*lt; 
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
    zmax=zmax+1;

    complex<float> *r1;
    r1=alloc1complex(zmax);
    for(is=0;is<zmax;is++)
    {
         r1[is].real()=0.0;
         r1[is].imag()=0.0;
    }
//      r1[is]=(0.0,0.0);

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
    {
         s[is].real()=0.0;
         s[is].imag()=0.0;
    }
//       s[is]=(0.0,0.0);

    complex<float> *r;
    r=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
    {
         r[is].real()=0.0;
         r[is].imag()=0.0;
    }
//       r[is]=(0.0,0.0);

    complex<float> *s_k;
    s_k=alloc1complex(nxb);
 
    complex<float> *s_k1;
    s_k1=alloc1complex(nxb);

    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nxb);

    complex <float> **gs_for;
    gs_for=alloc2complex(zmax,nxb);

    complex <float> **gs_for1;
    gs_for1=alloc2complex(zmax,nxb);
 
    complex <float> **gs_back;
    gs_back=alloc2complex(zmax,nxb); 

    complex <float> gtmp1;   
    gtmp1=(0.0,0.0);

//define the plans and arrays when fftw is complemented
    fftwf_complex *in4,*out4;
    fftwf_plan p4;
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
 
    p4=fftwf_plan_dft_1d(nxb,in4,out4,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in5,*out5,*in6,*out6;
    fftwf_plan p5,p6;
    in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    p5=fftwf_plan_dft_1d(nxb,in5,out5,FFTW_BACKWARD,FFTW_MEASURE);
    p6=fftwf_plan_dft_1d(nxb,in6,out6,FFTW_BACKWARD,FFTW_MEASURE);
    
    complex <float> *green_for;
    green_for=alloc1complex(nxb);
    complex <float> *green_back;
    green_back=alloc1complex(nxb);

//zeroing all arrays;
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
   
   for(is=0;is<nxb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dx*nxb); 
   for(is=nxb/2+1;is<nxb;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nxb/2)/float(dx*nxb);

    ofstream green1_real;
    green1_real.open("check_green_for_real.dat",ios::binary);
    if(!green1_real)
       {
          cout<<"cannot open check_green_for_real.dat"<<endl;
          abort();
       }
    ofstream green1_imag;
    green1_imag.open("check_green_for_imag.dat",ios::binary);
    if(!green1_imag)
       {
          cout<<"cannot open check_green_for_imag.dat"<<endl;
          abort();
       }
    ofstream green2_real;
    green2_real.open("check_green_back_real.dat",ios::binary);
    if(!green2_real)
       {
          cout<<"cannot open check_green_back_real.dat"<<endl;
          abort();
       }
    ofstream green2_imag;
    green2_imag.open("check_green_back_imag.dat",ios::binary);
    if(!green2_imag)
       {
          cout<<"cannot open check_green_back_imag.dat"<<endl;
          abort();
       }

    cout<<"Calculation of Green's Function Starts..."<<endl;

   for(it=0;it<ltt;it++)
        {
/*   
           cout<<omega[it]<<endl;
           return 0;
*/
//calculating Green's functions of each source (or equally receiver) in kx-w domain;    
           for(is2=0;is2<1;is2++)
            {
              for(is=0;is<nxb;is++)
                {
                 s[is].real()=0.0;
                 s[is].imag()=0.0;
                }
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
/*
              for(is=0;is<nxb;is++)
                 cout<<is<<"  "<<s_k[is]<<endl; 
              return 0;
*/
/*
              if((in6==NULL)||(out6==NULL))
                cout<<"memory insufficient"<<endl;
              else
                {
                  for(is=0;is<nxb;is++)
                    {
                       in6[is][0]=s_k[is].real();
                       in6[is][1]=s_k[is].imag();
                    }
                }

              fftwf_execute(p6);

              for(is=0;is<nxb;is++)
                {
                   s_k1[is].real()=out6[is][0]/nxb;
                   s_k1[is].imag()=out6[is][1]/nxb;
                }


              for(is=0;is<nxb;is++)
                 cout<<is<<"  "<<s_k1[is]<<endl; 
              return 0;
*/
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
                       for(iz=0;iz<zmax;iz++)
                          {
                             gs_for[is][iz].real()=0.0;
                             gs_for[is][iz].imag()=0.0;
                             gs_back[is][iz].real()=0.0;
                             gs_back[is][iz].imag()=0.0; 
                          }       
                     }
/*
                   else if(pow(omega[it]/v,2)-pow(kx[is],2)==0.0)
                     {
                       for(iz=0;iz<zmax;iz++)
                          {
                             gs_for[is][iz].real()=0.0;
                             gs_for[is][iz].imag()=0.0;
                             gs_back[is][iz].real()=0.0;
                             gs_back[is][iz].imag()=0.0;
                          }
                     }
*/
                   else
                    {
                     gtmp1.real()=0.0;//phase-shift extrapolating operator
                     gtmp1.imag()=-dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));

//forward extrapolation
                     gs_for[is][0].real()=s_k[is].real();
                     gs_for[is][0].imag()=s_k[is].imag();
                     
                     for(iz=1;iz<zmax;iz++)
                         gs_for[is][iz]=gs_for[is][iz-1]*exp(gtmp1);

//backward extrapolaion
                     for(iz=0;iz<zmax;iz++)
                        gs_for1[is][iz]=gs_for[is][iz]*r1[iz];

                     gs_back[is][0]=gs_for1[is][zmax-1];

                     for(iz=1;iz<zmax;iz++)
                         gs_back[is][iz]=gs_back[is][iz-1]*exp(gtmp1)+gs_for1[is][zmax-1-iz];
                    }
                 }//end of is, .ie. Green's function of a shot for a certain frequency
               }//end of is2 
             }// end of it 

           for(iz=0;iz<zmax;iz++)
             {   
              if((in5==NULL)||(out5==NULL))
                cout<<"memory insufficient"<<endl;
              else
                {
                  for(is=0;is<nxb;is++)
                    {
                       in5[is][0]=gs_for[is][iz].real();
                       in5[is][1]=gs_for[is][iz].imag();
                    }
                }

              fftwf_execute(p5);

              for(is=0;is<nxb;is++)
                {
                   green_for[is].real()=out5[is][0]/nxb;
                   green_for[is].imag()=out5[is][1]/nxb;
                }
  
              for(is=0;is<nxb;is++)
                {
                   green1_real.write((char*)&green_for[is].real(),sizeof(green_for[is].real()));
                   green1_imag.write((char*)&green_for[is].imag(),sizeof(green_for[is].imag()));
                }
 
              if((in6==NULL)||(out6==NULL))
                cout<<"memory insufficient"<<endl;
              else
                {
                  for(is=0;is<nxb;is++)
                    {
                       in6[is][0]=gs_back[is][iz].real();
                       in6[is][1]=gs_back[is][iz].imag();
                    }
                }

              fftwf_execute(p6);

              for(is=0;is<nxb;is++)
                {
                   green_back[is].real()=out6[is][0]/nxb;
                   green_back[is].imag()=out6[is][1]/nxb;
                }
              for(is=0;is<nxb;is++)
                {
                   green2_real.write((char*)&green_back[is].real(),sizeof(green_back[is].real()));
                   green2_imag.write((char*)&green_back[is].imag(),sizeof(green_back[is].imag()));

                }
              }
             cout<<"calculation of green's function for a certain frequency done!"<<endl;
             return 0;
}

