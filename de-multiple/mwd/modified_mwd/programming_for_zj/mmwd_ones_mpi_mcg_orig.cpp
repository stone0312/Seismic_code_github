#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#include "mpi.h"
#define pai 3.14159265

int main(int argc,char **argv)
{
   int np,myid;	 
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fn12[256],fn13[256],fn14[256];
   int nx,ns,nr,lt,ltt,nsr,zmax,zmin,bou,nxb,wid;
   float dt,dsx,drx,fmin,fmax,dx,dz,v,dip_max;
   int itmin,itmax,shot,len;

   int icount,icount1,nt,abs1;
   float st_beg,st_end;   

   int is,ir,it,iz,is1,ir1,is2,ir2,itmp,ix,tmp,tmp1;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&np);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   
   
    ifstream swq;
    swq.open("mmwd_ones_mpi_mcg.par");
    if(!swq)
      {
          cout<<"cannot open mmwd_ones_mpi_mcg.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>fn11>>fn12>>fn13>>fn14>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou>>shot>>dip_max>>wid>>len>>nt>>st_beg>>st_end;
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
    nxb=nx*(int)dsx/(int)dx+2*bou;

    icount=(1+(nr-(shot-1)))*(nr-(shot-1))/2+((1+((shot-1)+1))*((shot-1)+1)/2-1);
   cout<<"No. of total traces for MCG of shot "<<shot<<" is===="<<icount<<endl;

   if(nt==1)
     icount1=0;
   else if(nt<shot)
     icount1=(shot+(shot-(nt-2)))*(nt-1)/2;
   else if(nt==shot)
     icount1=(shot+2)*(shot-1)/2;
   else if (nt==shot+1)
     icount1=(shot+1)*shot/2;
   else
     icount1=(shot+1)*shot/2+(2+(shot-nt))*((shot-nt)-2+1)/2;


   abs1=abs((nt-1)-(shot-1))+1;

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
    r1=alloc1complex(zmax-zmin+1);
    for(is=0;is<zmax-zmin+1;is++)
      r1[is]=(0.0,0.0);

    if(zmin==zmax-1)
      {
         r1[zmax-zmin].real()=1.0;
         r1[zmax-zmin].imag()=0.0;
      }
    else
      {
         for(iz=1;iz<zmax-zmin+1;iz++)
           {
             for(is=0;is<nx;is++)
              {
                if(dep[is]==iz+zmin-1)
                 {
                   r1[iz].real()=1.0;
                   r1[iz].imag()=0.0;
                   is=nx;
                 }
                else
                 {
                   r1[iz].real()=0.0;
                   r1[iz].imag()=0.0;
                 }

              }
           }
      }
/*
   for(iz=0;iz<zmax-zmin+1;iz++)
      cout<<iz<<r1[iz]<<endl;
   return 0;
*/
   complex<float> mcg1;

   complex<float> *mcg;
   mcg=alloc1complex(icount);

   complex<float> *mcg2;
   mcg2=alloc1complex(ltt);

   float *mcg3;
   mcg3=alloc1float(ltt);

   float **mcg_cer;
   mcg_cer=alloc2float(lt,abs1);

   float *srm1;
   srm1=alloc1float(lt);

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

    float *abso;
    abso=alloc1float(nxb);
    
    for(is=0;is<nxb;is++)
       abso[is]=1.0;
    
    for(is=0;is<wid;is++)
       abso[is]=sqrt(sin(pai/2*is/(wid-1)));
    
    for(is=nxb-wid;is<nxb;is++)
        abso[is]=abso[nxb-is];
      
    complex <float> **gs_for;
    gs_for=alloc2complex(zmax-zmin+1,nxb);

    complex <float> **gs_for1;
    gs_for1=alloc2complex(zmax-zmin+1,nxb);
 
    complex <float> **gs_back;
    gs_back=alloc2complex(zmax-zmin+1,nxb); 

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
    for(iz=0;iz<zmax-zmin+1;iz++)
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

    fstream swq6;
    swq6.open(fn6,ios::binary|ios::out);
    if(!swq6)
       {
          cout<<"cannot open "<<fn6<<endl;
          abort();
       }

    fstream swq7;
    swq7.open(fn7,ios::binary|ios::out);
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

    fstream swq10;
    swq10.open(fn10,ios::binary|ios::out);
    if(!swq10)
       {
          cout<<"cannot open "<<fn10<<endl;
          abort();
       }

    fstream swq111;
    swq111.open(fn11,ios::binary|ios::out);
    if(!swq111)
       {
          cout<<"cannot open "<<fn11<<endl;
          abort();
       }

    fstream swq12;
    swq12.open(fn12,ios::binary|ios::out);
    if(!swq12)
       {
          cout<<"cannot open "<<fn12<<endl;
          abort();
       }

    cout<<"Prediction Starts..."<<endl;
	
    for(it=myid;it<ltt/2+1;it+=np)    
        {

          if(it<itmin)
            {
	       swq6.seekg(it*nr*4,ios::beg);
	       swq7.seekg(it*nr*4,ios::beg);
			   
               for(is=0;is<nr;is++)
                 {
                    mf=(0.0,0.0);  
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));  
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));                
                 }
 
               swq10.seekg(it*icount*4,ios::beg);
               swq111.seekg(it*icount*4,ios::beg);  
                   
               for(is=0;is<icount;is++)
                 {
                    mcg1=(0.0,0.0);
                    swq10.write((char*)&mcg1.real(),sizeof(mcg1.real()));
                    swq111.write((char*)&mcg1.imag(),sizeof(mcg1.imag()));
                 }
   
               cout<<it<<"frequency slices have been finished!"<<endl;
            }      
          else if(it>itmax&&it<ltt/2+1)
            {
		swq6.seekg(it*nr*4,ios::beg);
		swq7.seekg(it*nr*4,ios::beg);
					
                for(is=0;is<nr;is++)
                 {
                    mf=(0.0,0.0);
                    swq6.write((char*)&mf.real(),sizeof(mf.real()));
                    swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                 }

               swq10.seekg(it*icount*4,ios::beg);
               swq111.seekg(it*icount*4,ios::beg);

               for(is=0;is<icount;is++)
                 {
                    mcg1=(0.0,0.0);
                    swq10.write((char*)&mcg1.real(),sizeof(mcg1.real()));
                    swq111.write((char*)&mcg1.imag(),sizeof(mcg1.imag()));
                 }
 
             }
	  else
            { 
             if((it+1)%5==0)
                cout<<it<<"frequency slices have been finished!"<<endl;

             swq6.seekg(it*nr*4,ios::beg);
	     swq7.seekg(it*nr*4,ios::beg);
			 
             swq10.seekg(it*icount*4,ios::beg);
             swq111.seekg(it*icount*4,ios::beg);

//forming the data matrix for a certain frequcncy
               swq1.seekg(it*4,ios::beg);
               swq1.read((char*)&(dd[0].real()),sizeof(dd[0].real()));
               swq11.seekg(it*4,ios::beg);
               swq11.read((char*)&(dd[0].imag()),sizeof(dd[0].imag()));

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
                for(iz=0;iz<zmax-zmin+1;iz++)
                   {
                      gs_for[is][iz]=(0.0,0.0);
                      gs_for1[is][iz]=(0.0,0.0);
                      gs_back[is][iz]=(0.0,0.0);
                   } 

//forward extrapolation
                for(is=0;is<nxb;is++)
                   gs_for[is][0]=s_k[is];

//                for(is=nxb/2+1;is<nxb;is++)
//                   gs_for[is][0]=(0.0,0.0);

 
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
                          gtmp1.imag()=-zmin*dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
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
                      uxf[is].real()=out7[is][0];
                      uxf[is].imag()=out7[is][1];
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
                       gs_for[is][1].real()=out8[is][0]/nxb;
                       gs_for[is][1].imag()=out8[is][1]/nxb;
                     }

                for(iz=2;iz<zmax-zmin+1;iz++)
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

//backward extrapolaion
                 for(is=0;is<nxb;is++)
                     for(iz=0;iz<zmax-zmin+1;iz++)
                       gs_for1[is][iz]=gs_for[is][iz]*r1[iz];

                 for(is=0;is<nxb;is++)
                   gs_back[is][0]=gs_for1[is][zmax-zmin];

                 for(iz=1;iz<zmax-zmin;iz++)
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
                                   gs_back[is][iz]=gs_back[is][iz-1]*exp(gtmp1)+gs_for1[is][zmax-zmin-iz];
                              }
                          else
                              {
                    		   gs_back[is][iz].real()=0.0;
                                   gs_back[is][iz].imag()=0.0;
                              }
                          }//end of is
                }//end of iz for back extrapolation
             
                for(is=0;is<nxb;is++)
                 {
                    if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                      {
                         gs_back[is][zmax-zmin].real()=0.0;
                         gs_back[is][zmax-zmin].imag()=0.0;
                      } 

                     else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                      {
                        gtmp1.real()=0.0;               //phase-shift extrapolating operator
                        gtmp1.imag()=-zmin*dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                        gs_back[is][zmax-zmin]=gs_back[is][zmax-zmin-1]*exp(gtmp1);
                      }
                    else
                     {
                         gs_back[is][zmax-zmin].real()=0.0;
                         gs_back[is][zmax-zmin].imag()=0.0;
                     }
                 }
//final Green's functions in kx-w domain for a certain (kx,w) pair;
                for(is=0;is<nxb;is++)   
                   g_k[is2][is]=gs_back[is][zmax-zmin];

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
                      g[is][ir].real()=out3[ir][0]/nxb;
                      g[is][ir].imag()=out3[ir][1]/nxb;
                   }
                }
//calculating the predicted results in x-w domain;
                for(ir=0;ir<nx;ir++)
                   {
                      gtmp2.real()=0.0;
                      gtmp2.imag()=0.0;
                      gtmp3.real()=0.0;
                      gtmp3.imag()=0.0;

                      if(ir>shot-1)
                       {
                         for(ir1=shot-1;ir1<ir+1;ir1++)
                          {
                            mcg1=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou]+d[ir1][ir]*g[shot-1][ir1*int(dsx/dx)+bou];

//                            mcg1=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou];

                            swq10.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                            swq111.write((char*)&(mcg1.imag()),sizeof(mcg1.imag())); 
                            
                            gtmp2+=mcg1;

//                            gtmp2+=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou];
//                            gtmp3+=d[ir1][ir]*g[shot-1][ir1*int(dsx/dx)+bou];
                          }
//                         mf=gtmp2+gtmp3;

                         mf=gtmp2;

                         mf.real()=mf.real()/(ir-shot+2);
                         mf.imag()=mf.imag()/(ir-shot+2);

                         swq6.write((char*)&mf.real(),sizeof(mf.real()));
                         swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                       }
                     else if(ir<shot-1)
                       {
                         for(ir1=ir;ir1<shot;ir1++)
                           {
                            mcg1=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou]+d[ir1][ir]*g[shot-1][ir1*int(dsx/dx)+bou];

//                            mcg1=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou];

                            swq10.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                            swq111.write((char*)&(mcg1.imag()),sizeof(mcg1.imag())); 

//                            gtmp2+=d[shot-1][ir1]*g[ir][ir1*int(drx/dx)+bou];
//                            gtmp3+=d[ir1][ir]*g[shot-1][ir1*int(dsx/dx)+bou];

                            gtmp2+=mcg1;

                           }

                           mf=gtmp2; 
//                         mf=gtmp2+gtmp3;

                         mf.real()=mf.real()/(shot-ir);
                         mf.imag()=mf.imag()/(shot-ir);

                         swq6.write((char*)&mf.real(),sizeof(mf.real()));
                         swq7.write((char*)&mf.imag(),sizeof(mf.imag()));
                       }

                     else
                      {
                         mcg1=d[shot-1][shot-1]*g[ir][(shot-1)*int(drx/dx)+bou]+d[shot-1][ir]*g[shot-1][ir*int(dsx/dx)+bou];

//                         mcg1=d[shot-1][shot-1]*g[ir][(shot-1)*int(drx/dx)+bou];

                         swq10.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                         swq111.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));

//                         gtmp2+=d[shot-1][shot-1]*g[ir][(shot-1)*int(drx/dx)+bou];
//                         gtmp3+=d[shot-1][ir]*g[shot-1][(shot-1)*int(dsx/dx)+bou];

                         mf=mcg1;
//                         mf=gtmp2+gtmp3;

                         swq6.write((char*)&mf.real(),sizeof(mf.real()));
                         swq7.write((char*)&mf.imag(),sizeof(mf.imag()));

                      }

                   } //end of receiver for a certain shot and a certain frequency
            }//end of else,.ie. a certain frequency
         }//end of all positive frequency

        swq6.close();
        swq7.close();
        swq10.close();
        swq111.close(); 
	
	MPI_Finalize();

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

       ifstream swq88;
       swq88.open(fn10,ios::binary);
       if(!swq88)
         {
            cout<<"cannot open "<<fn10<<endl;
            abort();
         }

       ifstream swq99;
       swq99.open(fn11,ios::binary);
       if(!swq99)
         {
             cout<<"cannot open "<<fn11<<endl;
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

        cout<<"MCG Forming Starts..."<<endl;

        for(is=0;is<icount;is++)
            {
                swq88.seekg(is*4,ios::beg);
                swq88.read((char*)&(mcg2[0].real()),sizeof(mcg2[0].real()));
                swq99.seekg(is*4,ios::beg);
                swq99.read((char*)&(mcg2[0].imag()),sizeof(mcg2[0].imag()));

                for(it=1;it<lt+1;it++)
                   {
                      swq88.seekg((icount-1)*4,ios::cur);
                      swq88.read((char*)&(mcg2[it].real()),sizeof(mcg2[it].real()));
                      swq99.seekg((icount-1)*4,ios::cur);
                      swq99.read((char*)&(mcg2[it].imag()),sizeof(mcg2[it].imag()));
                   }
                for(it=lt+1;it<ltt;it++)
                   {
                      mcg2[it].real()=mcg2[ltt-it].real();
                      mcg2[it].imag()=-mcg2[ltt-it].imag();
                   }

                if((in2==NULL)||(out2==NULL))
                   cout<<"memory insufficient"<<endl;
                else
                  {
                     for(it=0;it<ltt;it++)
                       {
                         in2[it][0]=mcg2[it].real();
                         in2[it][1]=mcg2[it].imag();
                       }
                  }

                fftwf_execute(p2);

                for(it=0;it<ltt;it++)
                   mcg3[it]=out2[it][0];

                for(it=0;it<lt;it++)
                   swq12.write((char*)&(mcg3[it]),sizeof(mcg3[it]));

            }
 
        cout<<"MCG has been formed..."<<endl;
        cout<<"IFFT Ends..."<<endl;                               
       
        swq1.close(); 
        swq2.close();
        swq11.close();
        swq66.close();
        swq77.close();
        swq666.close();
        swq777.close();
        swq12.close();
        swq88.close();
        swq99.close();

        ifstream swq200;
          swq200.open(fn12,ios::binary);
          if(!swq200)
            {
              cout<<"cannot open "<<fn12<<endl;
              abort();
            }

           ofstream swq13;
           swq13.open(fn13,ios::binary);
           if(!swq13)
             {
                cout<<"cannot open "<<fn13<<endl;
                abort();
             }

           ifstream swq144;
           swq144.open(fn2,ios::binary);
           if(!swq144)
             {
                cout<<"cannot open "<<fn2<<endl;
                abort();
             }

           ofstream swq14;
           swq14.open(fn14,ios::binary);
           if(!swq14)
             {
                cout<<"cannot open "<<fn14<<endl;
                abort();
             }

           swq144.seekg(0,ios::beg);
           for(ir=0;ir<nt-1;ir++)
               swq144.seekg(lt*4,ios::cur);
           for(it=0;it<lt;it++)
              swq144.read((char*)&srm1[it],sizeof(srm1[it]));

           for(it=0;it<lt;it++)
              swq14.write((char*)&srm1[it],sizeof(srm1[it]));

           swq200.seekg(0,ios::beg);
           for(ir=0;ir<icount1;ir++)
               swq200.seekg(lt*4,ios::cur);

           for(ir=0;ir<abs1;ir++)
              for(it=0;it<lt;it++)
                 swq200.read((char*)&mcg_cer[ir][it],sizeof(mcg_cer[ir][it]));

           for(ir=0;ir<abs1;ir++)
              for(it=0;it<lt;it++)
                 swq13.write((char*)&mcg_cer[ir][it],sizeof(mcg_cer[ir][it]));


          swq200.close();
          swq13.close();
          swq14.close();
          swq144.close();

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






