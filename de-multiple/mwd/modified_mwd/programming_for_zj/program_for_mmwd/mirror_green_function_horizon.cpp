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
   int np,myid;	 
   char fn1[256],fn2[256],fn3[256],fn4[256];
   int nx,ns,nr,lt,ltt,nsr,zmax,zmin,bou,nxb,wid,nxx;
   float dt,dsx,drx,fmin,fmax,dx,dz,v,dip_max;
   int itmin,itmax,shot,len;

   int icount,icount1,nt,abs1;
   float st_beg,st_end;   

   int is,ir,it,iz,is1,ir1,is2,ir2,itmp,ix,tmp,tmp1;
   
    ifstream swq;
    swq.open("mirror_green_function.par");
    if(!swq)
      {
          cout<<"cannot open mirror_green_function.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou>>shot>>dip_max>>wid>>len>>nt>>st_beg>>st_end;
    swq.close();

    cout<<"fna of water layer depth is==== "<<fn1<<endl;
    cout<<"fna of real parts of green's function for positive frequency is==== "<<fn2<<endl;
    cout<<"fna of imaginary parts of green's function for positive frequency is==== "<<fn3<<endl;
    cout<<"fna of final green's function is==== "<<fn4<<endl;
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
    nx=2*nx;
    nxb=nx+2*bou;
    nsr=nx*nx;

//ordering the depth of water layer to find the maximum depth
    int *dep;
    dep=alloc1int(nx/2);

    int *mir_dep;
    mir_dep=alloc1int(nx/2);

    ifstream swq1;
    swq1.open(fn1,ios::binary);
    if(!swq1)
       {
          cout<<"cannot open "<<fn1<<endl;
          abort();
       } 
    for(is=0;is<nx/2;is++)
       swq1.read((char*)&dep[is],sizeof(dep[is]));
    swq1.close();

    for(ix=0;ix<nx/2;ix++)
      mir_dep[ix]=2*dep[ix];

    complex<float> *s;
    s=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
       s[is]=(0.0,0.0);

    complex<float> *s_k;
    s_k=alloc1complex(nxb);
 
    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nxb);

    complex <float> *g_tmp;
    g_tmp=alloc1complex(nxb);

    complex <float> *g_k; 
    g_k=alloc1complex(nxb);

    complex <float> **g;
    g=alloc2complex(nx/2,nx/2);

    complex <float> *gf;
    gf=alloc1complex(ltt);

    float *gt;
    gt=alloc1float(ltt);

    complex<float> gtmp1;

//define the plans and arrays when fftw is complemented
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3;
    fftwf_plan p1,p2,p3;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
 
    p1=fftwf_plan_dft_1d(nxb,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE);
    p2=fftwf_plan_dft_1d(nxb,in2,out2,FFTW_BACKWARD,FFTW_ESTIMATE);
    p3=fftwf_plan_dft_1d(ltt,in3,out3,FFTW_BACKWARD,FFTW_ESTIMATE);


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

   for(is=0;is<nxb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dx*nxb); 
   for(is=nxb/2+1;is<nxb;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nxb/2)/float(dx*nxb);
 
//read the original data

    ofstream swq2;
    swq2.open(fn2,ios::binary);
    if(!swq2)
       {
          cout<<"cannot open "<<fn2<<endl;
          abort();
       }
    ofstream swq3;
    swq3.open(fn3,ios::binary);
    if(!swq3)
       {
          cout<<"cannot open "<<fn3<<endl;
          abort();
       }

    cout<<"Prediction Starts..."<<endl;
	
    for(it=0;it<ltt/2+1;it++)    
        {

          if(it<itmin)
            {
               for(is=0;is<nx/2;is++)
                 for(ir=0;ir<nx/2;ir++) 
                 {
                    g[is][ir]=(0.0,0.0);  
                    swq2.write((char*)&g[is][ir].real(),sizeof(g[is][ir].real()));  
                    swq3.write((char*)&g[is][ir].imag(),sizeof(g[is][ir].imag()));                
                 }
               cout<<it<<"frequency slices have been finished!"<<endl;
            }      
          else if(it>itmax&&it<ltt/2+1)
            {
               for(is=0;is<nx/2;is++)
                 for(ir=0;ir<nx/2;ir++)
                 {
                    g[is][ir]=(0.0,0.0);
                    swq2.write((char*)&g[is][ir].real(),sizeof(g[is][ir].real()));
                    swq3.write((char*)&g[is][ir].imag(),sizeof(g[is][ir].imag())); 
                 }
               cout<<it<<"frequency slices have been finished!"<<endl;
            }

	  else
            { 
             if((it+1)%5==0)
                cout<<it<<"frequency slices have been finished!"<<endl;
//calculating Green's functions of each source (or equally receiver) in kx-w domain;    
           for(is2=0;is2<nx/2;is2++)
            {
              for(is=0;is<nxb;is++)
                 s[is]=(0.0,0.0);
              s[is2+bou].real()=1.0*ham[it];
              s[is2+bou].imag()=0.0;

              if((in1==NULL)||(out1==NULL))
                cout<<"memory insufficient"<<endl;
              else
                {
                  for(is=0;is<nxb;is++)
                    {
                       in1[is][0]=s[is].real();
                       in1[is][1]=s[is].imag();
                    }
                }

              fftwf_execute(p1);

              for(is=0;is<nxb;is++)
                {
                   s_k[is].real()=out1[is][0];
                   s_k[is].imag()=out1[is][1];
                }
//forward extrapolation

               for(ir=0;ir<nx/2;ir++) 
                {
                 for(is=0;is<nxb;is++)
                   {
                       if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                        {
                           g_k[is].real()=0.0;
                           g_k[is].imag()=0.0;
                        }

                       else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                        {
                          gtmp1.real()=0.0;               //phase-shift extrapolating operator
                          gtmp1.imag()=-mir_dep[ir]*dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                          g_k[is]=s_k[is]*exp(gtmp1);
                        }
                       else
                        {
                           g_k[is].real()=0.0;
                           g_k[is].imag()=0.0;
                        }
                   }
//final Green's functions in kx-w domain for a certain (kx,w) pair;
                 if((in2==NULL)||(out2==NULL))
                 cout<<"memory insufficient"<<endl;
                 else
                   {   
                      for(is=0;is<nxb;is++)
                        {   
                           in2[is][0]=g_k[is].real();
                           in2[is][1]=g_k[is].imag();
                        }   
                   }   
                 fftwf_execute(p2);

                 for(is=0;is<nxb;is++)
                   {   
                      g_tmp[is].real()=out2[is][0]/nxb;
                      g_tmp[is].imag()=out2[is][1]/nxb;
                   }  

                  g[is2][ir]=g_tmp[ir+bou];  
 
                }   

             }//end of is2 loop, .ie. Green's function of all shots for a certain frequency 

            for(is=0;is<nx/2;is++)
                 for(ir=0;ir<nx/2;ir++)
                 {
                    swq2.write((char*)&g[is][ir].real(),sizeof(g[is][ir].real()));
                    swq3.write((char*)&g[is][ir].imag(),sizeof(g[is][ir].imag()));
                 } 

            }//end of else,.ie. a certain frequency
         }//end of all positive frequency

        swq2.close();
        swq3.close();

        cout<<"Positive frequcies prediction Done!"<<endl;

       ifstream swq22;
       swq22.open(fn2,ios::binary);
       if(!swq22)
         {
           cout<<"cannot open "<<fn2<<endl;
           abort();
         }

       ifstream swq33;
       swq33.open(fn3,ios::binary);
       if(!swq33)
         {
           cout<<"cannot open "<<fn3<<endl;
           abort();
         }

       ofstream swq4;
       swq4.open(fn4,ios::binary);
       if(!swq4)
         {
             cout<<"cannot open "<<fn4<<endl;
             abort();
         }

//start to compute the predicticted model of the negative frequency
        cout<<"Negative frequcies prediction Done!"<<endl;
        cout<<"IFFT Starts..."<<endl;
        
        for(is=0;is<nsr/4;is++)
          {
           swq22.seekg(is*4,ios::beg); 
           swq22.read((char*)&gf[0].real(),sizeof(gf[0].real())); 
           swq33.seekg(is*4,ios::beg);  
           swq33.read((char*)&gf[0].imag(),sizeof(gf[0].imag()));

           for(it=1;it<ltt/2+1;it++)
             {
               swq22.seekg((nsr/4-1)*4,ios::cur);
               swq33.seekg((nsr/4-1)*4,ios::cur);
               swq22.read((char*)&gf[it].real(),sizeof(gf[it].real()));
               swq33.read((char*)&gf[it].imag(),sizeof(gf[it].imag()));
             } 
           for(it=ltt/2+1;it<ltt;it++)
               {
                 gf[it].real()=gf[ltt-it].real();
                 gf[it].imag()=-gf[ltt-it].imag();
               } 
              
           if((in3==NULL)||(out3==NULL))
              cout<<"memory insufficient"<<endl;
           else
              {
                for(it=0;it<ltt;it++)
                  {
                     in3[it][0]=gf[it].real();
                     in3[it][1]=gf[it].imag();
                  }
              }

              fftwf_execute(p3);

              for(it=0;it<ltt;it++)
                 gt[it]=out3[it][0]/ltt;

              for(it=0;it<lt;it++)
                swq4.write((char*)&gt[it],sizeof(gt[it]));
   
          }

       
        swq22.close(); 
        swq33.close();
        swq4.close();

        cout<<"disk written."<<endl;

        cout<<"all done!"<<endl;

		
        return 0;
 
   
}






