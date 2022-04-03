#include "iostream.h"
#include "fstream.h"
#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "stdlib.h"
#include "fftw3.h"
#define pai 3.14159265

using namespace std;

int main()
{
  char fn1[256],fn2[256];
  int ns,nr,lt,bou,wid,itmin,itmax;
  float dt,sz,rz,dz,dsx,drx,ds,dr,lamda,fmin,fmax,v;
  int it,is,ir,ix,irtmp,tmp;

  ifstream swq;
  swq.open("ghost_forward.par");
  if(!swq)
     {
        cout<<"cannot open ghost_forward.par"<<endl;
        abort();
     }

  swq>>fn1>>fn2>>ns>>nr>>lt>>dt>>sz>>rz>>dz>>dsx>>drx>>ds>>dr>>bou>>wid>>lamda>>fmin>>fmax>>v;
  swq.close();

  cout<<"Fna of input data without Ghosts is===="<<fn1<<endl; 
  cout<<"Fna of output data of Ghosts is===="<<fn2<<endl;
  cout<<"No. of shots is===="<<ns<<endl;
  cout<<"No. of receivers per shot is===="<<nr<<endl;
  cout<<"No. of temperal samplings is===="<<lt<<endl;
  cout<<"Temperal Interval (ms) is===="<<dt<<endl;
  cout<<"Depth of shots and receiver are===="<<sz<<" , "<<rz<<" respectively"<<endl;
  cout<<"Spatial Interval in z direction is===="<<dz<<endl;
  cout<<"Spatial Intervals of shots and receivers are===="<<dsx<<" , "<<drx<<" respectively"<<endl;
  cout<<"Spatial Interval of shots and receivers for calculation are===="<<ds<<" , "<<dr<<" respectively"<<endl;
  cout<<"Boundary width for calculation of Green's Function is===="<<bou<<endl;
  cout<<"Width for aborbing boundary is===="<<wid<<endl;
  cout<<"Minimum and Maximum frequency to be calculated are==== "<<fmin<<" , "<<fmax<<endl;
  cout<<"Velocity of water is===="<<v<<endl;

  int nsr,nxb;

  nsr=ns*nr;
  nxb=nr*int(drx)/int(dr)+2*bou;

  fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in5,*out5,*in6,*out6;
  
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
  out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
  in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
  out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
  in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

  fftwf_plan p1,p2,p3,p4;

  p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
  p2=fftwf_plan_dft_1d(nxb,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);
  p3=fftwf_plan_dft_1d(nxb,in3,out3,FFTW_FORWARD,FFTW_MEASURE);
  p4=fftwf_plan_dft_1d(lt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
    {
       cout<<"cannot open "<<fn1<<endl;
       abort();
    } 

  ofstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
    {
       cout<<"cannot open "<<fn2<<endl;
       abort();
    }

  float **u;
  u=alloc2float(lt,nr);

  float **uu;
  uu=alloc2float(lt,nxb);

  complex<float> **uf;
  uf=alloc2complex(lt,nxb);

  complex<float> **ufk;
  ufk=alloc2complex(lt,nxb);

  complex<float> *ps;
  ps=alloc1complex(nxb);
 
  complex<float> *dg;
  dg=alloc1complex(nxb); 

  complex<float> **ufk1;
  ufk1=alloc2complex(lt,nxb);

  complex<float> **uf1;
  uf1=alloc2complex(lt,nxb);

  float **u1;
  u1=alloc2float(lt,nxb);

  float **u_deg;
  u_deg=alloc2float(lt,nr);

  float *omega;
  omega=alloc1float(lt);

  float *kx;
  kx=alloc1float(nxb);

  float *abso;
  abso=alloc1float(nxb);  
 
//calculate omega and kx
  for(it=0;it<lt/2+1;it++)
     omega[it]=2*pai*it*1000/(dt*lt);
  for(it=lt/2+1;it<lt;it++)
     omega[it]=2*pai*(-1000/(2*dt)+(it-lt/2)*1000/(dt*lt));

  for(ix=0;ix<nxb/2+1;ix++)
     kx[ix]=2*pai*(float)ix/(float)(dr*nxb);    
  for(ix=nxb/2+1;ix<nxb;ix++) 
     kx[ix]=-2*pai/(2*dr)+2*pai*(float)(ix-nxb/2)/(float)(dr*nxb);

  itmin=int(fmin*dt*lt/1000);
  itmax=int(fmax*dt*lt/1000)+1; 
   
  cout<<"Totally "<<itmax-itmin+1<<" frequency slices needed to be calculated..."<<endl; 

  for(ix=0;ix<nxb;ix++)
     abso[ix]=1.0;
/*
  for(ix=0;ix<wid;ix++)
     abso[ix]=sqrt(sin(pai/2*ix/wid));

  for(ix=nxb-wid;ix<nxb;ix++)
     abso[ix]=abso[nxb-ix];
*/
  cout<<"Deghosting Starts..."<<endl;

  for(is=0;is<ns;is++)
    {
       cout<<is+1<<" shot deghosting starts..."<<endl;

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
           swq1.read((char*)&u[ir][it],sizeof(u[ir][it]));

       for(ir=0;ir<nxb;ir++)
         for(it=0;it<lt;it++)
           uu[ir][it]=0.0;
 
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
            {
              irtmp=bou+ir*int(drx)/int(dr);
              uu[irtmp][it]=u[ir][it];
            }


       for(ir=0;ir<nxb;ir++)
          {
             for(it=0;it<lt;it++)
                {
                  in1[it][0]=uu[ir][it];
                  in1[it][1]=0.0;
                }

             fftwf_execute(p1);

             for(it=0;it<lt;it++)
                {
                   uf[ir][it].real()=out1[it][0]/lt;
                   uf[ir][it].imag()=out1[it][1]/lt;
                } 
          }
/*
       for(it=0;it<lt;it++)
         {
           for(ir=0;ir<nxb;ir++)
             {
                 uf[ir][it].real()*=abso[ir];
                 uf[ir][it].imag()*=abso[ir];
             }   
         }
*/

       for(it=0;it<lt;it++)
          {
             for(ix=0;ix<nxb;ix++)
               {
                  in2[ix][0]=uf[ix][it].real();
                  in2[ix][1]=uf[ix][it].imag();
               }
             
             fftwf_execute(p2);

             for(ix=0;ix<nxb;ix++)
               {
                  ufk[ix][it].real()=out2[ix][0]/nxb;
                  ufk[ix][it].imag()=out2[ix][1]/nxb;
               }
          }

       for(it=0;it<lt/2+1;it++)
          {
             if(it<itmin)
                 {
                   for(ix=0;ix<nxb;ix++)
                     { 
                       ufk1[ix][it].real()=0.0;
                       ufk1[ix][it].imag()=0.0;   
                     }                   
                 }

             if(it>itmax+1&&it<lt/2+1)
                 {
                   for(ix=0;ix<nxb;ix++)
                     {
                       ufk1[ix][it].real()=0.0;
                       ufk1[ix][it].imag()=0.0;
                     }
                 }  
             else
                 {
   
                   cout<<it+1<<" frequency slices ... ";
                   tmp=0;
                    
                   for(ix=0;ix<nxb;ix++)
                     {
                       if(pow(omega[it]/v,2)-pow(kx[ix],2)<0)
                          {
                            ufk1[ix][it].real()=0.0;
                            ufk1[ix][it].imag()=0.0;
                          } 
                       else
                          {
                             tmp+=1;

                             ps[ix].real()=0.0;
                             ps[ix].imag()=-2*sqrt(pow(omega[it]/v,2)-pow(kx[ix],2))*rz;
                             ps[ix]=exp(ps[ix]);
                             
                             ufk1[ix][it]=-ufk[ix][it]*ps[ix];
                          }

                     }
                   
                  cout<<"  "<<tmp<<" kx are calculated..."<<endl;                  

                 } 
            }  

        for(it=0;it<lt/2+1;it++)
           {
              for(ix=0;ix<nxb;ix++)
                {
                   in3[ix][0]=ufk1[ix][it].real();
                   in3[ix][1]=ufk1[ix][it].imag();
                }

              fftwf_execute(p3);

              for(ix=0;ix<nxb;ix++)
                {
                   uf1[ix][it].real()=out3[ix][0];
                   uf1[ix][it].imag()=out3[ix][1];
                } 
           }

        for(ix=0;ix<nxb;ix++)
          for(it=lt/2+1;it<lt;it++)
            {
               uf1[ix][it].real()=uf1[ix][lt-it].real();
               uf1[ix][it].imag()=-uf1[ix][lt-it].imag();
            }      
        
        for(ix=0;ix<nxb;ix++)
           {
              for(it=0;it<lt;it++)
                {
                   in4[it][0]=uf1[ix][it].real();
                   in4[it][1]=uf1[ix][it].imag();
                }
              
              fftwf_execute(p4);   

              for(it=0;it<lt;it++)
                 u1[ix][it]=out4[it][0];    

           }  
        
        for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++) 
           {
              irtmp=bou+ir*int(drx)/int(dr);
              u_deg[ir][it]=u1[irtmp][it];    
           }          

        for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++)
            swq2.write((char*)&u_deg[ir][it],sizeof(u_deg[ir][it]));
 
        cout<<is+1<<" shot deghosting done!"<<endl;       
    }

   cout<<"ALL DONE!"<<endl;

   return 0;


}






































