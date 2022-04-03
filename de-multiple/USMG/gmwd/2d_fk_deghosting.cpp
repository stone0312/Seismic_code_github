#include "iostream.h"
#include "fstream.h"
#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "stdlib.h"
#include "fftw3.h"
#define pai 3.14159265

void tripd(float *d, float *e, float *b, int n);
int smooth2(float *vel, float r1, float *tripd_d, float *tripd_e, float *tripd_f, int n1);
int smooth2(float *vel, float r1, float *tripd_d, float *tripd_e, float *tripd_f, int n1);
int smooth2fb(float *vel, float r1, int nx);

int smooth2fb(float *vel, float r1, int nx)
{
   r1=r1*r1*0.25;
   
   float *dd;
   dd=alloc1float(nx);
   
   float *ee;
   ee=alloc1float(nx);
   
   float *ff;
   ff=alloc1float(nx);

   int status =0;
   status= smooth2(vel,r1,dd,ee,ff,nx);

   free(dd);
   free(ee);
   free(ff);
  
   return status;
}

void tripd(float *d, float *e, float *b, int n)
{
   int k;
   float temp;
   for(k=1;k<n;k++) 
    {
      temp=e[k-1];
      e[k-1]=temp/d[k-1];
      d[k]-=temp*e[k-1];
    }

   for(k=1;k<n;k++)
     b[k]-=e[k-1]*b[k-1];
   b[n-1]/=d[n-1];

   for(k=n-1;k>0;k--)
     b[k-1]=b[k-1]/d[k-1]-e[k-1]*b[k];

}

int smooth2(float *vel, float r1, float *tripd_d, float *tripd_e, float *tripd_f, int n1)
{
   int ix;
   for(ix=0;ix<n1-1;ix++)
     {
       tripd_d[ix]=1.0+r1*2.0;
       tripd_e[ix]=-r1;
       tripd_f[ix]=vel[ix];
     }

   tripd_d[0]=r1;
   tripd_d[n1-1]=1.0+r1;
   tripd_f[n1-1]=vel[n1-1];
   tripd(tripd_d,tripd_e,tripd_f,n1);
   for(ix=0;ix<n1;ix++)
      vel[ix]=tripd_f[ix];

   return 0;
}

using namespace std;

int main()
{
  char fn1[256],fn2[256];
  int ns,nr,lt,itmin,itmax,ltt;
  float dt,sz,rz,dz,dsx,drx,ds,dr,lamda,fmin,fmax,v,kz,gh_amp_max,degh_amp_max;
  int it,is,ir,ix,irtmp,tmp;

  ifstream swq;
  swq.open("2d_fk_deghosting.par");
  if(!swq)
     {
        cout<<"cannot open 2d_fk_deghosting.par"<<endl;
        abort();
     }

  swq>>fn1>>fn2>>ns>>nr>>lt>>dt>>sz>>rz>>dz>>dsx>>drx>>ds>>dr>>lamda>>fmin>>fmax>>v;
  swq.close();

  ltt=lt+400;

  cout<<"Fna of input data with Ghosts is===="<<fn1<<endl; 
  cout<<"Fna of output data without Ghosts is===="<<fn2<<endl;
  cout<<"No. of shots is===="<<ns<<endl;
  cout<<"No. of receivers per shot is===="<<nr<<endl;
  cout<<"No. of temperal samplings is===="<<lt<<endl;
  cout<<"Temperal Interval (ms) is===="<<dt<<endl;
  cout<<"Depth of shots and receiver are===="<<sz<<" , "<<rz<<" respectively"<<endl;
  cout<<"Spatial Interval in z direction is===="<<dz<<endl;
  cout<<"Spatial Intervals of shots and receivers are===="<<dsx<<" , "<<drx<<" respectively"<<endl;
  cout<<"Spatial Interval of shots and receivers for calculation are===="<<ds<<" , "<<dr<<" respectively"<<endl;
  cout<<"Minimum and Maximum frequency to be calculated are==== "<<fmin<<" , "<<fmax<<endl;
  cout<<"Velocity of water is===="<<v<<endl;

  fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in5,*out5,*in6,*out6;
  
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

  fftwf_plan p1,p2,p3,p4;

  p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
  p2=fftwf_plan_dft_1d(nr,in2,out2,FFTW_FORWARD,FFTW_MEASURE);
  p3=fftwf_plan_dft_1d(nr,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
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

 ofstream swq3;
 swq3.open("gh_amp.dat",ios::binary);

 ofstream swq4;
 swq4.open("orig_amp.dat",ios::binary);

 ofstream swq6;
 swq6.open("degh_amp.dat",ios::binary);

 ofstream swq7;
 swq7.open("degh_amp_smooth.dat",ios::binary);

  float **u;
  u=alloc2float(lt,nr);

  complex<float> **uf;
  uf=alloc2complex(lt,nr);

  complex<float> **ufk;
  ufk=alloc2complex(lt,nr);

  complex<float> *ps;
  ps=alloc1complex(lt);
 
  complex<float> **dg;
  dg=alloc2complex(lt,nr); 

  complex<float> **ufk1;
  ufk1=alloc2complex(lt,nr);

  complex<float> **uf1;
  uf1=alloc2complex(lt,nr);

  float **u1;
  u1=alloc2float(lt,nr);

  float **gh_amp;
  gh_amp=alloc2float(lt,nr);
  float **degh_amp;
  degh_amp=alloc2float(lt,nr);
  float **psr;
  psr=alloc2float(lt,nr);
  float **psi;
  psi=alloc2float(lt,nr);
  float *degh_amp_tmp;
  degh_amp_tmp=alloc1float(ltt);

  float **amp_orig;
  amp_orig=alloc2float(lt,nr);

  float **filt;
  filt=alloc2float(lt,nr);

  float *omega;
  omega=alloc1float(lt);

  float *kx;
  kx=alloc1float(nr);

//calculate omega and kx
  for(it=0;it<lt/2+1;it++)
     omega[it]=2*pai*it*1000/(dt*lt);
  for(it=lt/2+1;it<lt;it++)
     omega[it]=2*pai*(-1000/(2*dt)+(it-lt/2)*1000/(dt*lt));

  for(ix=0;ix<nr/2+1;ix++)
     kx[ix]=2*pai*(float)ix/(float)(dr*nr);    
  for(ix=nr/2+1;ix<nr;ix++) 
     kx[ix]=-2*pai/(2*dr)+2*pai*(float)(ix-nr/2)/(float)(dr*nr);

  itmin=int(fmin*dt*lt/1000);
  itmax=int(fmax*dt*lt/1000)+1; 
   
  cout<<"Totally "<<itmax-itmin+1<<" frequency slices needed to be calculated..."<<endl; 

  cout<<"Deghosting Starts..."<<endl;

  for(is=0;is<ns;is++)
    {
       cout<<is+1<<" shot deghosting starts..."<<endl;

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
           swq1.read((char*)&u[ir][it],sizeof(u[ir][it]));

       for(ir=0;ir<nr;ir++)
          {
             for(it=0;it<lt;it++)
                {
                  in1[it][0]=u[ir][it];
                  in1[it][1]=0.0;
                }

             fftwf_execute(p1);

             for(it=0;it<lt;it++)
                {
                   uf[ir][it].real()=out1[it][0]/lt;
                   uf[ir][it].imag()=out1[it][1]/lt;
                } 
          }

       for(it=0;it<lt;it++)
          {
             for(ix=0;ix<nr;ix++)
               {
                  in2[ix][0]=uf[ix][it].real();
                  in2[ix][1]=uf[ix][it].imag();
               }
             
             fftwf_execute(p2);

             for(ix=0;ix<nr;ix++)
               {
                  ufk[ix][it].real()=out2[ix][0]/nr;
                  ufk[ix][it].imag()=out2[ix][1]/nr;
               }
          }

       for(ix=0;ix<nr;ix++)
          for(it=0;it<lt;it++)
            {
              amp_orig[ix][it]=sqrt(pow(ufk[ix][it].real(),2)+pow(ufk[ix][it].imag(),2));
              swq4.write((char*)&amp_orig[ix][it],sizeof(amp_orig[ix][it])); 
            } 

       for(ix=0;ix<nr;ix++)
          {
            if((ix+1)%10==0)  
              cout<<ix+1<<" traces done... "<<endl;
              
            for(it=0;it<lt;it++)
              {
                 if(pow(omega[it]/v,2)-pow(kx[ix],2)<0)
                   {
                      gh_amp[ix][it]=0.0;
                      swq3.write((char*)&gh_amp[ix][it],sizeof(gh_amp[ix][it]));
                           
                      ufk[ix][it].real()=0.0;
                      ufk[ix][it].imag()=0.0;
                   } 
                  else
                   {
                      ps[it].real()=0.0;
                      ps[it].imag()=2*sqrt(pow(omega[it]/v,2)-pow(kx[ix],2))*rz;
                      ps[it]=exp(ps[it]);

                      dg[ix][it].real()=1.0-ps[it].real();
                      dg[ix][it].imag()=ps[it].imag();
                             
                      gh_amp[ix][it]=sqrt(pow(dg[ix][it].real(),2)+pow(dg[ix][it].imag(),2));
                      swq3.write((char*)&gh_amp[ix][it],sizeof(gh_amp[ix][it]));
                          
                   }
                }
            }  
        swq3.close();

        for(ix=0;ix<nr;ix++)
          {
             for(it=0;it<ltt;it++)
                degh_amp_tmp[it]=0.0;

             gh_amp_max=0.0;
             for(it=0;it<lt;it++)
                {
                   if(gh_amp[ix][it]>gh_amp_max)
                      gh_amp_max=gh_amp[ix][it];
                   else
                      gh_amp_max=gh_amp_max;
                }
//             cout<<ix<<" maximum ghdg_amp is===="<<gh_amp_max<<endl;
//             cout<<ix<<" deghosting real and imag parts are===="<<endl;
           
             for(it=0;it<lt;it++)
                {
//                   cout<<it<<"  "<<gh_amp[ix][it]/gh_amp_max<<endl;
                   if(gh_amp[ix][it]/gh_amp_max>0.1)
                      {
                         ufk1[ix][it]=ufk[ix][it]/dg[ix][it];
                         degh_amp[ix][it]=sqrt(ufk1[ix][it].real()*ufk1[ix][it].real()+ufk1[ix][it].imag()*ufk1[ix][it].imag());
                      }
                   else
                      {
                         dg[ix][it].real()+=lamda;

                         ufk1[ix][it]=ufk[ix][it]/dg[ix][it];
                         degh_amp[ix][it]=sqrt(ufk1[ix][it].real()*ufk1[ix][it].real()+ufk1[ix][it].imag()*ufk1[ix][it].imag());
                      }

                } 
             for(it=0;it<lt;it++)
                swq6.write((char*)&degh_amp[ix][it],sizeof(degh_amp[ix][it]));              
             
             for(it=200;it<lt+200;it++)
               degh_amp_tmp[it]=degh_amp[ix][it-200];
             smooth2fb(degh_amp_tmp,20.0,ltt);
             for(it=0;it<lt;it++)
                degh_amp[ix][it]=degh_amp_tmp[it+200];
             for(it=0;it<lt;it++)
                swq7.write((char*)&degh_amp[ix][it],sizeof(degh_amp[ix][it]));              

             degh_amp_max=0.0;
             for(it=0;it<lt;it++)
                {
                   if(degh_amp[ix][it]>degh_amp_max)
                      degh_amp_max=degh_amp[ix][it];
                   else
                      degh_amp_max=degh_amp_max;
                }
             for(it=0;it<lt;it++)
                {
                   if(degh_amp[ix][it]/degh_amp_max>0.1)
                      {
                         psr[ix][it]=ufk1[ix][it].real()/degh_amp[ix][it];
                         psi[ix][it]=ufk1[ix][it].imag()/degh_amp[ix][it];
                      }
                   else
                     {
                         psr[ix][it]=ufk1[ix][it].real()/(degh_amp[ix][it]+lamda);
                         psi[ix][it]=ufk1[ix][it].imag()/(degh_amp[ix][it]+lamda);
                     }
                }
/*
             for(it=200;it<lt+200;it++)
               degh_amp_tmp[it]=degh_amp[ix][it-200];

             smooth2fb(degh_amp_tmp,20.0,ltt);

             for(it=0;it<lt;it++)
                degh_amp[ix][it]=degh_amp_tmp[it+200];

             for(it=0;it<lt;it++)
                swq7.write((char*)&degh_amp[ix][it],sizeof(degh_amp[ix][it]));              
*/
          }

          for(ix=0;ix<nr;ix++)
            for(it=0;it<lt;it++)
               {
                  ufk1[ix][it].real()=degh_amp[ix][it]*psr[ix][it];
                  ufk1[ix][it].imag()=degh_amp[ix][it]*psi[ix][it];
               }   

          for(it=0;it<lt;it++)
           {
              for(ix=0;ix<nr;ix++)
                {
                  in3[ix][0]=ufk1[ix][it].real();
                  in3[ix][1]=ufk1[ix][it].imag();
                }

              fftwf_execute(p3);

              for(ix=0;ix<nr;ix++)
                {
                   uf1[ix][it].real()=out3[ix][0];
                   uf1[ix][it].imag()=out3[ix][1];
                } 
           }
        
        for(ix=0;ix<nr;ix++)
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
            swq2.write((char*)&u1[ir][it],sizeof(u1[ir][it]));
 
        cout<<is+1<<" shot deghosting done!"<<endl;       
    }

   cout<<"ALL DONE!"<<endl;

   return 0;


}






































