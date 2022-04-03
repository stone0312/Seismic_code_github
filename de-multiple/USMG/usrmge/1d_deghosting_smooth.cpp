#include "iostream.h"
#include "fstream.h"
#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "stdlib.h"
#include "fftw3.h"
#define pai 3.14159265

using namespace std;

void tripd(float *d, float *e, float *b, int n);
int smooth2 (float *vel, float r1,  float *tripd_d, float *tripd_e, float *tripd_f, int n1);
int smooth2 (float *vel, float r1,  float *tripd_d, float *tripd_e, float *tripd_f, int n1);
int     smooth2fb( float *vel, float r1, int nx);

int     smooth2fb( float *vel, float r1, int nx)
{

        /* scale the smoothing parameter */
        r1 = r1*r1*0.25;

        float   *dd;
        dd=alloc1float(nx);
        float   *ee;
        ee=alloc1float(nx);
        float   *ff;
        ff=alloc1float(nx);

        int     status  = 0 ;
        status  = smooth2(vel, r1, dd, ee, ff,nx);

        free(dd);
        free(ee);
        free(ff);

        return status;
}
void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Input:
d[]     the diagonal of A 
e[]     the superdiagonal of A
b[]     the rhs of Ax=b

Output:
b[]     b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
*****************************************************************************/
{
        int k;
        float temp;

        /* decomposition */
        for(k=1; k<n; ++k){
           temp = e[k-1];
           e[k-1] = temp/d[k-1];
           d[k] -= temp*e[k-1];
        }

        /* substitution */
        for(k=1; k<n; ++k)  b[k] -= e[k-1]*b[k-1];

        b[n-1] /=d[n-1];
        for(k=n-1; k>0; --k)  b[k-1] = b[k-1]/d[k-1] - e[k-1]*b[k];

 }

 int smooth2 (float *vel, float r1, float *tripd_d, float *tripd_e, float *tripd_f, int n1)
 {
        int ix;

        for(ix=0; ix<n1-1; ++ix)
          {
             tripd_d[ix] = 1.0+r1*2.0;
             tripd_e[ix] = -r1;
             tripd_f[ix] = vel[ix];
          }
        tripd_d[0] -= r1;
        tripd_d[n1-1] = 1.0+r1;
        tripd_f[n1-1] = vel[n1-1];
        tripd(tripd_d,tripd_e,tripd_f,n1);
        for(ix=0; ix<n1; ++ix)
           vel[ix] = tripd_f[ix];

        return 0;
 }

int main()
{
   int lt,ltt;
   float dt,rz,v,df;

   lt=4000;
   ltt=lt+400;

   dt=0.5;
   rz=100.0;
   v=1500.0;
 
   df=1000.0/float(lt)/dt;

   float lamda;
   lamda=100;
  
   int it,tmp,ittmp;
 
   float degh_amp_max,ftmp;
 
   float *u;
   u=alloc1float(lt);

   complex<float> *uf;
   uf=alloc1complex(lt);

   complex<float> *ps;
   ps=alloc1complex(lt);

   complex<float> *dg;
   dg=alloc1complex(lt);

   complex<float> *udgf;
   udgf=alloc1complex(lt);

   float *phase;
   phase=alloc1float(lt);
  
   float *psr;
   psr=alloc1float(lt);
  
   float *psi;
   psi=alloc1float(lt);
  
   float *udg;
   udg=alloc1float(lt);

   float *omega;
   omega=alloc1float(lt);
   for(it=0;it<lt/2+1;it++)
     omega[it]=2*pai*it*1000/(dt*lt);
   for(it=lt/2+1;it<lt;it++)
     omega[it]=2*pai*(-1000/(2*dt)+(it-lt/2)*1000/(dt*lt));

  float *degh_amp;
  degh_amp=alloc1float(lt);

  float *amp_orig;
  amp_orig=alloc1float(lt);

  float *filt;
  filt=alloc1float(lt);
  float *filt_tmp;
  filt_tmp=alloc1float(ltt);
  for(it=0;it<ltt;it++)
     filt_tmp[it]=0.0;
  float *amp;
  amp=alloc1float(lt);
  
  fftwf_complex *in1,*out1,*in2,*out2;
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

  fftwf_plan p1,p2;

  p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
  p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

  ifstream swq1;
  swq1.open("./1st_trace_primary_ghost_100m.dat",ios::binary);

  ofstream swq2;
  swq2.open("./orig_amp_1st_trace_100m.dat",ios::binary);

  ofstream swq3;
  swq3.open("./deghost_operator_amp_1st_trace_100m.dat",ios::binary);

  ofstream swq4;
  swq4.open("./deghost_amp_1st_trace_100m_lamda100.dat",ios::binary);

  ofstream swq6;
  swq6.open("./deghost_amp_1st_trace_100m_lamda100_inter.dat",ios::binary);
 
  ofstream swq5;
  swq5.open("./1st_trace_100m_deghosting_lamda100.dat",ios::binary);

  for(it=0;it<lt;it++)
    swq1.read((char*)&u[it],sizeof(u[it]));
  swq1.close();      

  for(it=0;it<lt;it++)
    {
       in1[it][0]=u[it];
       in1[it][1]=0.0;
    }

   fftwf_execute(p1);

   for(it=0;it<lt;it++)
     {
        uf[it].real()=out1[it][0]/lt;
        uf[it].imag()=out1[it][1]/lt;
     }

   for(it=0;it<lt;it++)
     {
       amp_orig[it]=sqrt(pow(uf[it].real(),2)+pow(uf[it].imag(),2));    
       swq2.write((char*)&amp_orig[it],sizeof(amp_orig[it]));
     }
   swq2.close();
/*
   for(it=0;it<lt;it++)
      cout<<it<<"  "<<uf[it].real()<<"   "<<uf[it].imag()<<endl;
   return 0;
*/
   for(it=0;it<lt;it++)
     {
         ps[it].real()=0.0;
         ps[it].imag()=2*rz*omega[it]/v;
 
         ps[it]=exp(ps[it]);

         dg[it].real()=1-ps[it].real();
         dg[it].imag()=ps[it].imag();         

         degh_amp[it]=sqrt(pow(dg[it].real(),2)+pow(dg[it].imag(),2)); 
         swq3.write((char*)&degh_amp[it],sizeof(degh_amp[it]));
     }
    swq3.close();

   degh_amp_max=0.0;
   for(it=0;it<lt;it++)
     {
        if(degh_amp[it]>degh_amp_max)
           degh_amp_max=degh_amp[it];
        else
           degh_amp_max=degh_amp_max;           
     }
   cout<<"maximum amp of deghosting operator is===="<<degh_amp_max<<endl;

    tmp=0; 

    for(it=0;it<lt;it++)
     {  
        if(degh_amp[it]/degh_amp_max>0.1)
          {
            udgf[it]=uf[it]/dg[it];
            psr[it]=udgf[it].real()/sqrt(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag()); 
            psi[it]=udgf[it].imag()/sqrt(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag());
          }
        else
          {
             tmp+=1;

             dg[it].real()+=lamda;
             udgf[it]=uf[it]/dg[it];
             
            psr[it]=udgf[it].real()/(sqrt(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag())+lamda);
            psi[it]=udgf[it].imag()/(sqrt(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag())+lamda);
          }

         filt[it]=sqrt(pow(udgf[it].real(),2)+pow(udgf[it].imag(),2));
         swq4.write((char*)&filt[it],sizeof(filt[it]));
     }
    swq4.close();

    for(it=0;it<lt;it++)
       phase[it]=acos(psr[it]);
    ofstream swq22;
    swq22.open("./phase_deghost.dat",ios::binary);
    for(it=0;it<lt;it++)
      swq22.write((char*)&phase[it],sizeof(phase[it]));
    swq22.close();

/*    
    for(it=0;it<lt;it++)
       {
          if(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag()==0.0)
             {
                psr[it]=0.0;
                psi[it]=0.0;
             }
          else
             {
                psr[it]=udgf[it].real()/sqrt(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag());
                psi[it]=udgf[it].imag()/sqrt(udgf[it].real()*udgf[it].real()+udgf[it].imag()*udgf[it].imag());
             }
       }
*/
/*
   for(it=0;it<lt;it++)
     cout<<it<<"   "<<psr[it]<<"  "<<psi[it]<<endl;
   return 0;
*/

   for(it=200;it<lt+200;it++) 
      filt_tmp[it]=filt[it-200];

   smooth2fb(filt_tmp, 20.0, ltt);
 
   for(it=0;it<lt/2+1;it++)
      amp[it]=filt_tmp[it+200];
   for(it=lt/2+1;it<lt;it++)
      amp[it]=amp[lt-it];

   for(it=0;it<lt;it++)
      swq6.write((char*)&amp[it],sizeof(amp[it]));
   swq6.close();


/*
//smoothing the phase spectrum 
   float psr_tmp[ltt]; 
   float psi_tmp[ltt];
   for(it=0;it<ltt;it++)
      { 
         psr_tmp[it]=0.0;
         psi_tmp[it]=0.0;
      }
   for(it=200;it<lt+200;it++) 
      {
        psr_tmp[it]=psr[it-200];
        psi_tmp[it]=psi[it-200];
      }   
   
   smooth2fb(psr_tmp, 20.0, ltt);
   smooth2fb(psi_tmp, 20.0, ltt);

   for(it=0;it<lt;it++)
     {
       psr[it]=psr_tmp[it+200];
       psi[it]=psi_tmp[it+200];
     }
*/       




/*
   for(it=0;it<lt;it++)
//       cout<<it<<"   "<<amp[it]<<endl;
     cout<<it<<"   "<<amp[it]*psr[it]<<"  "<<amp[it]*psi[it]<<endl;
   return 0;
*/
    
   for(it=0;it<lt;it++)
    {
       in2[it][0]=amp[it]*psr[it];
       in2[it][1]=amp[it]*psi[it];
    }

   fftwf_execute(p2);
   
   for(it=0;it<lt;it++)
     udg[it]=out2[it][0];
   
   for(it=0;it<lt;it++)
     swq5.write((char*)&udg[it],sizeof(udg[it]));
   swq5.close();

   return 0;

}



































