#include "iostream.h"
#include "math.h"
#include "stdio.h"
#include "fstream.h"
#include "stdlib.h"
#include "complex.h"
#include "alloc.c"
#include "fftw3.h"
#include "mpi.h"
#define pai 3.14159265

extern "C"
{
   void LS_TAU_P__ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void LS_TAU_P_ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void ls_tau_p (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void ls_tau_p_ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

}

using namespace std;

int main(int argc,char **argv)
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256];
   int nx,lt,np_h,nx_h;
   float dx,dt,v,theta_max,f1,f2,f3,f4,threshold,amp;    

   int is,ir,it,ix,tmp,tmp1,ir1,ir2,iw,ip;

   float  coe,phase,w,p,x;
   complex<float> Atmp,sum1;  
 
   int npp,myid;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&npp);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   ifstream swq;
   swq.open("local_tau_p_jump.par");
   if(!swq)
     {
        cout<<"Cannot open local_tau_p_jump.par"<<endl;
        abort();
     } 

   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>fn11>>nx>>dx>>lt>>dt>>v>>np_h>>nx_h>>theta_max>>threshold>>f1>>f2>>f3>>f4;

   swq.close();

   cout<<"fna of input CSG is===="<<fn1<<endl;
   cout<<"fna of input CRG of water-layer Green's Function is===="<<fn2<<endl;
   cout<<"fna of output WLRM to be predicted is===="<<fn3<<endl;
   cout<<"fna of output csg local tau-p in time domain is===="<<fn4<<endl;
   cout<<"fna of output crg local tau-p in time domain is===="<<fn5<<endl;
   cout<<"fna of output real parts of csg local tau-p in frequency domain is===="<<fn6<<endl;
   cout<<"fna of output imaginary parts of csg local tau-p in frequency domain is===="<<fn7<<endl;
   cout<<"fna of output real parts of crg local tau-p in frequency domain is===="<<fn8<<endl;
   cout<<"fna of output imaginary parts of crg local tau-p in frequency domain is===="<<fn9<<endl;
   cout<<"fna of output csg inverse local tau-p in time domain is===="<<fn10<<endl;
   cout<<"fna of output crg inverse local tau-p in time domain is===="<<fn11<<endl;
   cout<<"No. of shots or traces per shot is===="<<nx<<endl;
   cout<<"Reveiver or shot interval is===="<<dx<<"m"<<endl;
   cout<<"Temperal sampling length is===="<<lt<<endl;
   cout<<"Temperal interval is===="<<dt<<"s"<<endl;
   cout<<"Velocity of water is===="<<v<<"m/s"<<endl;
   cout<<"Half of the number of ray parameter is===="<<np_h<<endl;
   cout<<"Half of the number of width of local window is===="<<nx_h<<endl;  
   cout<<"Maximum angle to be calculated is===="<<theta_max<<endl;
   cout<<"Threshold for Local Tau-P Transform is===="<<threshold<<endl;
   cout<<"Range of frequency to be calculated is===="<<f1<<"Hz,  "<<f2<<"Hz,  "<<f3<<"Hz,  "<<f4<<"Hz"<<endl;  

    float *ham;
    ham=alloc1float(2*nx_h+1);
    for(ix=0;ix<2*nx_h+1;ix++)
      ham[ix]=0.5-0.5*cos(2*pai*ix/(2*nx_h));
  
    for(ix=0;ix<2*nx_h+1;ix++)
      cout<<ix<<"  "<<ham[ix]<<endl;
 //   return 0;

  int ltt,np,nx_l;
  ltt=2*lt;
  np=2*np_h+1;
  nx_l=2*nx_h+1;

  int nw1,nw2,nw3,nw4,nw;
  float dw;

  dw=2*pai/ltt/dt;
  nw1=int(2*pai*f1/dw+0.5);
  nw2=int(2*pai*f2/dw+0.5);
  nw3=int(2*pai*f3/dw+0.5);
  nw4=int(2*pai*f4/dw+0.5);

  nw=nw4-nw1+1;

  cout<<"Frequency Slice Number is===="<<nw1<<" , "<<nw2<<" , "<<nw3<<" , "<<nw4<<endl;

  float pmax,dp;
  pmax=sin(theta_max*2*pai/360.0)/v;
  dp=2*pmax/(np-1);

  cout<<"Pmax and Dp is===="<<pmax<<" , "<<dp<<endl;

  float **us;
  us=alloc2float(ltt,nx);
  float **ugr;
  ugr=alloc2float(ltt,nx);

  float **us1;
  us1=alloc2float(ltt,nx+nx_l);

  float **us2;
  us2=alloc2float(ltt,nx_l);
  float **ugr2;
  ugr2=alloc2float(ltt,nx_l);

  complex<float> **ctps;
  ctps=alloc2complex(ltt,np);
  float **tps;
  tps=alloc2float(ltt,np);

  complex<float> **ctps1;
  ctps1=alloc2complex(ltt,np);

  float **tpgr;
  tpgr=alloc2float(ltt,np);

  complex<float> **us_p;
  us_p=alloc2complex(ltt,nx*np);

  float **itp;
  itp=alloc2float(ltt,nx+nx_l);

  float **itp_tmp;
  itp_tmp=alloc2float(ltt,nx_l);
 
  complex<float> **citp;
  citp=alloc2complex(ltt,nx_l);

  complex<float> **AA;
  AA=alloc2complex(nx_l,np);

//define the arrays used in the subroutine
  float **shot_all;
  shot_all=alloc2float(ltt,nx_l); 
  float **semb;
  semb=alloc2float(ltt,np); 
  float *energy;
  energy=alloc1float(np); 
  complex<float> *ctrace;
  ctrace=alloc1complex(ltt); 
  complex<float> **cshot;
  cshot=alloc2complex(ltt,nx_l);
  complex<float> **csemb;
  csemb=alloc2complex(ltt,np);

  complex<float> **A;
  A=alloc2complex(nx_l,np);
  complex<float> **AH;
  AH=alloc2complex(np,nx_l);
  complex<float> **AHA;
  AHA=alloc2complex(np,np);
  complex<float> *f;
  f=alloc1complex(np);
  complex<float> *x1;
  x1=alloc1complex(np);
  complex<float> *x2;
  x2=alloc1complex(np);
  complex<float> *P1;
  P1=alloc1complex(np);
  complex<float> *P2;
  P2=alloc1complex(np);
  complex<float> *R1;
  R1=alloc1complex(np);
  complex<float> *R2;
  R2=alloc1complex(np);
  complex<float> *AP;
  AP=alloc1complex(np);


   fftwf_complex *in1,*out1;
   fftwf_plan p1;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   //p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);
   p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_ESTIMATE);

  float mem=0.0;
  mem+=(2*nx*ltt+2*ltt*(nx+2*nx_h)+2*ltt*nx_l+ltt*np*2+ltt*np+3*ltt*nx*np*2+ltt*nx*np+lt*nx+ltt*2)*4.0/1024.0/1024.0;
  
  cout<<"Totally "<<mem<<"MB needed to be allocated..."<<endl;

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
      {
         cout<<"cannot open "<<fn1<<endl;
         abort();
      }

   fstream swq4;
   swq4.open(fn4,ios::binary|ios::out);
   if(!swq4)
      {
         cout<<"cannot open "<<fn4<<endl;
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
   fstream swq10;
   swq10.open(fn10,ios::binary|ios::out);
   if(!swq10)
      {
         cout<<"cannot open "<<fn10<<endl;
         abort();
      }

  cout<<"CSG Local Tau-p Begins..."<<endl;

  for(is=myid;is<nx;is+=npp)
     {//3333

        cout<<"CSG====SHOT  "<<is+1<<"==== Local Tau-p Begins..."<<endl; 
  
        for(ir=0;ir<nx;ir++)
          for(it=0;it<ltt;it++)
              us[ir][it]=0.0;

           swq1.seekg(0,ios::beg);
           for(ir=0;ir<is;ir++)
             swq1.seekg(nx*lt*4,ios::cur);
           for(ir=0;ir<nx;ir++)
            for(it=0;it<lt;it++)
              swq1.read((char*)&(us[ir][it]),sizeof(us[ir][it]));
         
        for(ir=0;ir<nx+nx_l;ir++)
          for(it=0;it<ltt;it++)        
            us1[ir][it]=0.0;

        for(ir=nx_h;ir<nx;ir++)
          for(it=0;it<ltt;it++)
            us1[ir][it]=us[ir-nx_h][it];

        for(ir=0;ir<nx+nx_l;ir++)
          for(it=0;it<ltt;it++)
            itp[ir][it]=0.0;          

        for(ir=0;ir<nx/(nx_l/2);ir++)
          {//5555
             tmp=ir*nx_l/2;

             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                   us2[ix][it]=us1[tmp+ix][it];         

             cout<<ir<<endl;

             amp=0.0;
             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                  amp+=us2[ix][it]*us2[ix][it];
             amp=sqrt(amp);
             if(amp<0.000001)
                {
                   for(ix=0;ix<np;ix++)
                     for(it=0;it<ltt;it++)
                       {
                          tps[ix][it]=0.0;
                          ctps[ix][it].real()=0.0;
                          ctps[ix][it].imag()=0.0;
                       }
                }
            else
              ls_tau_p_ (us2[0], ctps[0], tps[0], &ltt, &dx, &dt, &v, &np_h, &nx_h, &np, &nx_l, &theta_max, &threshold, &f1, &f2, &f3, &f4);            

            if(ir==0)
             {
               for(ix=0;ix<np;ix++)
                 for(it=0;it<ltt;it++)
                   {
                     tps[ix][it]=0.0;
                     ctps[ix][it].real()=0.0;
                     ctps[ix][it].imag()=0.0;
                   }
             }

            swq6.seekg(0,ios::beg);
            swq7.seekg(0,ios::beg);

            for(ir1=0;ir1<is*nx+ir;ir1++)
              {
                swq6.seekg(np*ltt*4,ios::cur);
                swq7.seekg(np*ltt*4,ios::cur);
              }
            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                {
                  swq6.write((char*)&(ctps[ix][it].real()),sizeof(ctps[ix][it].real()));
                  swq7.write((char*)&(ctps[ix][it].imag()),sizeof(ctps[ix][it].imag()));
                }


            if(is%50==0)
               {
                  swq4.seekg(0,ios::beg);
                  for(ir1=0;ir1<int(is/50)*nx+ir;ir1++)
                     swq4.seekg(np*ltt*4,ios::cur);
                  for(ix=0;ix<np;ix++)
                    for(it=0;it<ltt;it++)
                     swq4.write((char*)&(tps[ix][it]),sizeof(tps[ix][it]));

               }

            cout<<ir*nx_h+1<<" to "<<ir*nx_h+nx_l+1<<"  traces done!"<<endl;

//             back_tau_p_(tps[0], itp_tmp[0],ctps[0], citp_tmp[0], AA[0], &ltt, &dt, &dx,&nw1, &nw2, &nw3, &nw4, &nw, &dp, &dw, &nx_h, &np_h, &nx_l, &np);

            for(it=0;it<lt+1;it++)
              {
                  w=it*dw;
                  
                  for(ip=0;ip<np;ip++)
                   {
                      p=(ip-np_h)*dp;
                      for(ix=0;ix<nx_l;ix++)
                        {
                          x=(ix-nx_h)*dx;
                          phase=w*p*x;
                          coe=1.0/sqrt(nx_l);
                          Atmp.real()=cos(phase/360.0); 
                          Atmp.imag()=sin(phase/360.0); 
                          AA[ix][ip]=coe*Atmp;
                        } 
                   }

                 for(ix=0;ix<nx_l;ix++)
                  {
                     sum1=(0.0,0.0);
                     for(ip=0;ip<np;ip++)
                      sum1+=AA[ix][ip]*ctps[ip][it];
                     citp[ix][it]=sum1; 
                  } 
              }

            for(ix=0;ix<nx_l;ix++)
               for(it=lt+1;it<ltt;it++)
                 {
                    citp[ix][it].real()=citp[ix][ltt-it].real();
                    citp[ix][it].imag()=-citp[ix][ltt-it].imag();
                 }

        for(ix=0;ix<nx_l;ix++)
          {
            for(it=0;it<ltt;it++)
              {
                 in1[it][0]=citp[ix][it].real();
                 in1[it][1]=citp[ix][it].imag();
              }

            fftwf_execute(p1);

            for(it=0;it<ltt;it++)
               itp_tmp[ix][it]=out1[it][0];
          }
            
            for(it=0;it<ltt;it++)
              {
                for(ix=0;ix<nx_l;ix++)
                   itp_tmp[ix][it]*=ham[ix];
              }

            for(ix=ir*nx_h;ix<ir*nx_h+nx_l;ix++)
               for(it=0;it<ltt;it++)
                  itp[ix][it]+=itp_tmp[ix-ir*nx_h][it];

            cout<<ir*nx_h+1<<" to "<<ir*nx_h+nx_l<<"  traces done!"<<endl;


          }//5555

         swq10.seekg(0,ios::beg);
         for(ir1=0;ir1<is;ir1++)
             swq10.seekg(nx*lt*4,ios::cur);

         for(ix=0;ix<nx;ix++)
           for(it=0;it<lt;it++)
             swq10.write((char*)&(itp[ix+nx_h][it]),sizeof(itp[ix+nx_h][it]));
 
       cout<<"====SHOT  "<<is+1<<"====Prediction Done!"<<endl; 

     }//3333 end of source   

     MPI_Barrier(MPI_COMM_WORLD);

    cout<<"CSG Local Tau-p Done!"<<endl;
     
    MPI_Finalize();

    swq1.close();
    swq4.close();
    swq6.close();
    swq7.close();
    swq10.close();

    cout<<"CSG Local Tau-p Done!"<<endl;
    cout<<"Local Tau-P Transform Done!"<<endl;

   cout<<"ALL DONE!"<<endl;

   return 0;

}








































