using namespace std;
#include <iostream>
#include "math.h"
#include "stdio.h"
#include <fstream>
#include "stdlib.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

extern "C"
{
   void LS_TAU_P__ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void LS_TAU_P_ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void ls_tau_p (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void ls_tau_p_ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

}

using namespace std;

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256];
   int nx,lt,np_h,nx_h;
   float dx,dt,v,theta_max,f1,f2,f3,f4,threshold;    

   int is,ir,it,ix,tmp,tmp1,ir1,ir2;

   ifstream swq;
   swq.open("local_tau_p.par");
   if(!swq)
     {
        cout<<"Cannot open local_tau_p.par"<<endl;
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

//   return 0;

  int ltt,np,nx_l;
  ltt=2*lt;
  np=2*np_h+1;
  nx_l=2*nx_h+1;

  float **us;
  us=alloc2float(ltt,nx);
  float **ugr;
  ugr=alloc2float(ltt,nx);

  float **us1;
  us1=alloc2float(ltt,nx+2*nx_h+1);
  float **ugr1;
  ugr1=alloc2float(ltt,nx+2*nx_h+1);

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

  complex<float> **ctpgr;
  ctpgr=alloc2complex(ltt,np);
  float **tpgr;
  tpgr=alloc2float(ltt,np);

  complex<float> **us_p;
  us_p=alloc2complex(ltt,nx*np);
  complex<float> **ugr_p;
  ugr_p=alloc2complex(ltt,nx*np);

  complex<float> **wlrm_pf;
  wlrm_pf=alloc2complex(ltt,np);   
  
  float **wlrm_pt;
  wlrm_pt=alloc2float(ltt,nx*np);

  float *wlrm;
  wlrm=alloc1float(ltt); 

  float *itp;
  itp=alloc1float(ltt);

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
   p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);

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

   ifstream swq2;
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

   ofstream swq4;
   swq4.open(fn4,ios::binary);
   if(!swq4)
      {
         cout<<"cannot open "<<fn4<<endl;
         abort();
      }

   ofstream swq5;
   swq5.open(fn5,ios::binary);
   if(!swq5)
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

   ofstream swq8;
   swq8.open(fn8,ios::binary);
   if(!swq8)
      {
         cout<<"cannot open "<<fn8<<endl;
         abort();
      }
  
   ofstream swq9;
   swq9.open(fn9,ios::binary);
   if(!swq9)
      {
         cout<<"cannot open "<<fn3<<endl;
         abort();
      }


   ofstream swq10;
   swq10.open(fn10,ios::binary);
   if(!swq10)
      {
         cout<<"cannot open "<<fn10<<endl;
         abort();
      }
   
   ofstream swq11;
   swq11.open(fn11,ios::binary);
   if(!swq11)
      {
         cout<<"cannot open "<<fn11<<endl;
         abort();
      }


  cout<<"CSG Local Tau-p Begins..."<<endl;

  for(is=0;is<nx;is++)
     {//3333

        cout<<"CSG====SHOT  "<<is+1<<"==== Local Tau-p Begins..."<<endl; 
  
        for(ir=0;ir<nx;ir++)
          for(it=0;it<ltt;it++)
              us[ir][it]=0.0;

        for(ir=0;ir<nx;ir++)
          for(it=0;it<lt;it++)
            swq1.read((char*)&(us[ir][it]),sizeof(us[ir][it])); 

        for(ir=0;ir<nx+2*nx_h+1;ir++)
          for(it=0;it<ltt;it++)        
            us1[ir][it]=0.0;

        for(ir=nx_h;ir<nx+nx_h;ir++)
          for(it=0;it<ltt;it++)
            us1[ir][it]=us[ir-nx_h][it];


        for(ir=0;ir<nx;ir++)
          {
             tmp=ir+nx_h;

             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                   us2[ix][it]=us1[tmp-nx_h+ix][it];         

            ls_tau_p_ (us2[0], ctps[0], tps[0], &ltt, &dx, &dt, &v, &np_h, &nx_h, &np, &nx_l, &theta_max, &threshold, &f1, &f2, &f3, &f4);            
             cout<<ir<<"  2222"<<endl;

            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                swq4.write((char*)&(tps[ix][it]),sizeof(tps[ix][it]));

            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                swq6.write((char*)&(ctps[ix][it].real()),sizeof(ctps[ix][it].real()));

            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                swq7.write((char*)&(ctps[ix][it].imag()),sizeof(ctps[ix][it].imag()));


            for(it=0;it<ltt;it++)
               itp[it]=0.0;

            for(it=0;it<ltt;it++)
               {
                 for(ix=0;ix<np;ix++)
                    itp[it]+=tps[ix][it];
               } 

            for(it=0;it<lt;it++)
               swq10.write((char*)&(itp[it]),sizeof(itp[it]));


             for(ix=0;ix<np;ix++)
               for(it=0;it<ltt;it++)
                  {
                     ctps1[ix][it].real()=ctps[np-1-ix][it].real();
                     ctps1[ix][it].imag()=ctps[np-1-ix][it].imag();
                  }
       
             for(ix=ir*np;ix<(ir+1)*np;ix++)
                for(it=0;it<ltt;it++)
                   {
                      us_p[ix][it].real()=ctps1[ix-ir*np][it].real();
                      us_p[ix][it].imag()=ctps1[ix-ir*np][it].imag();
                   }

          }
    
       cout<<"====SHOT  "<<is+1<<"====Prediction Done!"<<endl; 
      
     }//3333 end of source   

    cout<<"CSG Local Tau-p Done!"<<endl;

    cout<<"CRG Local Tau-p Begins..."<<endl;

  for(is=0;is<nx;is++)
     {//3333

        cout<<"CRG====Receiver  "<<is+1<<"==== Local Tau-p Begins..."<<endl;

        for(ir=0;ir<nx;ir++)
          for(it=0;it<ltt;it++)
              us[ir][it]=0.0;

        for(ir=0;ir<nx;ir++)
          for(it=0;it<lt;it++)
            swq2.read((char*)&(us[ir][it]),sizeof(us[ir][it]));

        for(ir=0;ir<nx+2*nx_h+1;ir++)
          for(it=0;it<ltt;it++)
            us1[ir][it]=0.0;

        for(ir=nx_h;ir<nx+nx_h;ir++)
          for(it=0;it<ltt;it++)
            us1[ir][it]=us[ir-nx_h][it];

        for(ir=0;ir<nx;ir++)
          {
             tmp=ir+nx_h;

             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                   us2[ix][it]=us1[tmp-nx_h+ix][it];

            ls_tau_p_ (us2[0], ctps[0], tps[0], &ltt, &dx, &dt, &v, &np_h, &nx_h, &np, &nx_l, &theta_max, &threshold, &f1, &f2, &f3, &f4);     
             cout<<ir<<"  2222"<<endl;

            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                swq5.write((char*)&(tps[ix][it]),sizeof(tps[ix][it]));

            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                swq8.write((char*)&(ctps[ix][it].real()),sizeof(ctps[ix][it].real()));

            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                swq9.write((char*)&(ctps[ix][it].imag()),sizeof(ctps[ix][it].imag()));

            for(it=0;it<ltt;it++)
               itp[it]=0.0;

            for(it=0;it<ltt;it++)
               {
                 for(ix=0;ix<np;ix++)
                    itp[it]+=tps[ix][it];
               }

            for(it=0;it<lt;it++)
               swq11.write((char*)&(itp[it]),sizeof(itp[it]));

          }

       cout<<"CRG====RECEIVER  "<<is+1<<"====Local Tau-p  Done!"<<endl;

     }//3333 end of source   

    swq1.close();
    swq2.close();
    swq4.close();
    swq5.close();
    swq6.close();
    swq7.close();
    swq8.close();
    swq9.close();
//    swq10.close();
//    swq11.close();

    cout<<"CSG Local Tau-p Done!"<<endl;
    cout<<"Local Tau-P Transform Done!"<<endl;


    cout<<"WLRM Prediction in Tau-P Domain Begins..."<<endl;

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
   swq88.open(fn8,ios::binary);
   if(!swq88)
      {
         cout<<"cannot open "<<fn8<<endl;
         abort();
      }

   ifstream swq99;
   swq99.open(fn9,ios::binary);
   if(!swq99)
      {
         cout<<"cannot open "<<fn9<<endl;
         abort();
      }


   for(is=0;is<nx;is++)
     {//3333
        cout<<"ISHOT===="<<is+1<<"  WLRM Prediction Begins..."<<endl;

        for(ir=0;ir<np*nx;ir++)
           for(it=0;it<ltt;it++)
              {
                 swq66.read((char*)&(us_p[ir][it].real()),sizeof(us_p[ir][it].real()));
                 swq77.read((char*)&(us_p[ir][it].imag()),sizeof(us_p[ir][it].imag()));
              }       
//transpose the local tau-p of CSG 

        for(ir=0;ir<nx;ir++)
           {
              for(ir1=ir*np;ir1<ir*np+np;ir1++)
                 for(it=0;it<ltt;it++)
                       ctps[ir1-ir*np][it]=us_p[ir1][it];
  
              for(ir1=0;ir1<np;ir1++)
                 for(it=0;it<ltt;it++)
                     ctps1[ir1][it]=ctps[np-1-ir1][it];

              for(ir1=ir*np;ir1<ir*np+np;ir1++)
                 for(it=0;it<ltt;it++)
                   us_p[ir1][it]=ctps1[ir1-ir*np][it];
           }          
      
      for(ir=0;ir<nx;ir++)
       {//2222 
        for(ir1=0;ir1<np*nx;ir1++)
           for(it=0;it<ltt;it++)
              {
                 swq88.read((char*)&(ugr_p[ir1][it].real()),sizeof(ugr_p[ir1][it].real()));
                 swq99.read((char*)&(ugr_p[ir1][it].imag()),sizeof(ugr_p[ir1][it].imag()));
              } 
         
        for(ir1=0;ir1<np;ir1++)
           for(it=0;it<ltt;it++)
             {
                wlrm_pf[ir1][it].real()=0.0;
                wlrm_pf[ir1][it].imag()=0.0;
             }
        
        for(ir1=0;ir1<np;ir1++)
           {
             for(ir2=0;ir2<nx;ir2++)
              {
                for(it=0;it<ltt;it++)
                    wlrm_pf[ir1][it]+=us_p[ir2*np+ir1][it]*ugr_p[ir2*np+ir1][it];
              } 
           }
//inverse fft
        for(ir1=0;ir1<np;ir1++)
         {
            for(it=0;it<ltt;it++)
              {
                 in1[it][0]=wlrm_pf[ir1][it].real();
                 in1[it][1]=wlrm_pf[ir1][it].imag();
              }

            fftwf_execute(p1);

            for(it=0;it<ltt;it++)
              wlrm_pt[ir1][it]=out1[it][0];
         }
//inverse tau-p transform to get the (is,ir) WLRM Prediction
        for(it=0;it<ltt;it++)
           wlrm[it]=0.0;

         for(it=0;it<ltt;it++)
           {
              for(ix=0;ix<np;ix++)
                wlrm[it]+=wlrm_pt[ix][it];
           }

            for(it=0;it<lt;it++)
               swq3.write((char*)&(wlrm[it]),sizeof(wlrm[it]));

       }//2222

       cout<<"ISHOT===="<<is+1<<"  WLRM Prediction Done!"<<endl;

     } //3333  

   cout<<"ALL DONE!"<<endl;

   return 0;

}








































