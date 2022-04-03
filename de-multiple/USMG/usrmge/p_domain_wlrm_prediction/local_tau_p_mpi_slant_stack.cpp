using namespace std;
#include <iostream>
#include "math.h"
#include "stdio.h"
#include <fstream>
#include "stdlib.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"
#include "mpi.h"
#define pai 3.14159265

int forward_tau_p(float **u,float **b,float **c,float *p,float *x,int *interc,float **sum1,float **sem,int ntrace,int lt,int np,int N,float theata_max,float dt,float dx,float v,float threhold,int win_time)
{
   int i,j,k,l,ll,ii,jj,kk;

         float **sum3=new float *[np];
         for(i=0;i<np;i++)
                 sum3[i]=new float [lt];

         float **sum4=new float *[np];
         for(i=0;i<np;i++)
                 sum4[i]=new float [lt];

         float *tmp=new float [win_time*2+1];
         float *tmp1=new float [win_time*2+1];

         for(i=0;i<win_time*2+1;i++)
           {
                          tmp[i]=0.0;
                                tmp1[i]=0.0;
           }
         for(i=0;i<np;i++)
                  for(j=0;j<lt;j++)
                        {
                                sum3[i][j]=0.0;
                                sum4[i][j]=0.0;
                        }

   for(i=0;i<np;i++)
      {
            for(j=0;j<lt;j++)
                {
                    for(k=0;k<N+1;k++)
                       interc[k]=int(j+1000*p[i]*x[k]/dt+0.5);

                    for(ii=0;ii<N+1;ii++)
                       {
                           kk=interc[ii];
                           if(kk<0)
                                 sum1[i][j]+=0.0;
                           else if(kk>lt)
                                 sum1[i][j]+=0.0;
                           else
                                 sum1[i][j]+=u[ii][kk];
                        }
                }
            for(j=win_time;j<lt-win_time;j++)
                {
                   k=j-win_time;
                   for(l=k;l<k+win_time*2+1;l++)
                      {
                         for(ll=0;ll<N+1;ll++)
                            interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);

                         for(ii=0;ii<N+1;ii++)
                             {
                                kk=interc[ii];
                                if(kk<0)
                                   {
                                     tmp[l-k]+=0.0;
                                     tmp1[l-k]+=0.0;
                                   }
                                else if(kk>lt)
                                   {
                                        tmp[l-k]+=0.0;
                                        tmp1[l-k]+=0.0;
                                   }
                                else
                                   {
                                        tmp[l-k]+=u[ii][kk];
                                        tmp1[l-k]+=pow(u[ii][kk],2);
                                   }
                             }
                       }
                   for(jj=0;jj<win_time*2+1;jj++)
                      {
                        sum3[i][j]+=tmp[jj];
                        sum4[i][j]+=tmp1[jj];
                      }
                   for(jj=0;jj<win_time*2+1;jj++)
                      {
                        tmp[jj]=0.0;
                        tmp1[jj]=0.0;
                      }
               }
          for(j=0;j<win_time;j++)
              {
                for(l=0;l<win_time*2+1;l++)
                  {
                    for(ll=0;ll<N+1;ll++)
                      interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                    for(ii=0;ii<N+1;ii++)
                      {
                         kk=interc[ii];
                         if(kk<0)
                           {
                             tmp[l]+=0.0;
                             tmp1[l]+=0.0;
                           }
                         else if(kk>lt)
                           {
                              tmp[l]+=0.0;
                              tmp1[l]+=0.0;
                           }
                        else
                           {
                              tmp[l]+=u[ii][kk];
                              tmp1[l]+=pow(u[ii][kk],2);
                           }
                      }
                   }
                  for(jj=0;jj<win_time*2+1;jj++)
                   {
                     sum3[i][j]+=tmp[jj];
                     sum4[i][j]+=tmp1[jj];
                   }
                 for(jj=0;jj<win_time*2+1;jj++)
                   {
                     tmp[jj]=0.0;
                     tmp1[jj]=0.0;
                   }
              }
           for(j=lt-win_time;j<lt;j++)
              {
                 for(l=lt-2*win_time-1;l<lt;l++)
                    {
                        for(ll=0;ll<N+1;ll++)
                          interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                        for(ii=0;ii<N+1;ii++)
                           {
                              kk=interc[ii];
                              if(kk<0)
                                {
                                  tmp[l-lt+2*win_time+1]+=0.0;
                                  tmp1[l-lt+2*win_time+1]+=0.0;
                                }
                              else if(kk>lt)
                                {
                                   tmp[l-lt+2*win_time+1]+=0.0;
                                   tmp1[l-lt+2*win_time+1]+=0.0;
                                }
                              else
                                {
                                   tmp[l-lt+2*win_time+1]+=u[ii][kk];
                                   tmp1[l-lt+2*win_time+1]+=pow(u[ii][kk],2);
                                }
                            }
                     }
                 for(jj=0;jj<win_time*2+1;jj++)
                     {
                        sum3[i][j]+=tmp[jj];
                        sum4[i][j]+=tmp1[jj];                                                                                         
                     }
                for(jj=0;jj<win_time*2+1;jj++)
                     {
                        tmp[jj]=0.0;
                        tmp1[jj]=0.0;
                     }
              }
      }

   for(i=0;i<np;i++)
      for(j=0;j<lt;j++)
        sem[i][j]=pow(sum3[i][j],2)/(sum4[i][j]+0.01);
   for(i=0;i<np;i++)
      for(j=0;j<lt;j++)
        sem[i][j]/=((N+1)*(2*win_time+1));
     for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
         c[i][j]=sem[i][j];

         for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          {
              if(sem[i][j]<threhold)
                 c[i][j]=0.0;
          }

   for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          b[i][j]=sum1[i][j]*c[i][j]/sqrt(N+1);

         return 0;
}

int inverse_tau_p(float **u,float *ufinal,int lt,int np)
{
    int i,j;
    for(i=0;i<lt;i++)
       {
           for(j=0;j<np;j++)
               ufinal[i]+=u[j][i];
       }
    return 0;
}



int main(int argc,char **argv)
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256];
   int nx,lt,np_h,nx_h,ns,nr;
   float dx,dt,v,theta_max,f1,f2,f3,f4,threshold,amp;    

   int is,ir,it,ix,tmp,tmp1,ir1,ir2;

   int npp,myid;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&npp);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   ifstream swq;
   swq.open("local_tau_p_mpi_slant_stack.par");
   if(!swq)
     {
        cout<<"Cannot open local_tau_p_mpi_slant_stack.par"<<endl;
        abort();
     } 

   swq>>fn1>>fn4>>fn6>>fn7>>fn10>>ns>>nr>>dx>>lt>>dt>>v>>np_h>>nx_h>>theta_max>>threshold>>f1>>f2>>f3>>f4;

   swq.close();

   cout<<"fna of input CSG is===="<<fn1<<endl;
   cout<<"fna of output csg local tau-p in time domain is===="<<fn4<<endl;
   cout<<"fna of output real parts of csg local tau-p in frequency domain is===="<<fn6<<endl;
   cout<<"fna of output imaginary parts of csg local tau-p in frequency domain is===="<<fn7<<endl;
   cout<<"fna of output csg inverse local tau-p in time domain is===="<<fn10<<endl;
   cout<<"No. of shots is===="<<ns<<endl;
   cout<<"No. of traces per shot is===="<<nr<<endl;
   cout<<"Reveiver or shot interval is===="<<dx<<"m"<<endl;
   cout<<"Temperal sampling length is===="<<lt<<endl;
   cout<<"Temperal interval is===="<<dt<<"s"<<endl;
   cout<<"Velocity of water is===="<<v<<"m/s"<<endl;
   cout<<"Half of the number of ray parameter is===="<<np_h<<endl;
   cout<<"Half of the number of width of local window is===="<<nx_h<<endl;  
   cout<<"Maximum angle to be calculated is===="<<theta_max<<endl;
   cout<<"Threshold for Local Tau-P Transform is===="<<threshold<<endl;
   cout<<"Range of frequency to be calculated is===="<<f1<<"Hz,  "<<f2<<"Hz,  "<<f3<<"Hz,  "<<f4<<"Hz"<<endl;  

  int ltt,np,nx_l;
  ltt=2*lt;
  np=2*np_h+1;
  nx_l=2*nx_h+1;

  float **us;
  us=alloc2float(ltt,nr);
  float **ugr;
  ugr=alloc2float(ltt,nr);

  float **us1;
  us1=alloc2float(ltt,nr+2*nx_h+1);
  float **ugr1;
  ugr1=alloc2float(ltt,nr+2*nx_h+1);

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
  us_p=alloc2complex(ltt,nr*np);
  complex<float> **ugr_p;
  ugr_p=alloc2complex(ltt,nr*np);

  complex<float> **wlrm_pf;
  wlrm_pf=alloc2complex(ltt,np);   
  
  float **wlrm_pt;
  wlrm_pt=alloc2float(ltt,nr*np);

  float *wlrm;
  wlrm=alloc1float(ltt); 

  float *itp;
  itp=alloc1float(ltt);

  complex<float> **mcg_pf;
  mcg_pf=alloc2complex(ltt,np*nr);

  float **mcg_pt;
  mcg_pt=alloc2float(ltt,nr*np);

  float *mcg;
  mcg=alloc1float(ltt);


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
  mem+=(2*nr*ltt+2*ltt*(nr+2*nx_h)+2*ltt*nx_l+ltt*np*2+ltt*np+3*ltt*nr*np*2+ltt*nr*np+lt*nr+ltt*2)*4.0/1024.0/1024.0;
  
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

  fstream swq888;
  swq888.open("test_io.dat",ios::binary|ios::out);


  for(is=myid;is<ns;is+=npp)
     {//3333

        cout<<"CSG====SHOT  "<<is+1<<"==== Local Tau-p Begins..."<<endl; 
  
        for(ir=0;ir<nr;ir++)
          for(it=0;it<ltt;it++)
              us[ir][it]=0.0;

           swq1.seekg(0,ios::beg);
           for(ir=0;ir<is;ir++)
             swq1.seekg(nr*lt*4,ios::cur);
           for(ir=0;ir<nr;ir++)
            for(it=0;it<lt;it++)
              swq1.read((char*)&(us[ir][it]),sizeof(us[ir][it]));
         
        for(ir=0;ir<nr+2*nx_h+1;ir++)
          for(it=0;it<ltt;it++)        
            us1[ir][it]=0.0;

        for(ir=nx_h;ir<nr+nx_h;ir++)
          for(it=0;it<ltt;it++)
            us1[ir][it]=us[ir-nx_h][it];

        for(ir=0;ir<nr;ir++)
//        for(ir=71;ir<72;ir++)
          {//5555
             tmp=ir+nx_h;

             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                   us2[ix][it]=us1[tmp-nx_h+ix][it];         


             swq888.seekg(0,ios::beg);
             for(ix=0;ix<ir;ix++)
               swq888.seekg(nx_l*ltt*4,ios::cur);

             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                  swq888.write((char*)&(us2[ix][it]),sizeof(us2[ix][it]));

             amp=0.0;
             for(ix=0;ix<nx_l;ix++)
               for(it=0;it<ltt;it++)
                  amp+=us2[ix][it]*us2[ix][it];
  
             amp=amp/(nx_l*ltt);
             amp=sqrt(amp); 
           
             cout<<is<<" , "<<ir<<", "<<amp<<endl;           


             if(amp<0.00001)
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
/*
            fstream swq888;
             swq888.open("test_tp.dat",ios::binary|ios::out);
             if(!swq888)
                {
                   cout<<"cannot open 888"<<endl;
                   abort();
                }

             for(ix=0;ix<np;ix++)
               for(it=0;it<lt;it++)
                  swq888.write((char*)&(tps[ix][it]),sizeof(tps[ix][it]));
             return 0; 
*/

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

            for(it=0;it<ltt;it++)
               itp[it]=0.0;

            for(it=0;it<ltt;it++)
               {
                 for(ix=0;ix<np;ix++)
                    itp[it]+=tps[ix][it];
               }

            swq4.seekg(0,ios::beg);
            swq6.seekg(0,ios::beg);
            swq7.seekg(0,ios::beg);

            for(ir1=0;ir1<is*nr+ir;ir1++)
              {
                swq4.seekg(np*ltt*4,ios::cur);
                swq6.seekg(np*ltt*4,ios::cur);
                swq7.seekg(np*ltt*4,ios::cur);
              }
            for(ix=0;ix<np;ix++)
              for(it=0;it<ltt;it++)
                {
                  swq4.write((char*)&(tps[ix][it]),sizeof(tps[ix][it]));
                  swq6.write((char*)&(ctps[ix][it].real()),sizeof(ctps[ix][it].real()));
                  swq7.write((char*)&(ctps[ix][it].imag()),sizeof(ctps[ix][it].imag()));
                }

            swq10.seekg(0,ios::beg);
            for(ir1=0;ir1<is*nr+ir;ir1++)
              swq10.seekg(lt*4,ios::cur);
            for(it=0;it<lt;it++)
              swq10.write((char*)&(itp[it]),sizeof(itp[it]));


          }//5555
   
 
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

   return 0;

}







































