#include <iostream>
#include "math.h"
#include "stdio.h"
#include <fstream>
#include "stdlib.h"
#include <complex>

using namespace std;

#include "alloc.c"
#include "fftw3.h"
#include "mpi.h"
#define pai 3.14159265


//================================================================================================================//
//                                               Forward_LS_tp                                                    //
//================================================================================================================//

extern "C"
{
   void LS_TAU_P__ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void LS_TAU_P_ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void ls_tau_p (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

   void ls_tau_p_ (float *shot_all, complex<float> *ctp, float *tp, int *lt, float *dx,float *dt, float *vmin, int *np_h, int *nx_h, int *np, int *nx_l, float *theta_max, float *threshold, float *f1, float *f2, float *f3, float *f4);

}


//================================================================================================================//
//                                                  Inverse_tp                                                    //
//================================================================================================================//

int inverse_tau_p(float **u,float *ufinal,int lt,int np)
{
    int i,j;
    for(i=0;i<lt;i++)
      ufinal[i]=0.0;
    for(i=0;i<lt;i++)
       {
           for(j=0;j<np;j++)
               ufinal[i]+=u[j][i];
       }
    return 0;
}


//================================================================================================================//
//                                                  Cal_Semblance                                                 //
//================================================================================================================//

int cal_semblance(float **u,float **b,float **cc, float **c,float *p,float *x,int *interc,float **sum1,float **sem,int lt,int np,int N,float theata_max,float theata_min,float dt,float dx,float v,float threhold,int win_time)
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
//               interc[k]=int(j+1000*p[i]*x[k]/dt+0.5);
               interc[k]=int(j+p[i]*x[k]/dt+0.5);

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
//                    interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                    interc[ll]=int(j+p[i]*x[ll]/dt+0.5);

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
//                      interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                      interc[ll]=int(j+p[i]*x[ll]/dt+0.5);
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
//                          interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                          interc[ll]=int(j+p[i]*x[ll]/dt+0.5);
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
//              if(sem[i][j]<threhold)
              if(sem[i][j]<0.0)
                 c[i][j]=0.0;
          }   

    for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
         cc[i][j]=c[np-i-1][j];


   for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          b[i][j]=sum1[i][j]*cc[i][j]/sqrt(N+1);
    
   return 0;  
}


//================================================================================================================//
//                                           Butterworth-type Filter                                              //
//================================================================================================================//

int Butterworth_filter(float **Sp, float **Sd, float **Sm, int np, int lt, int n, float  epi, float lamda) 
{

   int ip, it;
   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
         if(Sd[ip][it]<0.0001)
            Sp[ip][it]=0.0;
         else
            Sp[ip][it]=1/sqrt(1+pow(Sm[ip][it]/(epi*Sd[ip][it]+lamda),n));
       }
   return 0;
}

//================================================================================================================//               
//                                              Cal L2-Norm Error                                                 //
//================================================================================================================//

int l2norm_err(float l2norm, float * R, int lt)
{
   
   int it;
   float l2norm1=0.0;
   for(it=0;it<lt;it++)
     l2norm1+=R[it]*R[it];

   l2norm=l2norm1;

   return 0;
}


//================================================================================================================//
//                                   MCA-Based Primary and Multiple Separation                                    //
//================================================================================================================//

int MCA_Based_PM_Separation(float *u, float *p, float *m, float *r, float *p_final, float *m_final, float **ud, float **um, float **up, float **Sd, float **Sm, float **Sp, float     **Sd1, float **Sm1, float **udtp, float **umtp, float **uptp, float **r1, float **r2, float **umtp1, complex<float> **cudtp, complex<float> **cumtp, int lt, float dx, float dt, float v, int  np_h, int nx_h, int np, int nx_l, float theta_max, float theta_min, float threshold, float f1, float f2, float f3, float f4, float *pr, float *x, int *interc, float **sum1, int win_time, int it_max, float err, float l2err, int n, float epi, float lamda, float **tp_tmp1, float **tp_tmp2, float **sem)
{
  int ix,ip,it, ite_no;
  ite_no=1;
  float amp1,amp2;
  amp1=0.0;
  amp2=0.0; 

  int nx_l1;
  nx_l1=nx_l-1;

//initializatioin
 
  cout<<"Initializatioin ..."<<endl;

//ls tau-p transform of the original data and the predicted multiple

  for(ix=0;ix<nx_l;ix++)
    for(it=0;it<lt;it++)
      amp1+=ud[ix][it]*ud[ix][it]; 

  amp1=amp1/(nx_l*lt);
  amp1=sqrt(amp1);

  for(ix=0;ix<nx_l;ix++)
    for(it=0;it<lt;it++)
      amp2+=um[ix][it]*um[ix][it];             
  
  amp2=amp2/(nx_l*lt);
  amp2=sqrt(amp2);

  if(amp1<0.001)
    {
       for(ix=0;ix<np;ix++)
         for(it=0;it<lt;it++)
           udtp[ix][it]=0.0;

       cout<<"No need to LS-TP of original shot gather ..."<<endl;


    }
  else 
    ls_tau_p_ (ud[0], cudtp[0], udtp[0], &lt, &dx, &dt, &v, &np_h, &nx_h, &np, &nx_l, &theta_max, &threshold, &f1, &f2, &f3, &f4);
/*
   ofstream swq5;
   swq5.open("local_lstp_original_shot_gather.dat",ios::binary);
   if(!swq5)
      {
         cout<<"cannot open local_lstp_original_shot_gather.dat"<<endl;
         abort();
      }
 
   for(ix=0;ix<np;ix++)
     for(it=0;it<lt;it++) 
         swq5.write((char*)&udtp[ix][it],sizeof(udtp[ix][it]));        
 
   swq5.close();
*/ 
  if(amp2<0.001)
    {
       for(ix=0;ix<np;ix++)
         for(it=0;it<lt;it++)
           umtp[ix][it]=0.0;

       cout<<"No need to LS-TP of predicited multiple ..."<<endl;
    }
  else
    ls_tau_p_ (um[0], cumtp[0], umtp[0], &lt, &dx, &dt, &v, &np_h, &nx_h, &np, &nx_l, &theta_max, &threshold, &f1, &f2, &f3, &f4);

/*
   ofstream swq6;
   swq6.open("local_lstp_predicted_multiple.dat",ios::binary);
   if(!swq6)
      {
         cout<<"cannot open local_lstp_predicted_multiple.dat"<<endl;
         abort();
      }

   for(ix=0;ix<np;ix++)
     for(it=0;it<lt;it++)
         swq6.write((char*)&umtp[ix][it],sizeof(umtp[ix][it]));   

   swq6.close();
*/

//  cout<<"22222222"<<endl;

//calculate the semblance of the original data and the predicted multiple

//int cal_semblance(float **u,float **b,float **c,float *p,float *x,int *interc,float **sum1,float **sem,int lt,int np,int N,float theata_max,float theata_min,float dt,float dx,float v,float threhold,int win_time)

   cal_semblance(ud, tp_tmp1, Sd, Sd1, pr, x, interc, sum1, sem, lt,np, nx_l1, theta_max,theta_min,dt,dx,v,threshold,win_time);

   cal_semblance(um, tp_tmp2, Sm, Sm1, pr, x, interc, sum1, sem, lt,np, nx_l1, theta_max,theta_min,dt,dx,v,threshold,win_time);
/*
   ofstream swq7;
   swq7.open("Semblance_original_shot.dat",ios::binary);
   if(!swq7)
      {
         cout<<"cannot open Semblance_original_shot.dat"<<endl;
         abort();
      }
   for(ix=0;ix<np;ix++)
     for(it=0;it<lt;it++)
         swq7.write((char*)&Sd[ix][it],sizeof(Sd[ix][it]));
   swq7.close();

   ofstream swq8;
   swq8.open("Semblance_predicted_multiple.dat",ios::binary);
   if(!swq8)
      {
         cout<<"cannot open Semblance_predicted_multiple.dat"<<endl;
         abort();
      }
   for(ix=0;ix<np;ix++)
     for(it=0;it<lt;it++)
         swq8.write((char*)&Sm[ix][it],sizeof(Sm[ix][it]));
   swq8.close();
*/

   Butterworth_filter(Sp, Sd, Sm,  np,  lt,  n, epi, lamda);   
/*
    ofstream swq9;
   swq9.open("Butterworth_Semblance_predicted_multiple_n4.dat",ios::binary);
   if(!swq9)
      {
         cout<<"cannot open Butterworth_Semblance_predicted_multiple.dat"<<endl;
         abort();
      }
   for(ix=0;ix<np;ix++)
     for(it=0;it<lt;it++)
         swq9.write((char*)&Sp[ix][it],sizeof(Sp[ix][it]));
   swq9.close();
*/

//  cout<<"22222222"<<endl;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
          uptp[ip][it]=udtp[ip][it]*Sp[ip][it];
          umtp1[ip][it]=udtp[ip][it]*Sm[ip][it]; 
       }

//inverse tau-p transform to get the initialization of primary and the matched multiples
   inverse_tau_p(uptp, p, lt, np);
   inverse_tau_p(umtp1, m, lt, np);

   cout<<"Initializatioin Done ..."<<endl; 

   cout<<"Iteration Begins ..."<<endl; 

   while(ite_no<it_max)
    {
      cout<<"Iterative Time is===="<<ite_no<<endl;

      for(it=0;it<lt;it++)
        r[it]=u[it]-p[it]-m[it];
      l2norm_err(l2err, r,lt);
      
      if(l2err>err)
        {
//fix the matched multiple and update the primary
          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
              r1[ip][it]=ud[ip][it]-uptp[ip][it]-umtp1[ip][it];            
                
          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
              uptp[ip][it]+=r1[ip][it];

          for(ip=0;ip<np;ip++)
           for(it=0;it<lt;it++)
              {
               if(Sp[ip][it]==0)
                   uptp[ip][it]=0.0;
               else
                {
                   if(Sp[ip][it]<Sm[ip][it])
                      uptp[ip][it]=0.0;
                   else
                      uptp[ip][it]*=(1-Sm[ip][it]/Sp[ip][it]);  
                }
              }
         inverse_tau_p(uptp, p, lt, np);

         for(it=0;it<lt;it++)
           p_final[it]=p[it];

//fix the primary and update the matched multiple
         for(ip=0;ip<np;ip++)
           for(it=0;it<lt;it++)
             r2[ip][it]=ud[ip][it]-uptp[ip][it]-umtp1[ip][it];

         for(ip=0;ip<np;ip++)
           for(it=0;it<lt;it++)
             umtp1[ip][it]+=r2[ip][it];

         for(ip=0;ip<np;ip++)
           for(it=0;it<lt;it++)
             {
              if(Sm[ip][it]==0)
                   umtp1[ip][it]=0.0;
               else
                {
                   if(Sm[ip][it]<Sp[ip][it])
                      umtp1[ip][it]=0.0;
                   else
                      umtp1[ip][it]*=(1-Sp[ip][it]/Sm[ip][it]);
                }
             }

        inverse_tau_p(umtp1, m, lt, np);

        for(it=0;it<lt;it++)
           m_final[it]=m[it];

        ite_no=ite_no+1;
      }
    else
     {
//       cout<<"Maximum Iterative Time Reaches ..."<<endl;

       for(it=0;it<lt;it++)
         {
            p_final[it]=p[it];
            m_final[it]=m[it];
         }
       ite_no=it_max;
     }
   }

   return 0;
}

//================================================================================================================//
//                                               main function                                                    //
//================================================================================================================//

int main()
{
   char fn1[256], fn2[256], fn3[256],fn4[256];
   int ns, nt, lt;
   float dx, dt, v;
   int nx_h, np_h;
   float err;
   int it_max;
   float theta_max, theta_min, threshold, f1, f2, f3, f4;
   int n;
   float epi,lamda;
   int win_time;

   ifstream swq;
   swq.open("MCA_Based_PM_Separation.par");
   if(!swq)
     {
        cout<<"Cannot open MCA_Based_PM_Separation.par"<<endl;
        abort();
     }
   
   swq>>fn1>>fn2>>fn3>>fn4>>ns>>nt>>lt>>dx>>dt>>v>>nx_h>>np_h>>err>>it_max>>theta_max>>theta_min>>threshold>>f1>>f2>>f3>>f4>>n>>epi>>lamda>>win_time;
   swq.close();

   int np,nx_l;
   np=2*np_h+1;
   nx_l=2*nx_h+1;

   float l2err;

   int is,ir,ix,ip,it;

//int MCA_Based_PM_Separation(float *u, float *p, float *m, float *r, float *p_final, float *m_final, float **ud, float **um, float **up, float **Sd, float **Sm, float **Sp, float  **Sd1, float **Sm1, float **udtp, float **umtp, float **uptp, float **r1, float **r2, float **umtp1, complex<float> **cudtp, complex<float> **cumtp, int lt, float dx, float dt, float v, int  np_h, int nx_h, int np, int nx_l, float theta_max, float theta_min, float threshold, float f1, float f2, float f3, float f4, float *p, float *x, int *interc, float **sum1, int win_time, int it_max, float err, float l2err, int n, float epi, float lamda, float **sem)

   float *u;
   u=alloc1float(lt);
   float *p;
   p=alloc1float(lt);
   float *m;
   m=alloc1float(lt);

   float *r;
   r=alloc1float(lt);

   float *p_final;
   p_final=alloc1float(lt);
   float *m_final;
   m_final=alloc1float(lt);

   float **us;
   us=alloc2float(lt,nt);
   float **umpre;
   umpre=alloc2float(lt,nt);
 
   float **us1;
   us1=alloc2float(lt,nt+nx_l);
   float **umpre1;
   umpre1=alloc2float(lt,nt+nx_l);

   float **ud;
   ud=alloc2float(lt,nx_l);
   float **um;
   um=alloc2float(lt,nx_l); 
   float **up;
   up=alloc2float(lt,nx_l);

   float **tp_tmp1;
   tp_tmp1=alloc2float(lt,np);
   float **tp_tmp2;
   tp_tmp2=alloc2float(lt,np); 

   float **Sd;
   Sd=alloc2float(lt,np);
   float **Sm;
   Sm=alloc2float(lt,np);   
   float **Sp;
   Sp=alloc2float(lt,np);

   float **Sd1;
   Sd1=alloc2float(lt,np);
   float **Sm1;
   Sm1=alloc2float(lt,np);

   float **sem;
   sem=alloc2float(lt,np); 

   float **udtp;
   udtp=alloc2float(lt,np);
   float **umtp;
   umtp=alloc2float(lt,np);  
   float **uptp;
   uptp=alloc2float(lt,np);
 
   float **r1;
   r1=alloc2float(lt,np);
   float **r2;
   r2=alloc2float(lt,np);

   float **umtp1;
   umtp1=alloc2float(lt,np);

   complex<float> **cudtp;
   cudtp=alloc2complex(lt,np);
   complex<float> **cumtp;
   cumtp=alloc2complex(lt,np);

   float **sum1;
   sum1=alloc2float(lt,np);

   float pmax,pmin,dp;
   pmin=sin(theta_min*2*pai/360.0)/v;
   pmax=sin(theta_max*2*pai/360.0)/v;
   dp=(pmax-pmin)/(np-1);

   float *pr=new float [np];
   for(ip=0;ip<np;ip++)
        pr[ip]=pmin+ip*dp;

   float *x=new float [nx_l];
   for(ix=0;ix<nx_l;ix++)
        x[ix]=(ix-(nx_l-1)/2)*dx;

   int *interc=new int [nx_l];

   float mem;
   mem=0.0;
   mem+=(lt*4+nt*lt*2+(nt+nx_l)*lt*2+lt*nx_l*3+lt*np*16+np+nt*2);
   mem=mem*4/1024.0/1024.0;
   cout<<"Memory needed to be allocated is====  "<<mem<<"MB"<<endl;
 
   int tmp;

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

   for(is=0;is<ns;is++)
      {
         for(ir=0;ir<nt;ir++)
           for(it=0;it<lt;it++)
             {
                swq1.read((char*)&us[ir][it],sizeof(us[ir][it]));
                swq2.read((char*)&umpre[ir][it],sizeof(umpre[ir][it]));
             }         

         for(ir=0;ir<nt+2*nx_h+1;ir++)
          for(it=0;it<lt;it++)
            us1[ir][it]=0.0;
        for(ir=nx_h;ir<nt+nx_h;ir++)
          for(it=0;it<lt;it++)
            us1[ir][it]=us[ir-nx_h][it];
 
         for(ir=0;ir<nt+2*nx_h+1;ir++)
          for(it=0;it<lt;it++)
            umpre1[ir][it]=0.0;
        for(ir=nx_h;ir<nt+nx_h;ir++)
          for(it=0;it<lt;it++)
            umpre1[ir][it]=umpre[ir-nx_h][it];

        for(ir=0;ir<nt;ir++) 
          {
              tmp=ir+nx_h;
              for(ix=0;ix<nx_l;ix++)
               for(it=0;it<lt;it++)
                   ud[ix][it]=us1[tmp-nx_h+ix][it];
              
              for(ix=0;ix<nx_l;ix++)
               for(it=0;it<lt;it++)
                   um[ix][it]=umpre1[tmp-nx_h+ix][it];


//int MCA_Based_PM_Separation(float *u, float *p, float *m, float *r, float *p_final, float *m_final, float **ud, float **um, float **up, float **Sd, float **Sm, float **Sp, float **Sd1, float **Sm1, float **udtp, float **umtp, float **uptp, float **r1, float **r2, float **umtp1, complex<float> **cudtp, complex<float> **cumtp, int lt, float dx, float dt, float v, int  np_h, int nx_h, int np, int nx_l, float theta_max, float theta_min, float threshold, float f1, float f2, float f3, float f4, float *pr, float *x, int *interc, float **sum1, int win_time, int it_max, float err, float l2err, int n, float epi, float lamda, float **tp_tmp1, float **tp_tmp2, float **sem)
  
            MCA_Based_PM_Separation(u, p, m, r, p_final, m_final, ud, um, up, Sd, Sm, Sp, Sd1, Sm1, udtp, umtp, uptp, r1, r2, umtp1, cudtp, cumtp, lt,  dx,  dt, v, np_h, nx_h, np, nx_l, theta_max, theta_min, threshold, f1, f2, f3, f4, pr, x,interc, sum1, win_time, it_max, err, l2err, n, epi, lamda, tp_tmp1,tp_tmp2,sem);

              for(it=0;it<lt;it++)
                {            
                  swq3.write((char*)&p_final[it],sizeof(p_final[it]));
                  swq4.write((char*)&m_final[it],sizeof(m_final[it]));
                }

               cout<<is+1<<" shots, "<<is*nt+ir+1<<" traces done..."<<endl;

          }

      }  

     cout<<"ALL DONE!"<<endl; 
 
   return 0;

}







































