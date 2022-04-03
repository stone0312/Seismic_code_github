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

int forward_hr_local_tau_p(float **hrtp, float **dxt, complex<float> **dxf, complex<float> **dpf, float **dpt,  complex <float> *base, float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax, float **sem, float threhold, int ite_max, int pw_half, int tw_half, int thalf, float **sem_tmp, float *window, int hrmin, int hrmax)
{
   int ix, it, ip;
   int ite;
   float sem_max;
   int p_idx, t_idx;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       dpf[ip][it]=(0.0,0.0);

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        sem_tmp[ip][it]=sem[ip][it];

   fftwf_complex *in1,*out1,*in2,*out2;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
   in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
   out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

   fftwf_plan p1, p2;

   p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
   p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

   if((in1==NULL)||(out1==NULL))
      cout<<"memory insufficient"<<endl;
   else
     {
        for(ix=0;ix<nx;ix++)
           {
             for(it=0;it<lt;it++)
               {
                  in1[it][0]=dxt[ix][it];
                  in1[it][1]=0.0;
               }

              fftwf_execute(p1);

              for(it=0;it<lt;it++)
                {
                   dxf[ix][it].real()=out1[it][0]/(float)lt;
                   dxf[ix][it].imag()=out1[it][1]/(float)lt;
                }
           }
     }

   for(it=ifmin;it<ifmax;it++)
     {
        for(ip=0;ip<np;ip++)
           {
              for(ix=0;ix<nx;ix++)
               {
                 base[ix].real()=0.0;
                 base[ix].imag()=p[ip]*omega[it]*x[ix];
                 base[ix]=exp(base[ix]);
               }
              for(ix=0;ix<nx;ix++)
                 dpf[ip][it]+=base[ix]*dxf[ix][it];
           }
     }

   for(ip=0;ip<np;ip++)
     for(it=lt/2+1;it<lt;it++)
       {
          dpf[ip][it].real()=dpf[ip][lt-it].real();
          dpf[ip][it].imag()=-dpf[ip][lt-it].imag();
       }

  for(ip=0;ip<np;ip++)
    for(it=0;it<lt;it++) 
       {
         dpf[ip][it].real()/=nx;
         dpf[ip][it].imag()/=nx;
       }

   if((in2==NULL)||(out2==NULL))
      cout<<"memory insufficient"<<endl;
   else
     {
        for(ip=0;ip<np;ip++)
           {
             for(it=0;it<lt;it++)
               {
                  in2[it][0]=dpf[ip][it].real();
                  in2[it][1]=dpf[ip][it].imag();
               }

              fftwf_execute(p2);

              for(it=0;it<lt;it++)
                 dpt[ip][it]=out2[it][0];
           }
     }

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        hrtp[ip][it]=0.0;

   for(ite=0;ite<ite_max;ite++)
     {
        sem_max=0.0;

        for(ip=0;ip<np;ip++)
          for(it=0;it<lt;it++)
            {
              if(sem_tmp[ip][it]>sem_max)
                {
                   sem_max=sem_tmp[ip][it];
                   p_idx=ip;
                   t_idx=it;
                }
            }

         cout<<"Iteration Time is ==== "<<ite+1<<" , "<<"Maximum Semblance is ==== "<<sem_max<<endl;
         cout<<"Tau and P Index are ==== "<<"( "<<t_idx<<" , "<<p_idx<<" ) "<<endl;

         if(sem_max<threhold)
           ite=ite_max;

         else
           {//xxxx
            if(t_idx-tw_half<0)       
             {  
/*
               if(p_idx-hrmin<0)
                 {
                  for(ip=0;ip<p_idx+hrmax+1;ip++)
                    for(it=0;it<t_idx+thalf+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                 }

               else if(p_idx+hrmax+1>np)
                 {
                  for(ip=p_idx-hrmin;ip<np;ip++)
                    for(it=0;it<t_idx+thalf+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                 }
 
              else
               {
               for(ip=p_idx-hrmin;ip<p_idx+hrmax+1;ip++)
                for(it=0;it<t_idx+thalf+1;it++)
                  hrtp[ip][it]=dpt[ip][it]*window[it];
               }
*/
            if(p_idx-pw_half<0)
              {
                  for(ip=0;ip<p_idx+hrmax+1;ip++)
                    for(it=0;it<t_idx+thalf+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];

                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=0;it<t_idx+tw_half+1;it++)
                    sem_tmp[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                for(ip=p_idx-hrmin;ip<np;ip++)
                  for(it=0;it<t_idx+thalf+1;it++)
                    hrtp[ip][it]=dpt[ip][it]*window[it];

                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=0;it<t_idx+tw_half+1;it++)
                    sem_tmp[ip][it]=0.0;
              }
            else
             {
               for(ip=p_idx-hrmin;ip<p_idx+hrmax+1;ip++)
                for(it=0;it<t_idx+thalf+1;it++)
                  hrtp[ip][it]=dpt[ip][it]*window[it];

              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=0;it<t_idx+tw_half+1;it++)
                  sem_tmp[ip][it]=0.0;
             }
            }
 
           else if(t_idx+tw_half+1>lt) 
            {
/*  
               if(p_idx-1<0)
                 {
                  for(ip=0;ip<p_idx+hrmax+1;ip++)
                    for(it=t_idx-thalf;it<lt;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                 }

               else if(p_idx+hrmax+1>np)
                 {
                  for(ip=p_idx-hrmin;ip<np;ip++)
                    for(it=t_idx-thalf;it<lt;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                 }
 
              else
               {
               for(ip=p_idx-hrmin;ip<p_idx+hrmax+1;ip++)
                for(it=t_idx-thalf;it<lt;it++)
                  hrtp[ip][it]=dpt[ip][it]*window[it];
               }
*/

            if(p_idx-pw_half<0)
              {
                  for(ip=0;ip<p_idx+hrmax+1;ip++)
                    for(it=t_idx-thalf;it<lt;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];

                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=t_idx-tw_half;it<lt;it++)
                    sem_tmp[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                  for(ip=p_idx-hrmin;ip<np;ip++)
                    for(it=t_idx-thalf;it<lt;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=t_idx-tw_half;it<lt;it++)
                    sem_tmp[ip][it]=0.0;
              }
            else
             {
               for(ip=p_idx-hrmin;ip<p_idx+hrmax+1;ip++)
                for(it=t_idx-thalf;it<lt;it++)
                  hrtp[ip][it]=dpt[ip][it]*window[it];
              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=t_idx-tw_half;it<lt;it++)
                  sem_tmp[ip][it]=0.0;
             }
           } 
          
          else
          {
/*
               if(p_idx-hrmin<0)
                 {
                  for(ip=0;ip<p_idx+hrmax+1;ip++)
                    for(it=t_idx-thalf;it<t_idx+tw_half+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                 }

               else if(p_idx+hrmax+1>np)
                 {
                  for(ip=p_idx-hrmin;ip<np;ip++)
                    for(it=t_idx-thalf;it<t_idx+tw_half+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                 }
 
              else
               {
               for(ip=p_idx-hrmin;ip<p_idx+hrmax+1;ip++)
                for(it=t_idx-thalf;it<t_idx+tw_half+1;it++)
                  hrtp[ip][it]=dpt[ip][it]*window[it];
               }
*/

            if(p_idx-pw_half<0)
              {
                  for(ip=0;ip<p_idx+hrmax+1;ip++)
                    for(it=t_idx-thalf;it<t_idx+tw_half+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem_tmp[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                  for(ip=p_idx-hrmin;ip<np;ip++)
                    for(it=t_idx-thalf;it<t_idx+tw_half+1;it++)
                      hrtp[ip][it]=dpt[ip][it]*window[it];
                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem_tmp[ip][it]=0.0;
              }
            else
             {
               for(ip=p_idx-hrmin;ip<p_idx+hrmax+1;ip++)
                for(it=t_idx-thalf;it<t_idx+tw_half+1;it++)
                  hrtp[ip][it]=dpt[ip][it]*window[it];
              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                  sem_tmp[ip][it]=0.0;
             }
           }

         }//xxxx

     }
     
//      cout<<"Iteration Time is ==== "<<ite+1<<endl;

         return 0;

   }
        
int inverse_local_tau_p(float *uitp, float **utp, int np, int lt)
{
   int ip, it;
   for(it=0;it<lt;it++)
     uitp[it]=0.0;

   for(it=0;it<lt;it++)
     {
        for(ip=0;ip<np;ip++)
          uitp[it]+=utp[ip][it];
     }

  return 0;

}
 
      
int calc_sem(float **u, float *p,float *x,int *interc,float **sem, float **sum3, float **sum4, int lt,int np,int N,float theata_max,float theata_min,float dt,float dx,float v,float threhold,int win_time)
{
   int i,j,k,l,ll,ii,jj,kk;


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
	    for(j=win_time;j<lt-win_time;j++)
		{
		   k=j-win_time; 
                   for(l=k;l<k+win_time*2+1;l++)
		      {
			 for(ll=0;ll<N;ll++)
			    interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
												 
			 for(ii=0;ii<N;ii++)
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
		    for(ll=0;ll<N;ll++)
		      interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
		    for(ii=0;ii<N;ii++)
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
			for(ll=0;ll<N;ll++)
			  interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
			for(ii=0;ii<N;ii++)
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
      {
        if(sum4[i][j]==0.0)
          sem[i][j]=0.0;
        else
	  sem[i][j]=pow(sum3[i][j],2)/(sum4[i][j]+0.01);
      }


   for(i=0;i<np;i++)
      for(j=0;j<lt;j++)
	sem[i][j]/=(N*(2*win_time+1));

     for(i=0;i<np;i++)
        for(j=0;j<lt;j++)
         {
          if(sem[i][j]<threhold)
	   sem[i][j]=0.0;
          else
           sem[i][j]=sem[i][j];				
         }

	 return 0; 
}


int Butterworth(float **semd, float **semm, float **semb, float epsilon, int n, int np, int lt,float threhold)
{
   int ip,it;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       semb[ip][it]=0.0;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
         {
            if(semd[ip][it]>threhold)
               semb[ip][it]=1/sqrt(1+pow(semm[ip][it]/semd[ip][it]/epsilon,n));
         }

   return 0;

}

int Smoothing_Semblance(float **semb, float **semb1, int nhalf, int np, int lt)
{
   int ip,it,it1;
   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       semb1[ip][it]=0.0;

   for(ip=0;ip<np;ip++)
    {
      for(it=nhalf;it<lt-nhalf;it++)
         {
            for(it1=it-nhalf;it1<it+nhalf+1;it1++)
              semb1[ip][it]+=semb[ip][it1];
            semb1[ip][it]/=(2*nhalf+1);
         }
    }
 
  return 0;

}

int main()
{
    int is,ir,it,ip,ix,ishot,ir1,k,kk;
    char input1[256],output1[256],output2[256],output3[256],input11[256],output11[256],output22[256],output33[256],output4[256],output9[256],output10[256],output111[256],output12[256],output13[256],output14[256],output15[256],output16[256],output17[256],outputwp[256],outputwm[256],outputsnap1[256],outputsnap2[256],outputsnap3[256],outputsnap4[256],output555[256], output5555[256];
    int nshot,ntrace,np,lt,nx_half,n,nhalf, ite_max, pw_half, tw_half, thalf,taper, N, a1, a2;
    float theata_max,theata_min,dt,dx,v,threhold1,threhold,epsilon, lamda1, lamda2, eta, err,fmin, fmax;
    float  mem;
    int win_time;
    int tlen_tmp, p_idx, t_idx, ite, hrmin, hrmax;
    float sem_max;

     ifstream swq;
     swq.open("sparse_primary_separation_tp.par");
      if(!swq)
           cout<<"cannot open para card! "<<endl;
    
       swq>>input1>>input11>>output1>>output2>>output3>>output11>>output22>>output33>>output4>>output9>>output10>>output111>>output12>>output13>>output14>>output15>>output16>>output17>>outputwp>>outputwm>>outputsnap1>>outputsnap2>>outputsnap3>>outputsnap4>>output555>>output5555>output5555>>nshot>>ntrace>>lt>>np>>nx_half>>theata_max>>theata_min>>fmin>>fmax>>dt>>dx>>v>>threhold1>>threhold>>win_time>>epsilon>>n>>nhalf>>lamda1>>lamda2>>eta>>ite_max>>pw_half>>tw_half>>thalf>>taper>>a1>>a2>>hrmin>>hrmax;
     swq.close();
   
   float *win_off=new float [ntrace];
   for(ir=0;ir<ntrace;ir++)
      win_off[ir]=1.0;

   for(ir=0;ir<a1;ir++)
      win_off[ir]=sin(pai*ir/(2*a1));
   for(ir=a2;ir<ntrace;ir++)
      win_off[ir]=sin(pai*(ntrace-1-ir)/(2*(ntrace-1-a2)));

   float *tlen=new float [ntrace];
   for(ir=0;ir<ntrace;ir++)
     tlen[ir]=(int)(thalf*win_off[ir]+0.5);

/*
    for(ir=0;ir<ntrace;ir++)
      cout<<ir<<"   "<<win_off[ir]<<endl;

    for(ir=0;ir<ntrace;ir++)
      cout<<ir<<"   "<<tlen[ir]<<endl;
    return 0;
*/

   float *omega;
   omega=alloc1float(lt);
   for(it=0;it<lt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*lt);
   for(it=lt/2+1;it<lt;it++)
      omega[it]=-2*pai*1000/(2*dt)+2*pai*(it-lt/2)*1000/(dt*lt);

   int ifmin, ifmax;
   ifmin=int(fmin*dt*lt/1000);
   ifmax=int(fmax*dt*lt/1000)+1;

   cout<<ifmin<<" ,  "<<ifmax<<endl;
   cout<<"totally "<<ifmax-ifmin<<" frequency slices need to be calculated"<<endl;

   float *window=new float [2*thalf+1];
   for(it=0;it<2*thalf+1;it++)
      window[it]=1.0;
/*
   for(it=0;it<taper;it++)
      window[it]=sqrt(sin(pai*it/(2*taper)));
   for(it=2*thalf+1-taper;it<2*thalf+1;it++)
      window[it]=window[2*thalf-it];
*/

/*
    for(it=0;it<2*thalf+1;it++)
      cout<<it<<"   "<<window[it]<<endl;
*/
    float pmax,pmin,dp;
    pmin=sin(theata_min*2*pai/360.0)/v;
    pmax=sin(theata_max*2*pai/360.0)/v;
    dp=(pmax-pmin)/(np-1);
   float *p=new float [np];
   for(ip=0;ip<np;ip++)
        p[ip]=pmin+ip*dp;
   float *x=new float [2*nx_half+1];
   for(ix=0;ix<2*nx_half+1;ix++)
        x[ix]=(ix-nx_half)*dx;
   int *interc=new int [2*nx_half+1];

    float **u_ini=new float *[ntrace];//read initial data
    for(ir=0;ir<ntrace;ir++)
         u_ini[ir]=new float [lt];
    float **u_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         u_ini_for[ip]=new float [lt];
    complex<float> **u_ini_for_f=new complex<float> *[np];
    for(ip=0;ip<np;ip++)
         u_ini_for_f[ip]=new complex<float> [lt];

    float **hr_u_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         hr_u_ini_for[ip]=new float [lt];
    float **hr_m_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         hr_m_ini_for[ip]=new float [lt];
    float **sem_ini=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_ini[ip]=new float [lt];
    float **sem_d=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_d[ip]=new float [lt];

    float **sem_tmp=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_tmp[ip]=new float [lt];

    float **sem_tmp1=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_tmp1[ip]=new float [lt];

   float  *ufinal=new float [lt];
    complex<float> **u_final_f=new complex<float> *[2*nx_half+1];
    for(ip=0;ip<2*nx_half+1;ip++)
         u_final_f[ip]=new complex<float> [lt];

   float **u=new float *[ntrace+2*nx_half+1];
   for(ir=0;ir<ntrace+2*nx_half+1;ir++)
         u[ir]=new float [lt];
   float **ini_local=new float *[2*nx_half+1];
   for(ir=0;ir<2*nx_half+1;ir++)
         ini_local[ir]=new float [lt];
   
    float **m_ini=new float *[ntrace];//read initial data
    for(ir=0;ir<ntrace;ir++)
         m_ini[ir]=new float [lt];
    float **m_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         m_ini_for[ip]=new float [lt];
    complex<float> **m_ini_for_f=new complex<float> *[np];
    for(ip=0;ip<np;ip++)
         m_ini_for_f[ip]=new complex<float> [lt];
    float **sem_m_ini=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_m_ini[ip]=new float [lt];
    float **sem_m=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_m[ip]=new float [lt];

    float **sem1=new float *[np];
    for(ip=0;ip<np;ip++)
         sem1[ip]=new float [lt];
    float **sem2=new float *[np];
    for(ip=0;ip<np;ip++)
         sem2[ip]=new float [lt];
    float **sem3=new float *[np];
    for(ip=0;ip<np;ip++)
         sem3[ip]=new float [lt];

   float  *mfinal=new float [lt];
    complex<float> **m_final_f=new complex<float> *[2*nx_half+1];
    for(ip=0;ip<2*nx_half+1;ip++)
         m_final_f[ip]=new complex<float> [lt];
  
   float **m=new float *[ntrace+2*nx_half+1];
   for(ir=0;ir<ntrace+2*nx_half+1;ir++)
         m[ir]=new float [lt];
   float **m_local=new float *[2*nx_half+1];
   for(ir=0;ir<2*nx_half+1;ir++)
         m_local[ir]=new float [lt];
   
   float **sum3=new float *[np];
   for(ip=0;ip<np;ip++)
       sum3[ip]=new float [lt];
   float **sum4=new float *[np];
   for(ip=0;ip<np;ip++)
       sum4[ip]=new float [lt];
   float **sem=new float *[np];
   for(ip=0;ip<np;ip++)
       sem[ip]=new float [lt];

   float **semb=new float *[np];
   for(ip=0;ip<np;ip++)
       semb[ip]=new float [lt];

   float **semb1=new float *[np];
   for(ip=0;ip<np;ip++)
       semb1[ip]=new float [lt];
   
   float **p0tp=new float *[np];
   for(ip=0;ip<np;ip++)
      p0tp[ip]=new float [lt];
   float **m0tp=new float *[np];
   for(ip=0;ip<np;ip++)
      m0tp[ip]=new float [lt];

   float *p0t=new float [lt];
   float *m0t=new float [lt];

   complex<float> **dxf=new complex<float> *[2*nx_half+1];
   for(ix=0;ix<2*nx_half+1;ix++)
     dxf[ix]=new complex<float> [lt];
   complex<float> **dpf=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     dpf[ix]=new complex<float> [lt];
   complex<float> **mxf=new complex<float> *[2*nx_half+1];
   for(ix=0;ix<2*nx_half+1;ix++)
     mxf[ix]=new complex<float> [lt];
   complex<float> **mpf=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     mpf[ix]=new complex<float> [lt];

   complex<float> **p0t_f=new complex<float> *[2*nx_half+1];
   for(ix=0;ix<2*nx_half+1;ix++)
     p0t_f[ix]=new complex<float> [lt];
   complex<float> **m0t_f=new complex<float> *[2*nx_half+1];
   for(ix=0;ix<2*nx_half+1;ix++)
     m0t_f[ix]=new complex<float> [lt];
  
   complex<float> **p0tp_f=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     p0tp_f[ix]=new complex<float> [lt];
   complex<float> **m0tp_f=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     m0tp_f[ix]=new complex<float> [lt];

   complex<float> *base= new complex<float> [2*nx_half+1];
   complex<float> *base_T= new complex<float> [np];
   
//array for soft threshold interation

   float **wp=new float *[np];
   for(ip=0;ip<np;ip++)
      wp[ip]=new float [lt];

   float **wm=new float *[np];
   for(ip=0;ip<np;ip++)
      wm[ip]=new float [lt];

   float **ptp_tmp=new float *[np];
   for(ip=0;ip<np;ip++)
      ptp_tmp[ip]=new float [lt];

   float **mtp_tmp=new float *[np];
   for(ip=0;ip<np;ip++)
      mtp_tmp[ip]=new float [lt];

   float **ptp=new float *[np];
   for(ip=0;ip<np;ip++)
      ptp[ip]=new float [lt];

   float **mtp=new float *[np];
   for(ip=0;ip<np;ip++)
      mtp[ip]=new float [lt];

   float **p1=new float *[np];
   for(ip=0;ip<np;ip++)
      p1[ip]=new float [lt];

   float **m1=new float *[np];
   for(ip=0;ip<np;ip++)
      m1[ip]=new float [lt];

   float *pt=new float [lt];
   float *mt=new float [lt];
  
   float **snapshot1=new float *[ite_max*np];
   for(ip=0;ip<ite_max*np;ip++)
      snapshot1[ip]=new float [lt];

   float **snapshot2=new float *[ite_max*np];
   for(ip=0;ip<ite_max*np;ip++)
      snapshot2[ip]=new float [lt];

   float **snapshot3=new float *[ite_max*np];
   for(ip=0;ip<ite_max*np;ip++)
      snapshot3[ip]=new float [lt];

   float **snapshot4=new float *[ite_max*np];
   for(ip=0;ip<ite_max*np;ip++)
      snapshot4[ip]=new float [lt];


 float *ptfinal=new float [lt];
 float *mtfinal=new float [lt];

 // open input data
   ifstream swq1;
   swq1.open(input1,ios::binary);
   if(!swq1)
     cout<<"cannot open"<<input1<<endl;
   fstream swq4;
   swq4.open(output2,ios::binary|ios::out);
   fstream swq5;
   swq5.open(output1,ios::binary|ios::out);
   fstream swq7;
   swq7.open(output3,ios::binary|ios::out);
	
   ifstream swq11;
   swq11.open(input11,ios::binary);
   if(!swq11)
     cout<<"cannot open"<<input11<<endl;
   fstream swq44;
   swq44.open(output22,ios::binary|ios::out);
   fstream swq55;
   swq55.open(output11,ios::binary|ios::out);
   fstream swq77;
   swq77.open(output33,ios::binary|ios::out);

   fstream swq8;
   swq8.open(output4,ios::binary|ios::out);

   fstream swq9;
   swq9.open(output9,ios::binary|ios::out);
   fstream swq10;
   swq10.open(output10,ios::binary|ios::out);
   fstream swq111;
   swq111.open(output111,ios::binary|ios::out);
   fstream swq12;
   swq12.open(output12,ios::binary|ios::out);
   fstream swq13;
   swq13.open(output13,ios::binary|ios::out);

   fstream swq14;
   swq14.open(output14,ios::binary|ios::out);
   fstream swq15;
   swq15.open(output15,ios::binary|ios::out);
   fstream swq16;
   swq16.open(output16,ios::binary|ios::out);
   fstream swq17;
   swq17.open(output17,ios::binary|ios::out);

   fstream swqwp;
   swqwp.open(outputwp,ios::binary|ios::out);
   fstream swqwm;
   swqwm.open(outputwm,ios::binary|ios::out);
   fstream swq555;
   swq555.open(output555,ios::binary|ios::out);
   fstream swq5555;
   swq5555.open(output5555,ios::binary|ios::out);

//       for(ishot=0;ishot<nshot;ishot++)
       for(ishot=0;ishot<1;ishot++)
         {
           cout<<ishot+1<<"th shot separation begins ..."<<endl;

           for(ir=0;ir<ntrace;ir++)
             for(it=0;it<lt;it++)
                 swq1.read((char*)&u_ini[ir][it],sizeof(u_ini[ir][it]));
          for(ir=0;ir<nx_half;ir++)
             for(it=0;it<lt;it++)
                  u[ir][it]=0.0;
          for(ir=nx_half;ir<ntrace+nx_half;ir++)
             for(it=0;it<lt;it++)
                   u[ir][it]=u_ini[ir-nx_half][it];
          for(ir=ntrace+nx_half;ir<ntrace+2*nx_half+1;ir++)
             for(it=0;it<lt;it++)
                  u[ir][it]=0.0;
							  
           for(ir=0;ir<ntrace;ir++)
             for(it=0;it<lt;it++)
                 swq11.read((char*)&m_ini[ir][it],sizeof(m_ini[ir][it]));

          for(ir=0;ir<nx_half;ir++)
             for(it=0;it<lt;it++)
                  m[ir][it]=0.0;
          for(ir=nx_half;ir<ntrace+nx_half;ir++)
             for(it=0;it<lt;it++)
                   m[ir][it]=m_ini[ir-nx_half][it];
          for(ir=ntrace+nx_half;ir<ntrace+2*nx_half+1;ir++)
             for(it=0;it<lt;it++)
                  m[ir][it]=0.0;

          N=2*nx_half+1;

          for(ir=0;ir<ntrace;ir++)
//          for(ir=9;ir<10;ir++)
             {//5555

                 cout<<"    "<<ir+1<<" trace begins ..."<<endl;
   
                 for(ir1=0;ir1<2*nx_half+1;ir1++)
                    {
                        for(it=0;it<lt;it++)
                            ini_local[ir1][it]=u[ir+ir1][it];
                    }              
 
                 calc_sem(ini_local, p, x, interc, sem_ini, sum3, sum4, lt,np, N, theata_max, theata_min, dt, dx, v, threhold, win_time);

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq7.write((char*)&sem_ini[ip][it],sizeof(sem_ini[ip][it]));

                 forward_hr_local_tau_p(hr_u_ini_for, ini_local, dxf, dpf, u_ini_for ,  base, omega, p, x, N , np, lt, ifmin, ifmax,sem_ini, threhold1, ite_max, pw_half, tw_half, thalf, sem_tmp, window, hrmin, hrmax);

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq5.write((char*)&hr_u_ini_for[ip][it],sizeof(hr_u_ini_for[ip][it]));

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq555.write((char*)&u_ini_for[ip][it],sizeof(u_ini_for[ip][it]));

                 inverse_local_tau_p(ufinal, hr_u_ini_for, np,  lt );

                for(it=0;it<lt;it++)
                   swq4.write((char*)&ufinal[it],sizeof(ufinal[it]));

                 for(ir1=0;ir1<N;ir1++)
                    {
                        for(it=0;it<lt;it++)
                            m_local[ir1][it]=m[ir+ir1][it];
                    }               

                 calc_sem(m_local, p, x, interc, sem_m_ini, sum3, sum4, lt,np, N, theata_max, theata_min, dt, dx, v, threhold, win_time);

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq77.write((char*)&sem_m_ini[ip][it],sizeof(sem_m_ini[ip][it]));

                 forward_hr_local_tau_p(hr_m_ini_for, m_local, mxf, mpf, m_ini_for ,  base, omega, p, x, N , np, lt, ifmin, ifmax, sem_m_ini, threhold1, ite_max, pw_half, tw_half, thalf, sem_tmp, window, hrmin, hrmax);

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq55.write((char*)&hr_m_ini_for[ip][it],sizeof(hr_m_ini_for[ip][it]));
                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq5555.write((char*)&m_ini_for[ip][it],sizeof(m_ini_for[ip][it]));

                 inverse_local_tau_p(mfinal, hr_m_ini_for,  np,  lt);

                for(it=0;it<lt;it++)
                   swq44.write((char*)&mfinal[it],sizeof(mfinal[it]));
              
              Butterworth(sem_ini, sem_m_ini, semb,  epsilon,  n, np, lt, threhold);

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  swq8.write((char*)&semb[ip][it],sizeof(semb[ip][it]));
         
              Smoothing_Semblance(semb,semb1, nhalf, np,lt);
             
              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                   {
                      if(semb1[ip][it]<threhold)
                         semb[ip][it]=0.0;
                      else
                         semb[ip][it]=1.0;
                   }

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  {
                     sem1[ip][it]=semb[ip][it];
                     sem2[ip][it]=1-semb[ip][it];
                  }              

             for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                   sem3[ip][it]=sem_ini[ip][it]*semb[ip][it];

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  swq13.write((char*)&sem3[ip][it],sizeof(sem3[ip][it]));

/*
              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  p0tp[ip][it]=hr_u_ini_for[ip][it]*semb[ip][it];

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  m0tp[ip][it]=hr_u_ini_for[ip][it]*(1-semb[ip][it]);
*/

//iteration for initial multiple adaptive subtraction
   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
      {
         p0tp[ip][it]=0.0;
         m0tp[ip][it]=0.0;
      }

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        sem_tmp1[ip][it]=sem3[ip][it];

   for(ite=0;ite<ite_max;ite++)
     {
        sem_max=0.0;

        for(ip=0;ip<np;ip++)
          for(it=0;it<lt;it++)
            {
              if(sem_tmp1[ip][it]>sem_max)
                {
                   sem_max=sem_tmp1[ip][it];
                   p_idx=ip;
                   t_idx=it;
                }
            }

         if(sem_max<threhold1)
           ite=ite_max;
         else
           {
             cout<<"Iteration Time is ==== "<<ite+1<<" , "<<"Maximum Semblance is ==== "<<sem_max<<endl;
             cout<<"Tau and P Index are ==== "<<"( "<<t_idx<<" , "<<p_idx<<" ) "<<endl;

             tlen_tmp=tlen[ir]; 

           if(t_idx-tw_half<0)
            {
             for(it=0;it<t_idx+thalf+1;it++)
                p0tp[p_idx][it]=hr_u_ini_for[p_idx][it];

            if(p_idx-pw_half<0)
              {
                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=0;it<t_idx+tw_half+1;it++)
                    sem_tmp1[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=0;it<t_idx+tw_half+1;it++)
                    sem_tmp1[ip][it]=0.0;
              }
            else
             {
              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=0;it<t_idx+tw_half+1;it++)
                  sem_tmp1[ip][it]=0.0;
             }
           } 

          else if(t_idx+tw_half+1>lt)
           {
             for(it=t_idx-thalf;it<lt;it++)
                p0tp[p_idx][it]=hr_u_ini_for[p_idx][it];

            if(p_idx-pw_half<0)
              {
                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=t_idx-tw_half;it<lt;it++)
                    sem_tmp1[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=t_idx-tw_half;it<lt;it++)
                    sem_tmp1[ip][it]=0.0;
              }
            else
             {
              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=t_idx-tw_half;it<lt;it++)
                  sem_tmp1[ip][it]=0.0;
             }
          }
         
          else
           {
             for(it=t_idx-thalf;it<t_idx+thalf+1;it++)
                p0tp[p_idx][it]=hr_u_ini_for[p_idx][it];

            if(p_idx-pw_half<0)
              {
                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem_tmp1[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem_tmp1[ip][it]=0.0;
              }
            else
             {
              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                  sem_tmp1[ip][it]=0.0;
             }
           }

        }

     }

    for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        m0tp[ip][it]=hr_u_ini_for[ip][it]-p0tp[ip][it];

    for(ip=0;ip<np;ip++)
      for(it=0;it<lt;it++)
        {
           swq9.write((char*)&p0tp[ip][it],sizeof(p0tp[ip][it]));
           swq10.write((char*)&m0tp[ip][it],sizeof(m0tp[ip][it]));
        }

    inverse_local_tau_p(p0t, p0tp,  np,  lt);
    inverse_local_tau_p(m0t, m0tp,  np,  lt);

    for(it=0;it<lt;it++)
       {
          swq111.write((char*)&p0t[it],sizeof(p0t[it]));
          swq12.write((char*)&m0t[it],sizeof(m0t[it]));
       }
 
    cout<<"Initialization Done ..."<<endl;

    cout<<"Primary Optimization Begins ..."<<endl;
 
    for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
          ptp[ip][it]=p0tp[ip][it];
          mtp[ip][it]=m0tp[ip][it];
       }

    for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        sem_tmp1[ip][it]=sem_m_ini[ip][it];

    for(ite=0;ite<ite_max;ite++)
     {
        sem_max=0.0;

        for(ip=0;ip<np;ip++)
          for(it=0;it<lt;it++)
            {
              if(sem_tmp1[ip][it]>sem_max)
                {
                   sem_max=sem_tmp1[ip][it];
                   p_idx=ip;
                   t_idx=it;
                }
            }

         if(sem_max<threhold1)
           ite=ite_max;
         else
           {
             cout<<"Iteration Time is ==== "<<ite+1<<" , "<<"Maximum Semblance is ==== "<<sem_max<<endl;
             cout<<"Tau and P Index are ==== "<<"( "<<t_idx<<" , "<<p_idx<<" ) "<<endl;

             tlen_tmp=tlen[ir]; 

            if(p_idx-pw_half<0)
              {
                if(t_idx-tlen_tmp<0)                
                 {
                   for(ip=0;ip<p_idx+pw_half+1;ip++)
                    {
                     for(it=0;it<t_idx;it++)
                      ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                     for(it=t_idx;it<t_idx+tlen_tmp;it++)
                      ptp[ip][it]=0.0;
                    }
                 }
                else if(t_idx+tlen_tmp>lt)                
                 {
                   for(ip=0;ip<p_idx+pw_half+1;ip++)
                    {
                     for(it=t_idx-tlen_tmp;it<t_idx;it++)
                      ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                     for(it=t_idx;it<lt;it++)
                      ptp[ip][it]=0.0;
                    }
                  }
                else            
                 {
                   for(ip=0;ip<p_idx+pw_half+1;ip++)
                    {
                     for(it=t_idx-tlen_tmp;it<t_idx;it++)
                      ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                     for(it=t_idx;it<t_idx+tlen_tmp;it++)
                      ptp[ip][it]=0.0;
                     }
                 }

                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem_tmp1[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
             {
              if(t_idx-tlen_tmp<0)
               {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 {
                  for(it=0;it<t_idx;it++)
                    ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                  for(it=t_idx;it<t_idx+tlen_tmp;it++)
                   ptp[ip][it]=0.0;
                 }
               }
              else if(t_idx+tlen_tmp>lt)
               {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 {
                  for(it=t_idx-tlen_tmp;it<t_idx;it++)
                    ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                  for(it=t_idx;it<lt;it++)
                   ptp[ip][it]=0.0;
                 }
               }
              else
               {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 {
                  for(it=t_idx-tlen_tmp;it<t_idx;it++)
                    ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                  for(it=t_idx;it<t_idx+tlen_tmp;it++)
                   ptp[ip][it]=0.0;
                 }
               }

                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem_tmp1[ip][it]=0.0;
              }
            else
             {
              if(t_idx-tlen_tmp<0)
               {
                for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
                 {
                 for(it=0;it<t_idx;it++)
                   ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                 for(it=t_idx;it<t_idx+tlen_tmp;it++)
                    ptp[ip][it]=0.0;
                 }
               }
              else if(t_idx+tlen_tmp>lt)
               {
                for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
                 {
                 for(it=t_idx-tlen_tmp;it<t_idx;it++)
                   ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                 for(it=t_idx;it<lt;it++)
                    ptp[ip][it]=0.0;
                 }
               }
              else
               {
                for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
                 {
                 for(it=t_idx-tlen_tmp;it<t_idx;it++)
                   ptp[ip][it]*=float(1.0)/(float)(it-t_idx+tlen_tmp+1);
                 for(it=t_idx;it<t_idx+tlen_tmp;it++)
                    ptp[ip][it]=0.0;
                 }
               }

              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                  sem_tmp1[ip][it]=0.0;
             }

           }

     }

   cout<<"Primary Optimization Done!"<<endl;

    for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        mtp[ip][it]=hr_u_ini_for[ip][it]-ptp[ip][it];

     inverse_local_tau_p(pt, ptp, np,  lt);
     inverse_local_tau_p(mt, mtp, np,  lt);

     for(ip=0;ip<np;ip++)
       for(it=0;it<lt;it++)
         {
            swq14.write((char*)&ptp[ip][it],sizeof(ptp[ip][it]));
            swq15.write((char*)&mtp[ip][it],sizeof(mtp[ip][it]));
         }

     for(ip=0;ip<np;ip++)
       for(it=0;it<lt;it++)
         {
            swqwp.write((char*)&wp[ip][it],sizeof(wp[ip][it]));
            swqwm.write((char*)&wm[ip][it],sizeof(wm[ip][it]));
         }
 
     for(it=0;it<lt;it++)
      {
        swq16.write((char*)&pt[it],sizeof(pt[it]));
        swq17.write((char*)&mt[it],sizeof(mt[it]));
      }

     cout<<"    "<<ir+1<<" trace done ..."<<endl;

    }//5555
            
    cout<<ishot+1<<"th shot separation done!"<<endl;

   }

   swq1.close();
   swq4.close();
   swq5.close();
   swq7.close();
   swq44.close();
   swq55.close();
   swq77.close();
   swq8.close();
	 

   cout<<"all done!"<<endl;


   return 0;       
}




