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

int forward_xf_local_tau_p(float **dxt, complex<float> **dxf, complex<float> **dpf, float **dpt,  complex <float> *base, float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax)
{
   int ix, it, ip;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       dpf[ip][it]=(0.0,0.0);

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

//   fftwf_destroy(p1);
//   fftwf_destroy(p2);

         return 0;

         }
         
int inverse_xf_local_tau_p(float **dpt, complex<float> **dxf, complex<float> **dpf, float **dxt, complex <float> *base_T, float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax)
{
   int ix, it, ip;

   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++)
       dxf[ix][it]=(0.0,0.0);

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
        for(ip=0;ip<np;ip++)
           {
             for(it=0;it<lt;it++)
               {
                  in1[it][0]=dpt[ip][it];
                  in1[it][1]=0.0;
               }

              fftwf_execute(p1);

              for(it=0;it<lt;it++)
                {
                   dpf[ip][it].real()=out1[it][0]/(float)lt;
                   dpf[ip][it].imag()=out1[it][1]/(float)lt;
                }
           }
     }

   for(it=ifmin;it<ifmax;it++)
     {
          for(ix=0;ix<nx;ix++)
           {
             for(ip=0;ip<np;ip++)
               {
                 base_T[ip].real()=0.0;
                 base_T[ip].imag()=-p[ip]*omega[it]*x[ix];
                 base_T[ip]=exp(base_T[ip]);
               }
             for(ip=0;ip<np;ip++)
               dxf[ix][it]+=base_T[ip]*dpf[ip][it];
           }
     }

   for(ix=0;ix<nx;ix++)
     for(it=lt/2+1;it<lt;it++)
       {
          dxf[ix][it].real()=dxf[ix][lt-it].real();
          dxf[ix][it].imag()=-dxf[ix][lt-it].imag();
       }

  for(ix=0;ix<nx;ix++)
    for(it=0;it<lt;it++) 
       {
         dxf[ip][it].real()/=np;
         dxf[ip][it].imag()/=np;
       }

   if((in2==NULL)||(out2==NULL))
      cout<<"memory insufficient"<<endl;
   else
     {
        for(ix=0;ix<nx;ix++)
           {
             for(it=0;it<lt;it++)
               {
                  in2[it][0]=dxf[ix][it].real();
                  in2[it][1]=dxf[ix][it].imag();
               }

              fftwf_execute(p2);

              for(it=0;it<lt;it++)
                dxt[ix][it]=out2[it][0];
           }
     }

   //   fftwf_destroy(p1);
   //   fftwf_destroy(p2);

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
	  sem[i][j]=pow(sum3[i][j],2)/sum4[i][j];
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

int sign(float v, int a)
{
   if(v>0.0)
     a=1;
   else if(v==0.0)
     a=0;
   else
     a=-1;

   return 0;

}

int Soft_Threhold_Interation( float **pt, float **mt, float **ptp, float **mtp, float **wp, float **wm, float **ptp_tmp, float **mtp_tmp, float **p1, float **m1,  float lamda1, float lamda2, float eta, int ite_max, float err, float **p0tp, float **m0tp, float **p0t, float **m0t, float **dtp, float **dorig, int nx, int np, int lt, float **snapshot1, float **snapshot2, float **snapshot3, float **snapshot4, complex<float> **atptf, complex<float> ** ptpf, float ** atpt, complex<float> **atmtf, complex<float> ** mtpf, float ** atmt, complex<float> *base, complex<float> *base_T, complex<float> **aatptpf, complex<float> **atptf1, float ** aatptp, complex<float> **aatmtpf, complex<float> **atmtf1, float **aatmtp, complex<float> **ptpf1,complex<float> **ptf1,complex<float> **mtpf1,complex<float> **mtf1, float *omega, float *p, float *x, int *interc, int ifmin, int ifmax, float **sem_d, float **sem_m)
{
   int ip,it, ite, ix;
   float amp1,amp2;
   int a;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
          ptp[ip][it]=0.0;
          mtp[ip][it]=0.0;
          ptp_tmp[ip][it]=0.0;
          mtp_tmp[ip][it]=0.0;
       } 

  for(ix=0;ix<nx;ix++) 
    for(it=0;it<lt;it++)
     {
       pt[ix][it]=p0t[ix][it];
       mt[ix][it]=m0t[ix][it];
     }
   
   amp1=0.0;
  for(ix=0;ix<nx;ix++) 
   for(it=0;it<lt;it++)
     amp1+=pow(dorig[ix][it],2);

   cout<<"amp1 is====  "<<amp1<<endl;

   if(amp1==0)
    {//1111
     for(ip=0;ip<np;ip++)
      for(it=0;it<lt;it++)
       {
        ptp[ip][it]=0.0;
        mtp[ip][it]=0.0;
       }
    for(ix=0;ix<nx;ix++) 
      for(it=0;it<lt;it++)
      {
        pt[ix][it]=0.0;
        mt[ix][it]=0.0;
      }

      cout<<" No Iteration is Needed ..."<<endl;

    }//1111

  else
   {//2222
  
    cout<<"Soft Threshold Iteration Begins ..."<<endl;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
         wp[ip][it]=fabs(m0tp[ip][it])*lamda1/2/eta;
         wm[ip][it]=fabs(p0tp[ip][it])*lamda2/2/(1+eta);
       }

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
          ptp[ip][it]=p0tp[ip][it];
          mtp[ip][it]=m0tp[ip][it];
       }
   for(ix=0;ix<nx;ix++) 
     for(it=0;it<lt;it++)
       {
         pt[ix][it]=p0t[ix][it];
         mt[ix][it]=m0t[ix][it];
       }   

    for(ite=0;ite<ite_max;ite++) 
      {//3333
        amp2=0.0;
        err=0.0; 
        for(ix=0;ix<nx;ix++) 
         for(it=0;it<lt;it++)
          amp2+=pow(pt[ix][it]+mt[ix][it]-dorig[ix][it],2);
        err=amp2/amp1; 

        cout<<"amp2 is====  "<<amp2<<endl;

        cout<<ite+1<<" iteration time, error is  "<<err<<endl;

        if(err<0.001)
         { 
           cout<<"Iteration Done! Iteration Time is"<<ite+1<<endl;
           ite=ite_max;
         }
        else
          {//4444
            for(ip=ite*np;ip<(ite+1)*np;ip++)
              for(it=0;it<lt;it++)
                {
                   snapshot1[ip][it]=ptp[ip-ite*np][it];
                   snapshot2[ip][it]=mtp[ip-ite*np][it];
                }
//calculate ATpt and ATmt, inverse pt and mt
            
             inverse_xf_local_tau_p(ptp, atptf, ptpf, atpt , base_T, omega, p, x,  nx,  np,  lt,  ifmin, ifmax);
             inverse_xf_local_tau_p(mtp, atmtf, mtpf, atmt , base_T, omega, p, x,  nx,  np,  lt,  ifmin, ifmax);

//calculate A(ATpt) and A(ATmt), forward of the inverse results
   
            forward_xf_local_tau_p(atpt, atptf1, aatptpf, aatptp ,  base, omega, p, x, nx , np, lt, ifmin, ifmax);
            forward_xf_local_tau_p(atmt, atmtf1, aatmtpf, aatmtp ,  base, omega, p, x, nx , np, lt, ifmin, ifmax);

            for(ip=0;ip<np;ip++)          
              for(it=0;it<lt;it++)
                {
//                   ptp_tmp[ip][it]=p0tp[ip][it]+m0tp[ip][it]-aatmtp[ip][it]-aatptp[ip][it]+ptp[ip][it];    
                   ptp_tmp[ip][it]=dtp[ip][it]-aatmtp[ip][it]-aatptp[ip][it]+ptp[ip][it];    
                   mtp_tmp[ip][it]=m0tp[ip][it]+mtp[ip][it]-aatmtp[ip][it]+eta/(1+eta)*dtp[ip][it]-eta/(1+eta)*aatptp[ip][it];    
                }          


/*
            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
               {
                 ptp_tmp[ip][it]*=sem_d[ip][it];
                 mtp_tmp[ip][it]*=sem_d[ip][it];
               }
*/


            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                 {
                    if(ptp_tmp[ip][it]>0)
                      p1[ip][it]=1.0;
                    else if(ptp_tmp[ip][it]==0)
                      p1[ip][it]=0.0;
                    else  
                      p1[ip][it]=-1.0;
                 }

            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                 {
                    if(mtp_tmp[ip][it]>0)
                      m1[ip][it]=1.0;
                    else if(mtp_tmp[ip][it]==0)
                      m1[ip][it]=0.0;
                    else  
                      m1[ip][it]=-1.0;
                 }

            for(ip=0;ip<np;ip++)          
              for(it=0;it<lt;it++)
               {
                  ptp[ip][it]=fabs(ptp_tmp[ip][it])-fabs(wp[ip][it]); 
                  mtp[ip][it]=fabs(mtp_tmp[ip][it])-fabs(wm[ip][it]); 
               }

            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                {
                  if(ptp[ip][it]<0.0)
                    ptp[ip][it]=0.0; 
                  else
                    ptp[ip][it]=ptp[ip][it];         
                }
 
            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                {
                  if(mtp[ip][it]<0.0)
                    mtp[ip][it]=0.0; 
                  else
                    mtp[ip][it]=mtp[ip][it];         
                }
 

            for(ip=ite*np;ip<(ite+1)*np;ip++)
              for(it=0;it<lt;it++)
                {
                   snapshot3[ip][it]=ptp[ip-ite*np][it];
                   snapshot4[ip][it]=mtp[ip-ite*np][it];
//                   snapshot3[ip][it]=ptp_tmp[ip-ite*np][it];
//                   snapshot4[ip][it]=mtp_tmp[ip-ite*np][it];
//                   snapshot3[ip][it]=aatptp[ip-ite*np][it];
//                   snapshot4[ip][it]=aatmtp[ip-ite*np][it];
                }

             for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                {
                   ptp[ip][it]*=p1[ip][it];
                   mtp[ip][it]*=m1[ip][it];
                } 

             inverse_xf_local_tau_p(ptp, ptf1, ptpf1, pt , base_T, omega, p, x,  nx,  np,  lt,  ifmin, ifmax);
             inverse_xf_local_tau_p(mtp, mtf1, mtpf1, mt , base_T, omega, p, x,  nx,  np,  lt,  ifmin, ifmax);

          }//4444

      }//3333  

    }//2222

  return 0;
  
}

int main()
{
    int is,ir,it,ip,ix,ishot,ir1,k,kk;
    char input1[256],output1[256],output2[256],output3[256],input11[256],output11[256],output22[256],output33[256],output4[256],output9[256],output10[256],output111[256],output12[256],output13[256],output14[256],output15[256],output16[256],output17[256],outputwp[256],outputwm[256],outputsnap1[256],outputsnap2[256],outputsnap3[256],outputsnap4[256];
    int nshot,ntrace,np,lt,N,n,nhalf, ite_max;
    float theata_max,theata_min,dt,dx,v,threhold1,threhold,epsilon, lamda1, lamda2, eta, err,fmin, fmax;
    float  mem;
    int win_time;


     ifstream swq;
     swq.open("sparse_primary_separation.par");
      if(!swq)
           cout<<"cannot open para card! "<<endl;
    
       swq>>input1>>input11>>output1>>output2>>output3>>output11>>output22>>output33>>output4>>output9>>output10>>output111>>output12>>output13>>output14>>output15>>output16>>output17>>outputwp>>outputwm>>outputsnap1>>outputsnap2>>outputsnap3>>outputsnap4>>nshot>>ntrace>>lt>>np>>N>>theata_max>>theata_min>>fmin>>fmax>>dt>>dx>>v>>threhold1>>threhold>>win_time>>epsilon>>n>>nhalf>>lamda1>>lamda2>>eta>>ite_max;
     swq.close();

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

    float pmax,pmin,dp;
    pmin=sin(theata_min*2*pai/360.0)/v;
    pmax=sin(theata_max*2*pai/360.0)/v;
    dp=(pmax-pmin)/(np-1);
   float *p=new float [np];
   for(ip=0;ip<np;ip++)
        p[ip]=pmin+ip*dp;
   float *x=new float [N];
   for(ix=0;ix<N;ix++)
        x[ix]=(ix-N/2)*dx;
   int *interc=new int [N];

    float **u_ini=new float *[ntrace];//read initial data
    for(ir=0;ir<ntrace;ir++)
         u_ini[ir]=new float [lt];
    float **u_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         u_ini_for[ip]=new float [lt];
    complex<float> **u_ini_for_f=new complex<float> *[np];
    for(ip=0;ip<np;ip++)
         u_ini_for_f[ip]=new complex<float> [lt];

    float **sem_ini=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_ini[ip]=new float [lt];
    float **sem_d=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_d[ip]=new float [lt];

   float  **ufinal=new float *[N];
   for(ix=0;ix<N;ix++)
       ufinal[ix]=new float [lt];
    complex<float> **u_final_f=new complex<float> *[N];
    for(ip=0;ip<N;ip++)
         u_final_f[ip]=new complex<float> [lt];

   float **u=new float *[ntrace+N];
   for(ir=0;ir<ntrace+N;ir++)
         u[ir]=new float [lt];
   float **ini_local=new float *[N];
   for(ir=0;ir<N;ir++)
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

   float  **mfinal=new float *[N];
   for(ix=0;ix<N;ix++)
       mfinal[ix]=new float [lt];
    complex<float> **m_final_f=new complex<float> *[N];
    for(ip=0;ip<N;ip++)
         m_final_f[ip]=new complex<float> [lt];
  
   float **m=new float *[ntrace+N];
   for(ir=0;ir<ntrace+N;ir++)
         m[ir]=new float [lt];
   float **m_local=new float *[N];
   for(ir=0;ir<N;ir++)
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

   float **p0t=new float *[N];
   for(ix=0;ix<N;ix++)
     p0t[ix]=new float [lt];
   float **m0t=new float *[N];
   for(ix=0;ix<N;ix++)
     m0t[ix]=new float [lt];

   complex<float> **dxf=new complex<float> *[N];
   for(ix=0;ix<N;ix++)
     dxf[ix]=new complex<float> [lt];
   complex<float> **dpf=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     dpf[ix]=new complex<float> [lt];
   complex<float> **mxf=new complex<float> *[N];
   for(ix=0;ix<N;ix++)
     mxf[ix]=new complex<float> [lt];
   complex<float> **mpf=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     mpf[ix]=new complex<float> [lt];

   complex<float> **p0t_f=new complex<float> *[N];
   for(ix=0;ix<N;ix++)
     p0t_f[ix]=new complex<float> [lt];
   complex<float> **m0t_f=new complex<float> *[N];
   for(ix=0;ix<N;ix++)
     m0t_f[ix]=new complex<float> [lt];
  
   complex<float> **p0tp_f=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     p0tp_f[ix]=new complex<float> [lt];
   complex<float> **m0tp_f=new complex<float> *[np];
   for(ix=0;ix<np;ix++)
     m0tp_f[ix]=new complex<float> [lt];

   complex<float> *base= new complex<float> [N];
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

   float **pt=new float *[N];
   for(ix=0;ix<N;ix++)
     pt[ix]=new float [lt];

   float **mt=new float *[N];
   for(ix=0;ix<N;ix++)
     mt[ix]=new float [lt];
  
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

//atptf, ptpf, atpt, atmtf, mtpf, atmt,  aatptpf, atptf1, aatptp, aatmtpf, atmtf1, aatmtp
//atptf, ptpf, atpt, atmtf, mtpf, atmt,  aatptpf, atptf1, aatptp, aatmtpf, atmtf1, aatmtp

  complex<float> **atptf=new complex<float> *[N];
  for(ix=0;ix<N;ix++)
    atptf[ix]=new complex<float> [lt];
  complex<float> **ptpf=new complex<float> *[np];
  for(ix=0;ix<np;ix++)
    ptpf[ix]=new complex<float> [lt];
  float **atpt =new float *[N];
   for(ix=0;ix<N;ix++)
    atpt[ix]=new float [lt];

  complex<float> **atmtf=new complex<float> *[N];
  for(ix=0;ix<N;ix++)
    atmtf[ix]=new complex<float> [lt];
  complex<float> **mtpf=new complex<float> *[np];
  for(ix=0;ix<np;ix++)
    mtpf[ix]=new complex<float> [lt];
  float **atmt =new float *[N];
   for(ix=0;ix<N;ix++)
    atmt[ix]=new float [lt];

  complex<float> **aatptpf=new complex<float> *[np];
  for(ix=0;ix<np;ix++)
    aatptpf[ix]=new complex<float> [lt];
  complex<float> **atptf1=new complex<float> *[N];
  for(ix=0;ix<N;ix++)
    atptf1[ix]=new complex<float> [lt];
  float **aatptp =new float *[np];
   for(ix=0;ix<np;ix++)
    aatptp[ix]=new float [lt];

  complex<float> **aatmtpf=new complex<float> *[np];
  for(ix=0;ix<np;ix++)
    aatmtpf[ix]=new complex<float> [lt];
  complex<float> **atmtf1=new complex<float> *[N];
  for(ix=0;ix<N;ix++)
    atmtf1[ix]=new complex<float> [lt];
  float **aatmtp =new float *[np];
   for(ix=0;ix<np;ix++)
    aatmtp[ix]=new float [lt];

 complex<float> **ptpf1=new complex<float> *[np];
 for(ix=0;ix<np;ix++)
    ptpf1[ix]=new complex<float> [lt];
 complex<float> **ptf1=new complex<float> *[N];
 for(ix=0;ix<N;ix++)
    ptf1[ix]=new complex<float> [lt];
 
 complex<float> **mtpf1=new complex<float> *[np];
 for(ix=0;ix<np;ix++)
    mtpf1[ix]=new complex<float> [lt];
 complex<float> **mtf1=new complex<float> *[N];
 for(ix=0;ix<N;ix++)
    mtf1[ix]=new complex<float> [lt];

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

   fstream swqsnap1;
   swqsnap1.open(outputsnap1,ios::binary|ios::out);
   fstream swqsnap2;
   swqsnap2.open(outputsnap2,ios::binary|ios::out);

   fstream swqsnap3;
   swqsnap3.open(outputsnap3,ios::binary|ios::out);
   fstream swqsnap4;
   swqsnap4.open(outputsnap4,ios::binary|ios::out);

//       for(ishot=0;ishot<nshot;ishot++)
       for(ishot=0;ishot<1;ishot++)
         {
           for(ir=0;ir<ntrace;ir++)
             for(it=0;it<lt;it++)
                 swq1.read((char*)&u_ini[ir][it],sizeof(u_ini[ir][it]));
          for(ir=0;ir<N/2;ir++)
             for(it=0;it<lt;it++)
                  u[ir][it]=0.0;
          for(ir=N/2;ir<ntrace+N/2;ir++)
             for(it=0;it<lt;it++)
                   u[ir][it]=u_ini[ir-N/2][it];
          for(ir=ntrace+N/2;ir<ntrace+N;ir++)
             for(it=0;it<lt;it++)
                  u[ir][it]=0.0;
							  
           for(ir=0;ir<ntrace;ir++)
             for(it=0;it<lt;it++)
                 swq11.read((char*)&m_ini[ir][it],sizeof(m_ini[ir][it]));
          for(ir=0;ir<N/2;ir++)
             for(it=0;it<lt;it++)
                  m[ir][it]=0.0;
          for(ir=N/2;ir<ntrace+N/2;ir++)
             for(it=0;it<lt;it++)
                   m[ir][it]=m_ini[ir-N/2][it];
          for(ir=ntrace+N/2;ir<ntrace+N;ir++)
             for(it=0;it<lt;it++)
                  m[ir][it]=0.0;

//          for(ir=0;ir<ntrace;ir++)
          for(ir=0;ir<1;ir++)
             {//5555   
                 for(ir1=0;ir1<N;ir1++)
                    {
                        for(it=0;it<lt;it++)
                            ini_local[ir1][it]=u[ir+ir1][it];
                    }              
 
                 forward_xf_local_tau_p(ini_local, dxf, dpf, u_ini_for ,  base, omega, p, x, N , np, lt, ifmin, ifmax);

                 calc_sem(ini_local, p, x, interc, sem_ini, sum3, sum4, lt,np, N, theata_max, theata_min, dt, dx, v, threhold, win_time);

                  for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                     {
                        if(sem_ini[ip][it]<threhold1)
                          sem_d[ip][it]=0.0;
                        else
                          sem_d[ip][it]=1.0;
                     }

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      u_ini_for[ip][it]*=sem_d[ip][it];

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq7.write((char*)&sem_ini[ip][it],sizeof(sem_ini[ip][it]));

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq5.write((char*)&u_ini_for[ip][it],sizeof(u_ini_for[ip][it]));

                 inverse_xf_local_tau_p(u_ini_for, u_final_f, u_ini_for_f, ufinal , base_T, omega, p, x,  N,  np,  lt,  ifmin, ifmax);

              for(ix=0;ix<N;ix++)
                for(it=0;it<lt;it++)
                   swq4.write((char*)&ufinal[ix][it],sizeof(ufinal[ix][it]));

                 for(ir1=0;ir1<N;ir1++)
                    {
                        for(it=0;it<lt;it++)
                            m_local[ir1][it]=m[ir+ir1][it];
                    }               

                 forward_xf_local_tau_p(m_local, mxf, mpf, m_ini_for ,  base, omega, p, x, N , np, lt, ifmin, ifmax);

                 calc_sem(m_local, p, x, interc, sem_m_ini, sum3, sum4, lt,np, N, theata_max, theata_min, dt, dx, v, threhold, win_time);

                  for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                     {
                        if(sem_m_ini[ip][it]<threhold1)
                          sem_d[ip][it]=0.0;
                        else
                          sem_d[ip][it]=1.0;
                     }

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      m_ini_for[ip][it]*=sem_m[ip][it];

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq77.write((char*)&sem_m_ini[ip][it],sizeof(sem_m_ini[ip][it]));

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq55.write((char*)&m_ini_for[ip][it],sizeof(m_ini_for[ip][it]));

                 inverse_xf_local_tau_p(m_ini_for, m_final_f, m_ini_for_f,  mfinal , base_T, omega, p, x,  N,  np,  lt,  ifmin, ifmax);

              for(ix=0;ix<N;ix++)
                for(it=0;it<lt;it++)
                   swq44.write((char*)&mfinal[ix][it],sizeof(mfinal[ix][it]));
              
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
                  swq13.write((char*)&semb1[ip][it],sizeof(semb1[ip][it]));

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  p0tp[ip][it]=u_ini_for[ip][it]*semb[ip][it];

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  m0tp[ip][it]=u_ini_for[ip][it]*(1-semb[ip][it]);

              inverse_xf_local_tau_p(p0tp, p0t_f, p0tp_f, p0t , base_T, omega, p, x,  N,  np,  lt,  ifmin, ifmax);
              inverse_xf_local_tau_p(m0tp, m0t_f, m0tp_f, m0t , base_T, omega, p, x,  N,  np,  lt,  ifmin, ifmax);

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                 {
                   swq9.write((char*)&p0tp[ip][it],sizeof(p0tp[ip][it]));
                   swq10.write((char*)&m0tp[ip][it],sizeof(m0tp[ip][it]));
                 }

              for(ix=0;ix<N;ix++)
              for(it=0;it<lt;it++)
                 {
                   swq111.write((char*)&p0t[ix][it],sizeof(p0t[ix][it]));
                   swq12.write((char*)&m0t[ix][it],sizeof(m0t[ix][it]));
                 }

//         Soft_Threhold_Interation( pt, mt, ptp, mtp, wp, wm, ptp_tmp,mtp_tmp, p1, m1,  lamda1, lamda2,  eta, ite_max, err, p0tp, m0tp, p0t, m0t, dtp, dorig, nx,np,  lt, snapshot1,snapshot2, atptf, ptpf, atpt,atmtf, mtpf, atmt, base, base_T,aatptpf, atptf1, aatptp, aatmtpf, atmtf1, aatmtp, ptpf1,ptf1,mtpf1,mtf1, omega, p, x, interc, ifmin, ifmax)
  
           Soft_Threhold_Interation(pt, mt, ptp, mtp, wp, wm, ptp_tmp, mtp_tmp, p1, m1,  lamda1, lamda2, eta, ite_max, err, p0tp, m0tp, p0t, m0t, u_ini_for, ini_local, N, np, lt, snapshot1, snapshot2, snapshot3, snapshot4, atptf, ptpf, atpt, atmtf, mtpf, atmt, base, base_T, aatptpf, atptf1, aatptp, aatmtpf, atmtf1, aatmtp, ptpf1, ptf1, mtpf1, mtf1, omega, p, x, interc, ifmin, ifmax, sem_d, sem_m);

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

              for(ip=0;ip<ite_max*np;ip++)
                for(it=0;it<lt;it++)
                 {
                   swqsnap1.write((char*)&snapshot1[ip][it],sizeof(snapshot1[ip][it]));
                   swqsnap2.write((char*)&snapshot2[ip][it],sizeof(snapshot2[ip][it]));
                   swqsnap3.write((char*)&snapshot3[ip][it],sizeof(snapshot3[ip][it]));
                   swqsnap4.write((char*)&snapshot4[ip][it],sizeof(snapshot4[ip][it]));
                 }

             for(ix=0;ix<N;ix++)
               for(it=0;it<lt;it++)
                 {
                   swq16.write((char*)&pt[ix][it],sizeof(pt[ix][it]));
                   swq17.write((char*)&mt[ix][it],sizeof(mt[ix][it]));
                 }

            }//5555
            
               cout<<ishot+1<<"   th shot separation  finished!"<<endl;
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









