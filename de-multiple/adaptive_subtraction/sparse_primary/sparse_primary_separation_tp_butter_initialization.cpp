using namespace std;

#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "math.h"


#define pai 3.14159265


int forward_tau_p(float **u,float **b,float **c,float *p,float *x,int *interc,float **sum1,float **sem,int ntrace,int lt,int np,int N,float theata_max,float theata_min,float dt,float dx,float v,float threhold,int win_time)
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
	 
       for(i=0;i<np;i++)
         for(j=0;j<lt;j++)
          {
	    c[i][j]=0.0;
	    sum1[i][j]=0.0;
	    sum3[i][j]=0.0;
	    sum4[i][j]=0.0;
	    sem[i][j]=0.0;
	  }			
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
/*	  
   for(i=0;i<np;i++)
      for(j=0;j<lt;j++)
	sem[i][j]=pow(sum3[i][j],2)/(sum4[i][j]+0.01);
*/

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
	sem[i][j]/=((N+1)*(2*win_time+1));

     for(i=0;i<np;i++)
        for(j=0;j<lt;j++)
         {
          if(sem[i][j]<threhold)
	   c[i][j]=0.0;
          else
           c[i][j]=sem[i][j];				
         }

     for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          {
              if(c[i][j]<threhold)
                 sem[i][j]=0.0;
              else
                 sem[i][j]=1.0;          
          }

   for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          b[i][j]=sum1[i][j]*sem[i][j];
 	 

	 return 0; 
}

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
    char input1[256],output1[256],output2[256],output3[256],input11[256],output11[256],output22[256],output33[256],output4[256],output9[256],output10[256],output111[256],output12[256],output13[256];
    int nshot,ntrace,np,lt,N,n,nhalf;
    float theata_max,theata_min,dt,dx,v,threhold1,threhold,epsilon;
    float  mem;
    int win_time;


     ifstream swq;
     swq.open("sparse_primary_separation_tp_butter_initialization.par");
      if(!swq)
           cout<<"cannot open para card! "<<endl;
    
       swq>>input1>>input11>>output1>>output2>>output3>>output11>>output22>>output33>>output4>>output9>>output10>>output111>>output12>>output13>>nshot>>ntrace>>lt>>np>>N>>theata_max>>theata_min>>dt>>dx>>v>>threhold1>>threhold>>win_time>>epsilon>>n>>nhalf;
     swq.close();

    float pmax,pmin,dp;
    pmin=sin(theata_min*2*pai/360.0)/v;
    pmax=sin(theata_max*2*pai/360.0)/v;
    dp=(pmax-pmin)/(np-1);
   float *p=new float [np];
   for(ip=0;ip<np;ip++)
        p[ip]=pmin+ip*dp;
   float *x=new float [N+1];
   for(ix=0;ix<N+1;ix++)
        x[ix]=(ix-N/2)*dx;
   int *interc=new int [N+1];

		 
    float **u_ini=new float *[ntrace];//read initial data
    for(ir=0;ir<ntrace;ir++)
         u_ini[ir]=new float [lt];
    float **u_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         u_ini_for[ip]=new float [lt];
    float **sem_ini=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_ini[ip]=new float [lt];
    float  *ufinal=new float [lt];
   float **u=new float *[ntrace+N];
   for(ir=0;ir<ntrace+N;ir++)
         u[ir]=new float [lt];
   float **ini_local=new float *[N+1];
   for(ir=0;ir<N+1;ir++)
         ini_local[ir]=new float [lt];
   
    float **m_ini=new float *[ntrace];//read initial data
    for(ir=0;ir<ntrace;ir++)
         m_ini[ir]=new float [lt];
    float **m_ini_for=new float *[np];
    for(ip=0;ip<np;ip++)
         m_ini_for[ip]=new float [lt];
    float **sem_m_ini=new float *[np];
    for(ip=0;ip<np;ip++)
         sem_m_ini[ip]=new float [lt];
    float  *mfinal=new float [lt];
   float **m=new float *[ntrace+N];
   for(ir=0;ir<ntrace+N;ir++)
         m[ir]=new float [lt];
   float **m_local=new float *[N+1];
   for(ir=0;ir<N+1;ir++)
         m_local[ir]=new float [lt];
   
   float **sum1=new float *[np];
   for(ip=0;ip<np;ip++)
       sum1[ip]=new float [lt];
   float **sum2=new float *[np];
   for(ip=0;ip<np;ip++)
       sum2[ip]=new float [lt];
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

          for(ir=0;ir<ntrace;ir++)
             {   
                 for(ir1=0;ir1<N+1;ir1++)
                    {
                        for(it=0;it<lt;it++)
                            ini_local[ir1][it]=u[ir+ir1][it];
                    }              
 
                 forward_tau_p(ini_local,u_ini_for,sem_ini,p,x,interc,sum1,sem,ntrace,lt,np,N, theata_max, theata_min,dt,dx,v,threhold1,win_time);  

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq7.write((char*)&sem_ini[ip][it],sizeof(sem_ini[ip][it]));

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq5.write((char*)&u_ini_for[ip][it],sizeof(u_ini_for[ip][it]));

	        for(it=0;it<lt;it++)
		    ufinal[it]=0.0; 

               inverse_tau_p(u_ini_for,ufinal,lt,np);

               for(it=0;it<lt;it++)
                  ufinal[it]/=np;								 

              for(it=0;it<lt;it++)
                   swq4.write((char*)&ufinal[it],sizeof(ufinal[it]));

                 for(ir1=0;ir1<N+1;ir1++)
                    {
                        for(it=0;it<lt;it++)
                            m_local[ir1][it]=m[ir+ir1][it];
                    }               

                 forward_tau_p(m_local,m_ini_for,sem_m_ini,p,x,interc,sum1,sem,ntrace,lt,np,N, theata_max, theata_min,dt,dx,v,threhold1,win_time);  

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq77.write((char*)&sem_m_ini[ip][it],sizeof(sem_m_ini[ip][it]));

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      swq55.write((char*)&m_ini_for[ip][it],sizeof(m_ini_for[ip][it]));

	        for(it=0;it<lt;it++)
		    mfinal[it]=0.0; 

               inverse_tau_p(m_ini_for,mfinal,lt,np);

               for(it=0;it<lt;it++)
                  mfinal[it]/=np;								 

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
                  swq13.write((char*)&semb1[ip][it],sizeof(semb1[ip][it]));

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  p0tp[ip][it]=u_ini_for[ip][it]*semb[ip][it];

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                  m0tp[ip][it]=u_ini_for[ip][it]*(1-semb[ip][it]);

              inverse_tau_p(p0tp,p0t,lt,np);
              inverse_tau_p(m0tp,m0t,lt,np);

              for(ip=0;ip<np;ip++)
                for(it=0;it<lt;it++)
                 {
                   swq9.write((char*)&p0tp[ip][it],sizeof(p0tp[ip][it]));
                   swq10.write((char*)&m0tp[ip][it],sizeof(m0tp[ip][it]));
                 }

              for(it=0;it<lt;it++)
                 {
                   swq111.write((char*)&p0t[it],sizeof(p0t[it]));
                   swq12.write((char*)&m0t[it],sizeof(m0t[it]));
                 }

            }
            
               cout<<ishot<<"   th shot tau-p finished!"<<endl;
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




















