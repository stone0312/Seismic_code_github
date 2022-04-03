using namespace std;

#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "mpi.h"


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
       {
           for(j=0;j<np;j++)
               ufinal[i]+=u[j][i];
       }
    return 0;  
}
  
int main(int argc,char **argv)
{
    int i,j,k,ii,jj,kk,is;
    char input1[256],inifor[256],inisem[256],iniinv[256];
    int nshot,ntrace,np,lt,N;
    float theata_max,theata_min,dt,dx,v,threhold,threhold1;
    float  mem,coff;
    int win_time;

    int npp,myid;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&npp);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);


     ifstream swq1;
     swq1.open("tp.par");
      if(!swq1)
       {
           cout<<"cannot open para card! "<<endl;
//           abort();
       }
            swq1>>input1>>inifor>>iniinv>>inisem>>nshot>>ntrace>>lt>>np>>N>>theata_max>>theata_min>>dt>>dx>>v>>threhold>>threhold1>>coff>>win_time;
     swq1.close();

		 
    float **u_ini=new float *[ntrace];//read initial data
    for(i=0;i<ntrace;i++)
         u_ini[i]=new float [lt];
    
    float **u_ini_for=new float *[np];
    for(i=0;i<np;i++)
         u_ini_for[i]=new float [lt];

    float **sem_ini=new float *[np];
    for(i=0;i<np;i++)
         sem_ini[i]=new float [lt];

    float  *ufinal=new float [lt];

     float pmax,pmin,dp;
    pmin=sin(theata_min*2*pai/360.0)/v;
    pmax=sin(theata_max*2*pai/360.0)/v;
    dp=(pmax-pmin)/(np-1);

//    cout<<"pmax,dp===="<<pmax<<","<<dp<<endl;

   float **u=new float *[ntrace+N];
   for(i=0;i<ntrace+N;i++)
         u[i]=new float [lt];
   float **ini_local=new float *[N+1];
   for(i=0;i<N+1;i++)
         ini_local[i]=new float [lt];
   float **uu=new float *[ntrace+N];
   for(i=0;i<ntrace+N;i++)
         uu[i]=new float [lt];
   float **pre_local=new float *[N+1];
   for(i=0;i<N+1;i++)
         pre_local[i]=new float [lt];
   
   float *p=new float [np];
   for(i=0;i<np;i++)
        p[i]=pmin+i*dp;

   float *x=new float [N+1];
   for(i=0;i<N+1;i++)
        x[i]=(i-N/2)*dx;

   int *interc=new int [N+1];

   float **sum1=new float *[np];
   for(i=0;i<np;i++)
       sum1[i]=new float [lt];
   float **sum2=new float *[np];
   for(i=0;i<np;i++)
       sum2[i]=new float [lt];
   float **sem=new float *[np];
   for(i=0;i<np;i++)
       sem[i]=new float [lt];

   for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          {
             sum1[i][j]=0.0;
             sum2[i][j]=0.0;
	     sem[i][j]=0.0;
          }


 // open input data
   ifstream swq2;
   swq2.open(input1,ios::binary);
   if(!swq2)
      {
        cout<<"cannot open"<<input1<<endl;
//        abort();
      }

    fstream swq4;
    swq4.open(iniinv,ios::binary|ios::out);
    if(!swq4)
       {
           cout<<"cannot open"<<iniinv<<endl;
//           abort();
       }

   fstream swq5;
   swq5.open(inifor,ios::binary|ios::out);

   fstream swq7;
   swq7.open(inisem,ios::binary|ios::out);

//	 ofstream swq888;
//	 swq888.open("ini_sem_yuedong.dat",ios::binary);

//   ofstream swq8888;
//   swq8888.open("ini_for_yuedong.dat",ios::binary);
	
       for(is=myid;is<nshot;is+=npp)
         {

           swq2.seekg(0,ios::beg);
           for(i=0;i<is;i++)
             swq2.seekg(ntrace*lt*4,ios::cur);

           for(i=0;i<ntrace;i++)
             for(j=0;j<lt;j++)
                 swq2.read((char*)&u_ini[i][j],sizeof(u_ini[i][j]));

          for(i=0;i<N/2;i++)
             for(j=0;j<lt;j++)
                  u[i][j]=0.0;
							  
          for(i=N/2;i<ntrace+N/2;i++)
             for(j=0;j<lt;j++)
                {
                   k=i-N/2;
                   u[i][j]=u_ini[k][j];
                }

          for(i=ntrace+N/2;i<ntrace+N;i++)
             for(j=0;j<lt;j++)
                  u[i][j]=0.0;
							  
          for(ii=0;ii<ntrace;ii++)
             {   
                 for(jj=0;jj<N+1;jj++)
                    {
                        for(kk=0;kk<lt;kk++)
                            ini_local[jj][kk]=u[ii+jj][kk];
                    }               
             
                 forward_tau_p(ini_local,u_ini_for,sem_ini,p,x,interc,sum1,sem,ntrace,lt,np,N, theata_max, theata_min,dt,dx,v,threhold,win_time);  

                 swq5.seekg(0,ios::beg);
                 for(i=0;i<is*ntrace+ii;i++)
                    swq5.seekg(np*lt*4,ios::cur);

                 swq7.seekg(0,ios::beg);
                 for(i=0;i<is*ntrace+ii;i++)
                    swq7.seekg(np*lt*4,ios::cur);

                 for(i=0;i<np;i++)
                   for(j=0;j<lt;j++)
                      swq7.write((char*)&sem_ini[i][j],sizeof(sem_ini[i][j]));
                 for(i=0;i<np;i++)
                   for(j=0;j<lt;j++)
                      swq5.write((char*)&u_ini_for[i][j],sizeof(u_ini_for[i][j]));

								 
		   for(i=0;i<np;i++)
                     for(j=0;j<lt;j++)
                        { 
                            sum1[i][j]=0.0;
                            sum2[i][j]=0.0;
                            sem[i][j]=0.0;
                        }

	        for(j=0;j<lt;j++)
		    ufinal[j]=0.0; 
               inverse_tau_p(u_ini_for,ufinal,lt,np);
               for(j=0;j<lt;j++)
                  ufinal[j]/=np;								 
             
              swq4.seekg(0,ios::beg);
              for(i=0;i<is*ntrace+ii;i++)
                swq4.seekg(lt*4,ios::cur);

              for(j=0;j<lt;j++)
                   swq4.write((char*)&ufinal[j],sizeof(ufinal[j]));

               for(j=0;j<lt;j++)
                   ufinal[j]=0.0;
               for(i=0;i<np;i++)
                  for(j=0;j<lt;j++)
		    sem_ini[i][j]=0.0;

               if((ii+1)%50==0)
                cout<<is+1<<"  shot, "<<ii+1<<"   th trace tau-p finished!"<<endl;
            }
            
               cout<<is<<"   th shot tau-p finished!"<<endl;
          }

    MPI_Barrier(MPI_COMM_WORLD);

   swq2.close();
   swq4.close();
   swq5.close();
   swq7.close();
	 
   for(i=0;i<ntrace;i++)
       delete[]   u_ini[i]; 
   delete[]   u_ini;  
   for(i=0;i<ntrace+N;i++)
         delete[] u[i];
   delete[] u;

   for(i=0;i<N+1;i++)
         delete[] ini_local[i];
   delete[] ini_local;


   for(i=0;i<np;i++)
         delete[] u_ini_for[i];
   delete[] u_ini_for;
   
   for(i=0;i<np;i++)
         delete[] sem_ini[i];
   delete[] sem_ini;
 

   for(i=0;i<np;i++)
         delete[] sum1[i];
   delete[] sum1;

   for(i=0;i<np;i++)
         delete[] sum2[i];
   delete[] sum2;

   for(i=0;i<np;i++)
         delete[] sem[i];
   delete[] sem;

   cout<<"all done!"<<endl;


   return 0;       
}




















