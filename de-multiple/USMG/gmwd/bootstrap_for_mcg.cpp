using namespace std;
#include <iostream>
#include "math.h"
#include <fstream>
#include "stdlib.h"
#include <complex>
#include "alloc.c"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#define pai 3.14159265
   
int bootstrap_for_mcg(float **u, int nr, int lt, int n,int B,float theta_max,float theta_min)
{
   int ix,ir,it,ib;
   int ir_beg,ir_end,ir_tmp;

   float *theta1;
   theta1=alloc1float(n);
   float *theta2;  
   theta2=alloc1float(B);
 
   float theta;
   theta=0.0;
   float sd;
   sd=0.0;

   float *tmp1;
   tmp1=alloc1float(lt);
   float *tmp2;
   tmp2=alloc1float(lt);
   float *tmp3;
   tmp3=alloc1float(lt);

   float **utmp;
   utmp=alloc2float(lt,nr);

   float **sem;
   sem=alloc2float(n,B);

//   usigned int seed1;        
   srand((int)time(0));
   for(ix=0;ix<n;ix++)
      theta1[ix]=(float)(rand()/float(RAND_MAX))*(theta_max-theta_min)+theta_min;
/*
   for(ix=0;ix<n;ix++)
      cout<<ix<<"  "<<theta1[ix]<<endl;
   cout<<endl;
   return 0;
*/


  for(ix=0;ix<n;ix++)
     {
        ir_beg=(int)(nr*theta1[ix]+0.5);
        ir_end=nr-ir_beg;

//        cout<<ix<<"  "<<ir_beg<<"  "<<ir_end<<endl;
//        return 0;

        for(ir=0;ir<ir_end-ir_beg+1;ir++)
          for(it=0;it<lt;it++)
             utmp[ir][it]=u[ir+ir_beg][it];

        for(it=0;it<lt;it++)
          {
             tmp1[it]=0.0;
             tmp2[it]=0.0;
             tmp3[it]=0.0;
          } 

        for(it=0;it<lt;it++)
          {
             for(ir=0;ir<ir_end-ir_beg+1;ir++)
                {
                   tmp1[it]+=utmp[ir][it];
                   tmp2[it]+=utmp[ir][it]*utmp[ir][it];
                }
             tmp1[it]=tmp1[it]*tmp1[it];
             if(tmp2[it]==0)
                tmp2[it]=tmp2[it]+10000;
             else
                tmp2[it]=tmp2[it];
             tmp3[it]=tmp1[it]/tmp2[it]/(ir_end-ir_beg+1);
          }

       for(it=0;it<lt;it++)
          sem[0][ix]+=tmp3[it];     
       sem[0][ix]/=lt;    

       for(ib=1;ib<B;ib++)
         {
           for(it=0;it<lt;it++)
             {
               tmp1[it]=0.0;
               tmp2[it]=0.0;
               tmp3[it]=0.0;
             }
           for(ir=0;ir<ir_end-ir_beg+1;ir++)
             {
//               srand((int)time(0));
//               ir_tmp=(int)(rand()/int(RAND_MAX))*(ir_end-ir_beg);
               ir_tmp=rand()%(ir_end-ir_beg)+ir_beg;
               
//               if(ib%5==0)
//                 cout<<ib<<"  "<<ir_beg<<"  "<<ir_end<<"  "<<ir<<"  "<<ir_tmp<<endl;

               for(it=0;it<lt;it++)
                  utmp[ir][it]=u[ir_tmp][it];
             } 

           for(it=0;it<lt;it++)
             {
              for(ir=0;ir<ir_end-ir_beg+1;ir++)
                {
                   tmp1[it]+=utmp[ir][it];
                   tmp2[it]+=utmp[ir][it]*utmp[ir][it];
                }

              tmp1[it]=tmp1[it]*tmp1[it];
              if(tmp2[it]==0)
                tmp2[it]=tmp2[it]+10000;
              else
                tmp2[it]=tmp2[it];
              tmp3[it]=tmp1[it]/tmp2[it]/(ir_end-ir_beg+1);
             }

           for(it=0;it<lt;it++)
              sem[ib][ix]+=tmp3[it];
           sem[ib][ix]/=lt;
         } 
     } 
    
    for(ib=0;ib<B;ib++)
      {
         for(ir=0;ir<n-1;ir++)
           { 
              if(sem[ib][ir+1]>sem[ib][ir])
                  theta2[ib]=theta1[ir+1];
              else
                  theta2[ib]=theta1[ir];      
           }
      } 

   cout<<"Thetas for each realization of Bootstrap is===="<<endl;
   for(ib=0;ib<B;ib++)
     cout<<ib<<"  "<<theta2[ib]<<endl;
   
   for(ib=0;ib<B;ib++)  
      theta+=theta2[ib];
   theta/=B;

   for(ib=0;ib<B;ib++)
      sd+=(theta-theta2[ib])*(theta-theta2[ib]);
   sd=sqrt(sd)/B;

   cout<<"Average and Standard Error is===="<<endl;
   cout<<theta<<"  "<<sd<<endl;

   for(ir=0;ir<(int)(nr*theta+0.5);ir++)
     for(it=0;it<lt;it++)
        u[ir][it]=0.0;
   for(ir=(int)(nr*(1-theta)+0.5);ir<nr;ir++)
     for(it=0;it<lt;it++)
        u[ir][it]=0.0;
   for(ir=(int)(nr*theta+0.5);ir<(int)(nr*(1-theta)+0.5);ir++)
        u[ir][it]/=((int)(nr*(1-theta)+0.5)-(int)(nr*theta+0.5)); 

   return 0;
}

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256];
   int nr,lt,win,n,B,nwin;
   float theta_max,theta_min;

   ifstream swq;
   swq.open("bootstrap_for_mcg.par");
   if(!swq)
      {
         cout<<"cannot open bootstrap_for_mcg.par"<<endl;
         abort();
      }
   swq>>fn1>>fn2>>fn3>>fn4>>nr>>lt>>win>>n>>B>>theta_max>>theta_min;
   swq.close();

   nwin=int(lt/win);

   cout<<"No. of windows is===="<<win<<endl;
   cout<<"Temperal Samples for each window is===="<<nwin<<endl;

   int is,ir,it,iwin;

   float **u;
   u=alloc2float(lt,nr);
 
   float *m;
   m=alloc1float(lt);

   float **u1;
   u1=alloc2float(lt,nr);

   float *m1;
   m1=alloc1float(lt);

   float **uu;
   uu=alloc2float(nwin,nr);

   ifstream swq1;
   swq1.open(fn1,ios::binary);
    if(!swq1)
      {
         cout<<"cannot open "<<fn1<<endl;
         abort();
      }
   for(ir=0;ir<nr;ir++)
     for(it=0;it<lt;it++)
       swq1.read((char*)&u[ir][it],sizeof(u[ir][it]));
   swq1.close();

   ofstream swq2;
   swq2.open(fn2,ios::binary);
    
   ofstream swq3;
   swq3.open(fn3,ios::binary);

   ofstream swq4;
   swq4.open(fn4,ios::binary);

   for(it=0;it<lt;it++)
     {
        m[it]=0.0;
        m1[it]=0.0;
     }     

   for(it=0;it<lt;it++)
     {
        for(ir=0;ir<nr;ir++)
          m[it]+=u[ir][it];
     }    

    for(it=0;it<lt;it++)
      swq3.write((char*)&m[it],sizeof(m[it])); 
    swq3.close();


   return 0;
//    cout<<"2222"<<endl;

    if(lt%win==0)
      {
         for(iwin=0;iwin<win;iwin++)
            {
               cout<<iwin<<" th window begins..."<<endl;
         
               for(ir=0;ir<nr;ir++)
                 for(it=0;it<nwin;it++)
                    uu[ir][it]=u[ir][iwin*nwin+it];

//               cout<<nr<<"  "<<nwin<<"  "<<n<<"  "<<B<<"  "<<theta_max<<"  "<<theta_min<<endl;


               bootstrap_for_mcg(uu, nr,  nwin,  n, B, theta_max, theta_min); 

               for(ir=0;ir<nr;ir++)
                 for(it=iwin*nwin;it<iwin*nwin+nwin;it++)
                     u1[ir][it]=uu[ir][it-iwin*nwin];

               cout<<iwin<<" th window done!"<<endl;
            }             
      }
    else
      {
         for(iwin=0;iwin<win;iwin++)
            {
               cout<<iwin<<" th window begins..."<<endl;
               
               for(ir=0;ir<nr;ir++)
                 for(it=0;it<nwin;it++)
                    uu[ir][it]=u[ir][iwin*nwin+it];

               bootstrap_for_mcg(uu, nr,  nwin,  n, B, theta_max, theta_min);                   

               for(ir=0;ir<nr;ir++)
                 for(it=iwin*nwin;it<iwin*nwin+nwin;it++)
                     u1[ir][it]=uu[ir][it-iwin*nwin];
               cout<<iwin<<" th window done!"<<endl;
            }
         for(ir=0;ir<nr;ir++)
           for(it=win*nwin;it<lt;it++)
             u1[ir][it]=u[ir][it];   
      }

  for(it=0;it<lt;it++)
     {
        for(ir=0;ir<nr;ir++)
          m1[it]+=u1[ir][it];
     }

   for(ir=0;ir<nr;ir++)
     for(it=0;it<lt;it++)
       swq2.write((char*)&u1[ir][it],sizeof(u1[ir][it]));
   swq2.close();

  for(it=0;it<lt;it++)
      swq4.write((char*)&m1[it],sizeof(m1[it]));
    swq4.close();
   
  cout<<"All Done!"<<endl;
  
  return 0;

}



























































