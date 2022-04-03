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

int main()
{
   char fn1[256],fn2[256],fn3[256];
   int ns,nr,np,lt,tw_h;
   float thre;

   ifstream swq;
   swq.open("neighbour_correlation_global.par");
   swq>>fn1>>fn2>>fn3>>ns>>nr>>np>>lt>>thre;
   swq.close();

   int is, ir, ix, ip, it, it1;
   float amp1,amp2,amp3;

   float **sem;
   sem=alloc2float(lt,np);

   float *amp;
   amp=alloc1float(np);

   float *coe;
   coe=alloc1float(np);

   float **sem_new;
   sem_new=alloc2float(lt,np);

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
      {
        cout<<"Cannot open "<<fn1<<endl;
        return 0;
      }

   ofstream swq2;
   swq2.open(fn2,ios::binary);
   if(!swq2)
      {
        cout<<"Cannot open "<<fn2<<endl;
        return 0;
      }


   for(is=0;is<ns;is++)
     {
       for(ir=0;ir<nr;ir++)
          {
             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                 swq1.read((char*)&sem[ip][it],sizeof(sem[ip][it]));

             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                  sem_new[ip][it]=sem[ip][it];

             for(ip=0;ip<np;ip++)
                coe[ip]=0.0;

             for(ip=0;ip<np-1;ip++)
                {

                           amp1=0.0;
                   	   amp2=0.0;
                   	   amp3=0.0;

                     for(it=0;it<lt;it++) 
                        {
                           amp1+=sem[ip][it]*sem[ip+1][it];
                           amp2+=sem[ip][it]*sem[ip][it]; 
                           amp3+=sem[ip+1][it]*sem[ip+1][it]; 

                           if(amp2!=0.0&&amp3!=0.0)
                             coe[ip]=amp1/sqrt(amp2*amp3);
                           else
                             coe[ip]=0.0;
                        }
                }

            for(ip=0;ip<np-1;ip++)
               cout<<ip+1<<"   "<<coe[ip]<<endl;


             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++) 
                  {
                    if(coe[ip]<thre)
                       sem_new[ip][it]=0.0;
                  }


             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                 swq2.write((char*)&sem_new[ip][it],sizeof(sem_new[ip][it]));

             cout<<ir+1<<" trace average filtering done..."<<endl;

          }
       cout<<is+1<<" shot average filtering done..."<<endl;

     }

   cout<<"All  Done!"<<endl;

   return 0;

}





















