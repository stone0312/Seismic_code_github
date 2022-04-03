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
   int ns,nr,np,lt,tnum;
   float thre;

   ifstream swq;
   swq.open("optimize_p.par");
   swq>>fn1>>fn2>>ns>>nr>>np>>lt>>tnum;
   swq.close();

   int is, ir, ix, ip, it, it1;
   int pmax;

   float **sem;
   sem=alloc2float(lt,np);

   float *amp1;
   amp1=alloc1float(np);

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

             for(it=0;it<tnum;it++)
                {
                  for(ip=0;ip<np;ip++)
                     amp1[ip]=0.0;

                  for(ip=0;ip<np;ip++)
                     for(it1=it*lt/tnum;it1<(it+1)*lt/tnum;it1++)               
                         sem_new[ip][it1]=0.0;

                  for(ip=0;ip<np;ip++)
                     {
                       for(it1=it*lt/tnum;it1<(it+1)*lt/tnum;it1++)               
                         amp1[ip]+=sem[ip][it1]*sem[ip][it1];
                     }

                   pmax=0;
                   for(ip=0;ip<np-1;ip++)
                      {
                         if(amp1[ip+1]>amp1[pmax])
                           pmax=ip+1;
                      }

                   if(amp1[pmax]!=0.0)
                     {

                       if(pmax>0&&pmax<np)
                        {
                          cout<<it+1<<" time window, pmax===="<<pmax<<endl;

                          for(it1=it*lt/tnum;it1<(it+1)*lt/tnum;it1++)
                             sem_new[pmax][it1]=sem[pmax][it1];
/*
                          for(ip=pmax-1;ip<pmax+2;ip++)
                            for(it1=it*lt/tnum;it1<(it+1)*lt/tnum;it1++)
                             sem_new[ip][it1]=sem[ip][it1];
*/
                        }

                     }

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





















