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
   char fn1[256],fn2[256];
   int ns,nr,np,lt,pw,tw;

   ifstream swq;
   swq.open("weighted_tp.par");
   swq>>fn1>>fn2>>ns>>nr>>np>>lt>>pw>>tw;
   swq.close();

   int is, ir, ix, ip, it;
   int pnum;
   float amp_all,amp_ave;

   float **sem;
   sem=alloc2float(lt,np);

   float *amp;
   amp=alloc1float(np);

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
                  sem_new[ip][it]=0.0;

             for(ip=pw;ip<np-pw;ip++)
               for(it=tw;it<lt-tw;it++)
                  sem_new[ip][it]+=sem[ip-pw][it-tw]*sem[ip+pw][it+tw]; 

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





















