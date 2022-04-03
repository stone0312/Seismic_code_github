#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"

int main()
{
   int nx,lt,nx1;
   int is,ix,it,tmp;
   int ns=479;
   nx=479;
   lt=3000;
   nx2=nx*2-1;

   float **u1;
   u1=alloc2float(lt,nx);

   float **u2;
   u2=alloc2float(lt,nx);

   ifstream swq1;
   swq1.open("/data2/swq/479shot_raw.dat",ios::binary);

   ofstream swq2;
   swq2.open("/data2/swq/new_479shot_raw.dat",ios::binary);

   for(is=0;is<ns;is++)
     {
        for(ix=0;ix<nx;ix++)
         {
          for(it=0;it<lt;it++)
             swq1.read((char*)&(u1[ix][it]),sizeof(u1[ix][it]));
         }
    
        for(ix=0;ix<nx;ix++)
          {
             tmp=nx-ix-1;
             for(it=0;it<lt;it++)
               u2[ix][it]=u1[tmp][it];             
          }
 
        for(ix=0;ix<nx;ix++)
          for(it=0;it<lt;it++)
             swq2.write((char*)&(u2[ix][it]),sizeof(u2[ix][it]));
 
        cout<<is<<" shot modification done!"<<endl;

     }

   swq1.close();
   swq2.close();


   return 0;

}
