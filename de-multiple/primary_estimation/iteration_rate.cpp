using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#define pai 3.14159265

int main()
{
   char fn1[256],fn2[256];
   int iter,tmax;

   ifstream swq;
   swq.open("iteration_rate.par");
   swq>>fn1>>fn2>>tmax;
   swq.close();

   float *w;
   w=alloc1float(tmax);

   ifstream swq1;
   swq1.open(fn1);
  
   for(iter=0;iter<tmax;iter++)
     swq1>>w[iter];

   for(iter=0;iter<tmax;iter++)
     cout<<iter+1<<" ,  "<<w[iter]<<endl;

   ofstream swq2;
   swq2.open(fn2,ios::binary);
   
   for(iter=0;iter<tmax;iter++)
     swq2.write((char*)&w[iter],sizeof(w[iter]));

   return 0;
}
