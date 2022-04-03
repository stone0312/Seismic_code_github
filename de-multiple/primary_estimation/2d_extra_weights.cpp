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
   char fn1[256],fn2[256],fn3[256], fn4[256],fn5[256];
   int ns, nr, ltt , np,lmcg,  r0, t0;
   
   int is,ir,it;

   ifstream swq1;
   swq1.open("2d_extra_weights.par");
   swq1>>fn1>>fn2>>fn3>>ns>>nr>>ltt>>np>>lmcg>>r0>>t0;
   swq1.close();

   float *w1;
   w1=alloc1float(lmcg*nr*np); 

   float *w2;
   w2=alloc1float(lmcg*np); 

   float **w3;
   w3=alloc2float(np,lmcg); 

   float **w4;
   w4=alloc2float(np,lmcg); 

   float **w5;
   w5=alloc2float(np,lmcg); 

   ifstream swq2;
   swq2.open(fn1,ios::binary);
   if(!swq2)
     {
       cout<<"Cannot Open "<<fn1<<endl;
       return 0;
     } 

   ofstream swq3;
   swq3.open(fn2,ios::binary);

   ofstream swq4;
   swq4.open(fn3,ios::binary);

   swq2.seekg(0,ios::beg);
   swq2.seekg(t0*4,ios::cur);
     swq2.read((char*)&w1[0],sizeof(w1[0]));

   for(ir=1;ir<lmcg*nr*np;ir++)
    {
     swq2.seekg(ltt*4,ios::cur);
     swq2.read((char*)&w1[ir],sizeof(w1[ir]));
    }

   for(ir=0;ir<lmcg*np;ir++)
     w2[ir]=w1[r0*lmcg*np+ir];

   for(ir=0;ir<lmcg;ir++)    
     for(it=0;it<np;it++)
       w3[ir][it]=w2[ir*np+it];

   float max=0.0;
   for(ir=0;ir<lmcg;ir++)
     for(it=0;it<np;it++)     
       {
         if(w3[ir][it]>max)
           max=w3[ir][it];
       }

   cout<<"Maximum ==== "<<max<<endl;

   if(max!=0.0)
    {
     for(ir=0;ir<lmcg;ir++)
       for(it=0;it<np;it++)     
       w3[ir][it]/=max;
    }
   else
    for(ir=0;ir<lmcg;ir++)
       for(it=0;it<np;it++)
       w3[ir][it]=0.0;

   for(ir=0;ir<lmcg;ir++)    
     for(it=0;it<np;it++)
     swq3.write((char*)&w3[ir][it],sizeof(w3[ir][it]));

   for(ir=0;ir<lmcg;ir++)    
     for(it=0;it<np;it++)
       w4[ir][it]=0.0;

   for(it=3;it<6;it++)
     w4[20][it]=w3[42][it+1];

//   w4[80][10]=(w4[80][9]+w4[80][11])/2.0;
   for(it=3;it<6;it++)
     {
       w4[21][it]=w4[20][it]/2.0;
       w4[22][it]=w4[21][it]/2.0;
       w4[23][it]=w4[22][it]/2.0;
       w4[19][it]=w4[20][it]/2.0;
       w4[18][it]=w4[19][it]/2.0;
       w4[17][it]=w4[18][it]/2.0;
     }
      
   for(it=2;it<5;it++)
     w4[28][it]=w3[42][it+2];  
 
//   w4[34][11]=(w4[34][10]+w4[34][12])/2.0;

   for(it=2;it<5;it++)
     {
       w4[29][it]=w4[28][it]/2.0;
       w4[30][it]=w4[29][it]/2.0;
       w4[31][it]=w4[30][it]/2.0;
       w4[27][it]=w4[28][it]/2.0;
       w4[26][it]=w4[27][it]/2.0;
       w4[25][it]=w4[26][it]/2.0;
     }

   for(it=1;it<4;it++)
     w4[15][it]=w3[42][it+3];  
 
//   w4[34][11]=(w4[34][10]+w4[34][12])/2.0;

   for(it=1;it<4;it++)
     {
       w4[16][it]=w4[15][it]/2.0;
       w4[17][it]=w4[16][it]/2.0;
       w4[18][it]=w4[17][it]/2.0;
       w4[14][it]=w4[15][it]/2.0;
       w4[13][it]=w4[14][it]/2.0;
       w4[12][it]=w4[13][it]/2.0;
     }
      
   for(ir=0;ir<lmcg;ir++)    
     for(it=0;it<np;it++)
     swq4.write((char*)&w4[ir][it],sizeof(w4[ir][it]));

  return 0;

}






















