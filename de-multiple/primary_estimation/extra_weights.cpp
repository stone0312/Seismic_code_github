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
   int ns, nr, ltt , t0, r0;
   
   int is,ir,it;

   ifstream swq1;
   swq1.open("extra_weights.par");
   swq1>>fn1>>fn2>>fn3>>fn4>>fn5>>ns>>nr>>ltt>>t0>>r0;
   swq1.close();

   float **w1;
   w1=alloc2float(ltt,nr); 

   float **w2;
   w2=alloc2float(nr,ns); 

   float *w3;
   w3=alloc1float(nr); 

   float *w4;
   w4=alloc1float(nr); //smoothing results of w3

   float *gauss;
   gauss=alloc1float(15); 

  for(ir=0;ir<15;ir++)
   gauss[ir]=0.5*(1-cos(2*pai*ir/14));

  for(ir=0;ir<15;ir++)
   cout<<gauss[ir]<<endl;

   float *w5;
   w5=alloc1float(nr); //final results 
   float *w6;
   w6=alloc1float(nr); //final results 
   float *w7;
   w7=alloc1float(nr); //final results 

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

   ofstream swq5;
   swq5.open(fn4,ios::binary);
 
   ofstream swq6;
   swq6.open(fn5,ios::binary);

   for(is=0;is<ns;is++)
     {
        for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
           swq2.read((char*)&w1[ir][it],sizeof(w1[ir][it]));

        for(ir=0;ir<nr;ir++)
          w2[is][ir]=w1[ir][t0];
     }

   for(is=0;is<ns;is++)
    {
      for(ir=0;ir<nr;ir++)
       swq3.write((char*)&w2[is][ir],sizeof(w2[is][ir]));
    }

  for(ir=0;ir<nr;ir++)
     w3[ir]=w2[r0][ir];  

  for(ir=0;ir<nr;ir++)
   swq4.write((char*)&w3[ir],sizeof(w3[ir]));

  for(ir=0;ir<nr;ir++)
   w4[ir]=0.0; 

  float max=0.0;
  for(ir=0;ir<nr;ir++)
    {
       if(w3[ir]>max)
         max=w3[ir];
    } 

  for(ir=0;ir<nr;ir++)
   w4[ir]=w3[ir]/max;

  w4[18]=(w4[17]+w4[19])/2.0;
  w4[21]=(w4[20]+w4[22])/2.0;
  w4[23]=(w4[22]+w4[24])/2.0;
  w4[26]=(w4[25]+w4[27])/2.0;
  w4[27]=(w4[26]+w4[28])/2.0;
  w4[31]=(w4[30]+w4[32])/2.0;
  w4[33]=(w4[32]+w4[34])/2.0;
  w4[35]=(w4[34]+w4[36])/2.0;
  w4[26]=(w4[25]+w4[27])/2.0;
  w4[25]=(w4[24]+w4[26])/2.0;
  w4[27]=(w4[26]+w4[28])/2.0;
  w4[26]=(w4[25]+w4[27])/2.0;
  w4[12]=(w4[11]+w4[14])*1/3.0;
  w4[13]=(w4[11]+w4[14])*2/3.0;
/*
  w4[28]=(w4[27]+w4[29])/2.0;
  w4[30]=(w4[29]+w4[31])/2.0;
  w4[31]=(w4[30]+w4[32])/2.0;
  w4[33]=(w4[32]+w4[34])/2.0;
  w4[35]=(w4[34]+w4[36])/2.0;

  w4[30]=(w4[29]+w4[31])/2.0;
  w4[31]=(w4[30]+w4[32])/2.0;
  w4[32]=(w4[31]+w4[33])/2.0;

  w4[16]=(w4[15]+w4[17])/2.0;
  w4[18]=(w4[17]+w4[19])/2.0;
  w4[20]=(w4[19]+w4[21])/2.0;
  w4[22]=(w4[21]+w4[23])/2.0;
  w4[24]=(w4[23]+w4[25])/2.0;

  w4[19]*=6.0;
  w4[18]*=4.0;
  w4[20]*=4.0;
  w4[17]*=3.0;
  w4[21]*=3.0;
  w4[24]*=0.0;
  w4[13]*=0.1;

  for(ir=37;ir<49;ir++)
    w4[ir]=w4[72-ir];
*/


/*
  w4[64]=(w4[63]+w4[65])/2.0;
  w4[66]=(w4[65]+w4[67])/2.0;
  w4[68]=(w4[67]+w4[69])/2.0;
  w4[71]=(w4[70]+w4[72])/2.0;
  w4[73]=(w4[72]+w4[74])/2.0;
  w4[75]=(w4[74]+w4[76])/2.0;
  w4[77]=(w4[76]+w4[78])/2.0;
  w4[81]=(w4[80]+w4[82])/2.0;
  w4[83]=(w4[82]+w4[84])/2.0;
  w4[85]=(w4[84]+w4[86])/2.0;

  for(ir=49;ir<62;ir++)
    w4[ir]=w4[124-ir];
  for(ir=85;ir<94;ir++)
    w4[ir]=w4[168-ir];
*/
  for(ir=0;ir<nr;ir++)
   swq5.write((char*)&w4[ir],sizeof(w4[ir]));   

  for(ir=0;ir<nr;ir++)
   {
     w5[ir]=0.0; 
     w6[ir]=0.0; 
     w7[ir]=0.0; 
   }

  for(ir=13;ir<37;ir++)
    w5[ir]=w4[ir+11];
  for(ir=28;ir<42;ir++)
    w6[ir]=w4[ir-4-2*(ir-28)];
 
  for(ir=0;ir<nr;ir++)
   w7[ir]=w5[ir]+w6[ir]; 
/*
  w7[16]-=0.2;
  w7[17]-=0.2;
  w7[18]-=0.2;
  w7[19]-=0.2;
  w7[19]=(w7[18]+w7[20])/2.0;
*/


/*
  w7[68]/=2.0;
  w7[69]/=2.0;
  w7[70]/=2.0;
  w7[71]/=2.0;

  w7[68]=(w7[67]+w7[69])/2.0;
  w7[69]=(w7[68]+w7[70])/2.0;
  w7[70]=(w7[69]+w7[71])/2.0;
*/
  max=0.0;
  for(ir=0;ir<nr;ir++)
    {
       if(w7[ir]>max)
         max=w7[ir];
    } 

  for(ir=0;ir<nr;ir++)
   w7[ir]=w7[ir]/max;

  for(ir=0;ir<nr;ir++)
   swq6.write((char*)&w7[ir],sizeof(w7[ir]));   

/*
  for(ir=0;ir<nr;ir++)
   w5[ir]=0.0; 

  w5[74]=w4[62];
  w5[64]=w4[84];

  for(ir=67;ir<82;ir++)
    w5[ir]=w5[74]*gauss[ir-67];

  for(ir=57;ir<72;ir++)
    w6[ir]=w5[64]*gauss[ir-57];

  for(ir=0;ir<nr;ir++)
   w7[ir]=w5[ir]+w6[ir]; 

  for(ir=0;ir<nr;ir++)
   swq6.write((char*)&w7[ir],sizeof(w7[ir]));   
*/
  return 0;

}






















