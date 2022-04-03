#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"


int main()
{
   int ns=500;
   int nr=500;
   int ltt=10000;
   int nsr=ns*nr;

   int is,ir,it;
   int itmp;

   it=100;

   float *d;
   d=alloc1float(nsr);

   float **dd;
   dd=alloc2float(nr,ns);

   ifstream swq1;
   swq1.open("/data2/swq/2d_modeling_data/fft_orig_shots_ltt_real.dat",ios::binary);
   if(!swq1)
     {
       cout<<"cannot open "<<endl;
       abort();
     }
    
   ofstream swq2;
   swq2.open("/data2/swq/2d_modeling_data/test_seekg_fft_orig_shots_one_frequency_real.dat",ios::binary);
   if(!swq2)
     {
       cout<<"cannot open "<<endl;
       abort();
     }
 
   swq1.seekg(it*4,ios::beg);
   swq1.read((char*)&(d[0]),sizeof(d[0]));
   
   streampos sp=swq1.tellg();
   cout<<"0 "<<sp<<endl;

   for(is=1;is<nsr;is++)
      {
         swq1.seekg((ltt-1)*4,ios::cur);

         streampos sp=swq1.tellg();
         cout<<is<<"    "<<sp<<endl;

         swq1.read((char*)&(d[is]),sizeof(d[is]));
//         if((is+1)%nr==0)
//            cout<<(is+1)/nr<<"th shot done!"<<endl;
      }

   for(is=0;is<ns;is++)
     {
       for(ir=0;ir<nr;ir++)
         {
            itmp=is*nr+ir;
            dd[is][ir]=d[itmp];
         }
     }

     for(is=0;is<ns;is++)
         for(ir=0;ir<nr;ir++)
           swq2.write((char*)&(dd[is][ir]),sizeof(dd[is][ir]));
     cout<<"ALL DONE!"<<endl;


     swq1.close();
     swq2.close();
 
    return 0;

}
