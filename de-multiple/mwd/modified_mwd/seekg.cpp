#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"


int main()
{
   int ns=4;
   int lt=3;

   int is,ir,it;

   float **d;
   d=alloc2float(lt,ns);

   float *dd;
   dd=alloc1float(lt*ns);
 
   d[0][0]=1;
   d[0][1]=2;
   d[0][2]=3;
   d[1][0]=4;
   d[1][1]=5;
   d[1][2]=6;
   d[2][0]=7;
   d[2][1]=8;
   d[2][2]=9;
   d[3][0]=10;
   d[3][1]=11; 
   d[3][2]=12;

   cout<<"Original Data is====="<<endl;
   for(is=0;is<ns;is++) 
      {
         for(it=0;it<lt;it++)
             cout<<d[is][it]<<"   ";
         cout<<endl;
      } 

   ofstream swq;
   swq.open("11.dat",ios::binary);
   for(is=0;is<ns;is++)
     for(it=0;it<lt;it++)
        swq.write((char*)&d[is][it],sizeof(d[is][it]));
   swq.close();

   ifstream swq1;
   swq1.open("11.dat",ios::binary);

   swq1.seekg(1*4,ios::beg);
   streampos sp=swq1.tellg();
   cout<<"0  "<<sp<<endl;

   swq1.read((char*)&dd[0],sizeof(dd[0]));

   for(is=1;is<ns;is++)
      {
         swq1.seekg((lt-1)*4,ios::cur);
 
         streampos sp=swq1.tellg();
         cout<<is<<"  "<<sp<<endl; 

         swq1.read((char*)&dd[is],sizeof(dd[is]));
      }

   for(is=0;is<ns;is++)
      cout<<dd[is]<<"  ";
   cout<<endl;

   swq1.close();
 
    return 0;

}
