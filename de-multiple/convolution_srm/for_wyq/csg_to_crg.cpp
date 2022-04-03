#include <iostream>
#include <fstream>
#include <iomanip>
#include "math.h"

using namespace std;

int main()
{
   char fn1[256],fn2[256];
   int ns,nr,lt;
   int is,ir,it;

   ifstream swq;
   swq.open("csg_to_crg.par");
   swq>>fn1>>fn2>>ns>>nr>>lt;
   swq.close();

   float **crg=new float *[ns];   
   for(is=0;is<ns;is++)
      crg[is]=new float [lt];

   ifstream swq1;
   swq1.open(fn1,ios::binary);

   ofstream swq2;
   swq2.open(fn2,ios::binary);

   for(ir=0;ir<nr;ir++)
      {
          swq1.seekg(ir*lt*4,ios::beg);
          for(it=0;it<lt;it++)
             swq1.read((char*)&crg[0][it],sizeof(crg[0][it]));

          for(is=1;is<ns;is++)
             {
                swq1.seekg((nr-1)*lt*4,ios::cur);
                for(it=0;it<lt;it++)
                   swq1.read((char*)&crg[is][it],sizeof(crg[is][it]));
             } 
          for(is=0;is<ns;is++)
             {
                for(it=0;it<lt;it++)
                    swq2.write((char*)&crg[is][it],sizeof(crg[is][it]));
             }

      }
 
    cout<<"all done!"<<endl;

    return 0;

}

