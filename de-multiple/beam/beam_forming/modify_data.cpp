

using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"

int main()
{
   char fn1[256],fn2[256],fn3[256];
   int ns,nr,lt;
   int thre;
   float amp;

   ifstream swq;
   swq.open("modify_data.par");
   swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>thre;
   swq.close();

   cout<<thre<<endl;


   int ir,is,it;
   float *u1=new float [lt];
   float *u2=new float [lt];

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
      {
         cout<<"cannot open "<<fn1<<endl;
         return 0;
      }
  
   ifstream swq2;
   swq2.open(fn2,ios::binary);
   if(!swq2)
      {
         cout<<"cannot open "<<fn2<<endl;
         return 0;
      }

   ofstream swq3;
   swq3.open(fn3,ios::binary);
   if(!swq3)
      {
         cout<<"cannot open "<<fn3<<endl;
         return 0;
      }

   for(is=0;is<ns*nr;is++)
     {
         for(it=0;it<lt;it++)
           swq1.read((char*)&u1[it],sizeof(u1[it]));
         for(it=0;it<lt;it++)
           swq2.read((char*)&u2[it],sizeof(u2[it]));
/*
         for(it=0;it<lt;it++)
           {
              if(u2[it]<pow(10,thre))
                 u1[it]=0.0;
           } 
*/
         amp=0.0;
         for(it=0;it<lt;it++)
            amp+=u1[it]*u1[it];

         if(amp==0.0)
           {
             for(it=0;it<lt;it++)
               u2[it]=0.0;
           }
         else
          {
           for(it=0;it<lt;it++)
               u2[it]*=(u1[it]*u1[it])/amp;
          }
         for(it=0;it<lt;it++)
           swq3.write((char*)&u2[it],sizeof(u2[it]));

             cout<<is+1<<" trace  done!"<<endl;

      }


   swq1.close();
   swq2.close();
   swq3.close();
   
   cout<<"done!"<<endl;
   return 0;
   
   
}
