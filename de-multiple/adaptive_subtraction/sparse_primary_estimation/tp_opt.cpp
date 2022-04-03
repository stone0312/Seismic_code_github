using namespace std;
#include <iostream>
#include <fstream>
#include <string.h>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"

#ifndef PAI
#define PAI (3.141592653589793)
#endif

int inverse_tau_p(float **u,float *ufinal,int lt,int np)
{
    int i,j;

    for(i=0;i<lt;i++)
      ufinal[i]=0.0;

    for(i=0;i<lt;i++)
       {
           for(j=0;j<np;j++)
               ufinal[i]+=u[j][i];
       }
    return 0; 
}


int main()
{
  char input1[256], input2[256], output1[256], output2[256],input3[256],output3[256],output4[256]; 
  int ns,nr,lt,np;
  float thres;
  int is,ir,ip,it;

  ifstream swq;
  swq.open("tp_opt.par");
  swq>>input1>>input2>>output2>>output3>>ns>>nr>>lt>>np>>thres;
  swq.close();

  ifstream swq1;
  swq1.open(input1,ios::binary);
  if(!swq1)
   {
     cout<<"Cannot Open "<<input1<<endl;
     return 0;
   }

  ifstream swq2;
  swq2.open(input2,ios::binary);
  if(!swq2)
   {
     cout<<"Cannot Open "<<input2<<endl;
     return 0;
   }

  ofstream swq4;
  swq4.open(output2,ios::binary);
  if(!swq4)
   {
     cout<<"Cannot Open "<<output2<<endl;
     return 0;
   }

  ofstream swq5;
  swq5.open(output3,ios::binary);
  if(!swq5)
   {
     cout<<"Cannot Open "<<output3<<endl;
     return 0;
   }

  float **tpraw;
  tpraw=alloc2float(lt,np);
  float **semraw;
  semraw=alloc2float(lt,np);
  float **sem;
  sem=alloc2float(lt,np);
  float **tp;
  tp=alloc2float(lt,np);
  float *uitp;
  uitp=alloc1float(lt);

  float max; 
 
  for(is=0;is<ns;is++)
    {
      cout<<"  "<<is+1<<" Shot, "<<" Begin ..."<<endl;

      for(ir=0;ir<nr;ir++)
        {
           max=0.0;
           for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
              {
                swq1.read((char*)&tpraw[ip][it],sizeof(tpraw[ip][it])); 
                swq2.read((char*)&semraw[ip][it],sizeof(semraw[ip][it])); 
              }

         for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
             {
               if(semraw[ip][it]>thres)
                 sem[ip][it]=semraw[ip][it];
               else 
                 sem[ip][it]=0.0;
             }

          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
             {
               if(sem[ip][it]>0.0)
                 tp[ip][it]=tpraw[ip][it];
               else 
                 tp[ip][it]=0.0;
             }
           for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
                swq4.write((char*)&tp[ip][it],sizeof(tp[ip][it]));


          inverse_tau_p(tp,uitp,lt,np);

          for(it=0;it<lt;it++)
             swq5.write((char*)&uitp[it],sizeof(uitp[it]));


          cout<<"     "<<is+1<<" Shot, "<<ir+1<<" Trace Done ..."<<endl;
        }

       cout<<"  "<<is+1<<" Shot, "<<" Done!"<<endl;

    }  

     cout<<"All Done!"<<endl;

  return 0; 

}
