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

int main()
{
  char input1[256], input2[256], output1[256], output2[256],input3[256],output3[256],output4[256]; 
  int ns,nr,lt,np;
  int is,ir,ip,it;
  float thres;

  ifstream swq;
  swq.open("sparse_subtraction.par");
  swq>>input1>>input2>>input3>>output1>>output2>>ns>>nr>>lt>>np>>thres;
  swq.close();

  float *win;
  win=alloc1float(nr);

  for(ir=0;ir<nr;ir++)
//     win[ir]=(int)(30*sqrt(sin(PAI/2*(nr-1-ir)/(nr-1))));
     win[ir]=(int)(30*exp(-(float)(4*ir)/(float)(nr-1)));

  for(ir=0;ir<nr;ir++)
    win[ir]=50;
/* 
  for(ir=0;ir<nr;ir++)
     cout<<ir+1<<" , "<<win[ir]<<endl;
*/

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

  ifstream swq3;
  swq3.open(input3,ios::binary);
  if(!swq3)
   {
     cout<<"Cannot Open "<<input3<<endl;
     return 0;
   }

  ofstream swq4;
  swq4.open(output1,ios::binary);
  if(!swq4)
   {
     cout<<"Cannot Open "<<output1<<endl;
     return 0;
   }

  ofstream swq5;
  swq5.open(output2,ios::binary);
  if(!swq5)
   {
     cout<<"Cannot Open "<<output2<<endl;
     return 0;
   }

  float max,maxraw,semmax;
  int t_idx,p_idx;

  float **data_tp;
  data_tp=alloc2float(lt,np);
  float **sem_d;
  sem_d=alloc2float(lt,np);
  float **sem_m;
  sem_m=alloc2float(lt,np);
  float **sem_p;
  sem_p=alloc2float(lt,np);
  float **pri_tp;
  pri_tp=alloc2float(lt,np);
  float **pri_t;
  pri_t=alloc2float(lt,nr);

  for(is=0;is<ns;is++)
    {
      cout<<is+1<<" Shot Begin ..."<<endl;
 
     for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
          pri_t[ir][it]=0.0;

      for(ir=0;ir<nr;ir++)
        {
           for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
              {
                swq1.read((char*)&data_tp[ip][it],sizeof(data_tp[ip][it])); 
                swq2.read((char*)&sem_d[ip][it],sizeof(sem_d[ip][it])); 
                swq3.read((char*)&sem_m[ip][it],sizeof(sem_m[ip][it])); 
              }

          for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               {
                  if(sem_m[ip][it]>thres)
                    sem_p[ip][it]=1.0-sem_d[ip][it]*sem_m[ip][it];
                  else
                    sem_p[ip][it]=1.0;
               }
          for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               swq4.write((char*)&sem_p[ip][it],sizeof(sem_p[ip][it]));  

          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
               pri_tp[ip][it]=data_tp[ip][it]*sem_p[ip][it];


            for(it=0;it<lt;it++)
               {
                 for(ip=0;ip<np;ip++)
                   pri_t[ir][it]+=pri_tp[ip][it];
               }   

            for(it=0;it<lt;it++)
                swq5.write((char*)&pri_t[ir][it],sizeof(pri_t[ir][it]));

            if((ir+1)%10==0)
               cout<<is+1<<" Shot, "<<ir+1<<" Trace Done ..."<<endl;

        }

      cout<<is+1<<" Shot Done !"<<endl;

    }  


  return 0; 

}
