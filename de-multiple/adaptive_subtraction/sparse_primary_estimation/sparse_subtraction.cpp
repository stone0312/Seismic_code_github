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

  ifstream swq;
  swq.open("sparse_subtraction.par");
  swq>>input1>>input2>>input3>>output1>>output2>>output3>>output4>>ns>>nr>>lt>>np;
  swq.close();

  float *win;
  win=alloc1float(nr);

  for(ir=0;ir<nr;ir++)
//     win[ir]=(int)(30*sqrt(sin(PAI/2*(nr-1-ir)/(nr-1))));
     win[ir]=(int)(30*exp(-(float)(4*ir)/(float)(nr-1)));

  for(ir=0;ir<nr;ir++)
    win[ir]=50;
 
  for(ir=0;ir<nr;ir++)
     cout<<ir+1<<" , "<<win[ir]<<endl;
//  return 0;

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

  ifstream swq5;
  swq5.open(input3,ios::binary);
  if(!swq5)
   {
     cout<<"Cannot Open "<<input3<<endl;
     return 0;
   }

  ofstream swq3;
  swq3.open(output1,ios::binary);
  if(!swq3)
   {
     cout<<"Cannot Open "<<output1<<endl;
     return 0;
   }

  ofstream swq4;
  swq4.open(output2,ios::binary);
  if(!swq4)
   {
     cout<<"Cannot Open "<<output2<<endl;
     return 0;
   }

  ofstream swq6;
  swq6.open(output3,ios::binary);
  if(!swq6)
   {
     cout<<"Cannot Open "<<output3<<endl;
     return 0;
   }

  ofstream swq7;
  swq7.open(output4,ios::binary);
  if(!swq7)
   {
     cout<<"Cannot Open "<<output4<<endl;
     return 0;
   }

  float max,maxraw,semmax;
  int t_idx,p_idx;

  float **data_tp;
  data_tp=alloc2float(lt,np);
  float **sem;
  sem=alloc2float(lt,np);
  float **semraw;
  semraw=alloc2float(lt,np);
  float **filter;
  filter=alloc2float(lt,np);
  float **filterraw;
  filterraw=alloc2float(lt,np);
  float **pri_tp;
  pri_tp=alloc2float(lt,np);
  float **pri_t;
  pri_t=alloc2float(lt,nr);
  

  for(is=0;is<ns;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
          pri_t[ir][it]=0.0;

      for(ir=0;ir<nr;ir++)
        {
           max=0.0;
           for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
              {
                swq1.read((char*)&data_tp[ip][it],sizeof(data_tp[ip][it])); 
                swq2.read((char*)&sem[ip][it],sizeof(sem[ip][it])); 
                swq5.read((char*)&semraw[ip][it],sizeof(semraw[ip][it])); 
              }

          for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               {
                  if(fabs(sem[ip][it]>max))
                    max=fabs(sem[ip][it]);
                  else
                    max=max;
               }
           if(max!=0.0)
            {
             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                sem[ip][it]=fabs(sem[ip][it])/max;
            }
          else
           {
             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                sem[ip][it]=0.0;
           }

          for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               swq3.write((char*)&sem[ip][it],sizeof(sem[ip][it]));  
/*
          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
               filter[ip][it]=1-sem[ip][it];    
*/
 
          maxraw=0.0;
          for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               {
                  if(fabs(semraw[ip][it]>maxraw))
                    maxraw=fabs(semraw[ip][it]);
                  else
                    maxraw=maxraw;
               }
           if(maxraw!=0.0)
            {
             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                semraw[ip][it]=fabs(semraw[ip][it])/maxraw;
            }
          else
           {
             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                semraw[ip][it]=0.0;
           }

          for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               swq6.write((char*)&semraw[ip][it],sizeof(semraw[ip][it]));  

           for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
                {
                  if(semraw[ip][it]>0.0000001)
                    semraw[ip][it]=1.0;
                  else
                    semraw[ip][it]=0.0;
                }

          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
               data_tp[ip][it]*=semraw[ip][it];

           for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
                swq7.write((char*)&data_tp[ip][it],sizeof(data_tp[ip][it]));


          t_idx=0;
          p_idx=0;
          semmax=0.0;

          for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
               {
                 if(sem[ip][it]>semmax)
                   {
                     semmax=sem[ip][it];
                     t_idx=it;
                     p_idx=ip;
                   }
               } 

          cout<<ir+1<<" trace , "<<"( "<<p_idx<<" , "<<t_idx<<" )"<<endl;   

            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                   filter[ip][it]=1.0;

           for(ip=p_idx-3;ip<p_idx+3;ip++)
             for(it=t_idx-win[ir];it<t_idx+20;it++)
               filter[ip][it]=0.0;

           for(ip=0;ip<np;ip++)
            for(it=0;it<lt;it++)
               pri_tp[ip][it]=data_tp[ip][it]*filter[ip][it];

           for(it=0;it<lt;it++)
            {
             for(ip=0;ip<np;ip++)
              pri_t[ir][it]+=pri_tp[ip][it];
            }

            
            for(it=0;it<lt;it++)
              swq4.write((char*)&pri_t[ir][it],sizeof(pri_t[ir][it]));

        }
    }  


  return 0; 

}
