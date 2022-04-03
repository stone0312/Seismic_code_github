using namespace std;

#include <iostream>
#include "math.h"
#include "stdio.h"
#include <fstream>
#include "stdlib.h"
#include <complex>
#include "alloc.c"

#define pai 3.14159265

int soft_thres_hrtp(float **hrtp, float **tp, float **sem, int np, int lt, int pw_half, int tw_half, int thalf,  float threshold, int ite_max, float *window)
{
   int ip, it, ite;
   float sem_max;
  
   int p_idx, t_idx;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
        hrtp[ip][it]=0.0;

   for(ite=0;ite<ite_max;ite++)
     {
        sem_max=0.0;

        for(ip=0;ip<np;ip++)
          for(it=0;it<lt;it++)
            {
              if(sem[ip][it]>sem_max)
                {
                   sem_max=sem[ip][it];
                   p_idx=ip;
                   t_idx=it;
                }
            }

         cout<<"Iteration Time is ==== "<<ite+1<<" , "<<"Maximum Semblance is ==== "<<sem_max<<endl;
         cout<<"Tau and P Index are ==== "<<"( "<<t_idx<<" , "<<p_idx<<" ) "<<endl;

         if(sem_max<threshold)
           ite=ite_max;

         else 
           {
             for(it=t_idx-thalf;it<t_idx+thalf+1;it++)
                hrtp[p_idx][it]=tp[p_idx][it]*window[it-t_idx+thalf];

            if(p_idx-pw_half<0)
              {
                for(ip=0;ip<p_idx+pw_half+1;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem[ip][it]=0.0;
              }

            else if(p_idx+pw_half+1>np)
              {
                for(ip=p_idx-pw_half;ip<np;ip++)
                 for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                    sem[ip][it]=0.0;
              }
            else
             {
              for(ip=p_idx-pw_half;ip<p_idx+pw_half+1;ip++)
               for(it=t_idx-tw_half;it<t_idx+tw_half+1;it++)
                  sem[ip][it]=0.0;
             }

           }

     }
   
    return 0;

}

int main()
{
   char input1[256], input2[256], output1[256];
   int np, lt, pw_half, tw_half, thalf, ite_max, taper;
   float threshold;
  
   int ip, it;
 
   ifstream swq;
   swq.open("soft_thres_hrtp.par");

   swq>>input1>>input2>>output1>>np>>lt>>pw_half>>tw_half>>thalf>>ite_max>>threshold>>taper;   
   swq.close();

   float **tp=new float *[np];
   for(ip=0;ip<np;ip++)
      tp[ip]=new float [lt]; 

   float **sem=new float *[np];
   for(ip=0;ip<np;ip++)
      sem[ip]=new float [lt]; 

   float **hrtp=new float *[np];
   for(ip=0;ip<np;ip++)
      hrtp[ip]=new float [lt]; 

   float *window=new float [2*thalf+1];
   for(it=0;it<2*thalf+1;it++)
      window[it]=1.0;
   for(it=0;it<taper;it++)
      window[it]=sqrt(sin(pai*it/(2*taper)));
   for(it=2*thalf+1-taper;it<2*thalf+1;it++)
      window[it]=window[2*thalf-it];

    for(it=0;it<2*thalf+1;it++)
      cout<<it<<"   "<<window[it]<<endl;
//    return 0;

   ifstream swq1;
   swq1.open(input1,ios::binary);
   if(!swq1)
     cout<<"cannot open"<<input1<<endl; 

   ifstream swq2;
   swq2.open(input2,ios::binary);
   if(!swq2)
     cout<<"cannot open"<<input2<<endl; 
 
   ofstream swq3;
   swq3.open(output1,ios::binary);
   if(!swq3)
     cout<<"cannot open"<<output1<<endl; 

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       {
           swq1.read((char*)&tp[ip][it],sizeof(tp[ip][it]));
           swq2.read((char*)&sem[ip][it],sizeof(sem[ip][it]));
       } 

   soft_thres_hrtp(hrtp, tp, sem, np, lt, pw_half, tw_half, thalf,  threshold, ite_max, window);

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
      swq3.write((char*)&hrtp[ip][it],sizeof(hrtp[ip][it]));

   return 0;

}













































