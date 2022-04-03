using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"

int main()
{
  char fn1[256],fn2[256],fn3[256],fn4[256];
  int ns,nr,np,lt;
  float thres;

  int is,ir,ip,it;

  ifstream swq;
  swq.open("primary_estimation_forward_tp_opt.par");
  swq>>fn1>>fn2>>fn3>>fn4>>ns>>nr>>np>>lt>>thres;
  swq.close();

  cout<<"Fna of CSG Real is===="<<fn1<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<endl;

  float **usf;
  usf=alloc2float(lt,np);
  float **sem;
  sem=alloc2float(lt,np);
  float **multp;
  multp=alloc2float(lt,np);

  float * mult;
  mult=alloc1float(lt); 

  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl; 

  ifstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;

  ofstream swq3;
  swq3.open(fn3,ios::binary);
  if(!swq3)
       cout<<"cannot open "<<fn3<<endl;

  ofstream swq4;
  swq4.open(fn4,ios::binary);
  if(!swq4)
       cout<<"cannot open "<<fn4<<endl;

  for(is=0;is<ns;is++)
//  for(is=0;is<1;is++)
    {
       for(ir=0;ir<nr;ir++)
         {
            for(it=0;it<lt;it++) 
              mult[it]=0.0;

            for(ip=0;ip<np;ip++) 
              for(it=0;it<lt;it++)
                {
                  swq1.read((char*)&usf[ip][it],sizeof(usf[ip][it])); 
                  swq2.read((char*)&sem[ip][it],sizeof(sem[ip][it])); 
                }
            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                {
                   if(sem[ip][it]<thres)
                     multp[ip][it]=0.0;
                   else
                     multp[ip][it]=usf[ip][it];
                }
            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
                  swq3.write((char*)&multp[ip][it],sizeof(multp[ip][it])); 


            for(it=0;it<lt;it++)
              {
                 for(ip=0;ip<np;ip++)
                    mult[it]+=multp[ip][it];
              }
              for(it=0;it<lt;it++)
                  swq4.write((char*)&mult[it],sizeof(mult[it])); 

         }

      cout<<is+1<<" Shot Optimization Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











