#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#define pai 3.14159265
#define A 100 


int main()
{
    int lt=4000;
    int nx=2000;
    float dx=5.0; 
    int  nz=300;
    float v=1500.0;
    float dz=5.0;
    int ltt=2*lt;
    float dt=0.5;

    int is,it,ix,iz;

    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nx);

    for(it=0;it<lt;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
    for(it=lt;it<ltt;it++)
      omega[it]=-2*pai*1000/(2*dt)+2*pai*(it-lt)*1000/(dt*ltt);

    for(is=0;is<nx/2;is++)
      kx[is]=2*pai*float(is)/float(dx*nx);
   for(is=nx/2;is<nx;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nx/2)/float(dx*nx);


   float **t;
   t=alloc2float(nz,nx);
   
   for(it=1;it<ltt;it++)
    {
     for(ix=0;ix<nx;ix++)
       t[ix][0]=0.0;     

     cout<<"No."<<it<<" frequency time-shift is===="<<endl;
   
     for(ix=0;ix<nx;ix++)
       {
         for(iz=1;iz<nz;iz++)
           {
               if(pow(omega[A]/v,2)-pow(kx[is],2)>0.0)
                   t[ix][iz]=t[ix][iz-1]+dz*sqrt(pow(omega[it]/v,2)-pow(kx[is],2))/omega[it];

               else
                   t[ix][iz]=0.0;
           }
       }

     for(ix=0;ix<nx;ix++)
      {
        cout<<"No."<<ix<<"  wavenumber time-shift is====="<<endl;
        for(iz=0;iz<nz;iz++)
          cout<<iz<<"   "<<t[ix][iz]<<endl;
      }
      return 0;
    } 
     return 0;
    
}
