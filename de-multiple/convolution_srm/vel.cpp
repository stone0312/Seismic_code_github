#define Nx 200
#define Nz 200
#include "math.h"
#include "iostream.h"
#include "fstream.h"

using namespace std;

int main()
{
  int i,j;

  float **v=new float *[Nx];
  for(i=0;i<Nx;i++)
     v[i]=new float [Nz];

  for(i=0;i<Nx;i++)
    for(j=0;j<5;j++) 
        v[i][j]=1500.0;

   for(i=0;i<Nx;i++)
     for(j=5;j<50;j++)
         v[i][j]=2000.0;

   for(i=0;i<Nx;i++)
     for(50;j<100;j++)
        v[i][j]=2500.0;

    for(i=0;i<Nx;i++)
     for(j=100;j<120;j++)
        v[i][j]=3000.0;
 
    for(i=0;i<Nx;i++)
     for(j=120;j<160;j++)
        v[i][j]=2700.0;

    for(i=0;i<Nx;i++)
     for(j=160;j<200;j++)
        v[i][j]=3500.0;

/*
  for(i=0;i<Nx;i++)
     for(j=0;j<Nz;j++)
        v[i][j]=1500.0;
*/ 
  ofstream swq;
    swq.open("shallow_srm_hor_vel.dat" ,ios::binary);
      
    for(i=0;i<Nx;i++)
       for(j=0;j<Nz;j++)
            swq.write((char*)&v[i][j],sizeof(v[i][j]));
      
    swq.close();
     return 0;


}


