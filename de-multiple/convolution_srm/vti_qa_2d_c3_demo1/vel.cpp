#define Nx 500
#define Nz 500
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
/*
  for(i=0;i<Nx;i++)
    for(j=0;j<30;j++) 
        v[i][j]=1500.0;

   for(i=0;i<Nx;i++)
//     for(j=60;j<100+i*0.1;j++)
     for(j=30;j<100;j++) 
        v[i][j]=2500.0;

   for(i=0;i<Nx;i++)
//     for(j=100+i*0.1;j<160;j++)
     for(j=100;j<150;j++) 
       v[i][j]=3000.0;

    for(i=0;i<Nx;i++)
     for(j=150;j<200;j++)
        v[i][j]=3500.0;
*/

  for(i=0;i<Nx;i++)
     for(j=0;j<Nz;j++)
        v[i][j]=1500.0;
 
  ofstream swq;
    swq.open("/home/swq/homo_new_sigsbee2a_velocity_500*1000.dat" ,ios::binary);
      
    for(i=0;i<Nx;i++)
       for(j=0;j<Nz;j++)
            swq.write((char*)&v[i][j],sizeof(v[i][j]));
      
    swq.close();
     return 0;


}


