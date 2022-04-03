#define Nx 167
#define Nz 200
using namespace std;
#include <iostream>
#include <fstream>

int main()
{
  int ix,iz;
  float v1;
  v1=1480.0;
  float dz;
  dz=30.0;

  float *dep=new float [Nx];
  for(ix=0;ix<Nx;ix++)
      dep[ix]=0;

  float **v=new float *[Nx];
  for(ix=0;ix<Nx;ix++)
      v[ix]=new float [Nz];

  float *time=new float [Nx];

  for(ix=0;ix<Nx;ix++)
      time[ix]=0.0;
/*
    ifstream swq1;
    swq1.open("/data200/swq/bohai_shallow_data/shallow_water_depth14m.txt");

    for(ix=0;ix<Nx;ix++)
       swq1>>time[ix];
    swq1.close();
    for(ix=0;ix<Nx;ix++)
       dep[ix]=(int)(v1*time[ix]/1000.0/2.0/dz+0.5);
*/
    ofstream swq;
    swq.open("/data200/swq/131_data/shallow_water_depth11m.dat" ,ios::binary);
 
    for(ix=0;ix<Nx;ix++)
       dep[ix]=11;
     
    for(ix=0;ix<Nx;ix++)
       swq.write((char*)&dep[ix],sizeof(dep[ix]));      
    swq.close();
/*
    for(ix=0;ix<Nx;ix++)
     {
       for(iz=0;iz<dep[ix];iz++)
          v[ix][iz]=1480.0;
       for(iz=dep[ix];iz<Nz;iz++)
          v[ix][iz]=3000.0;
     } 

    ofstream swq11;
    swq11.open("/io-raid0/swq/zj_slope_break_data/slope_vel_shot3679to3908_dz30m.dat" ,ios::binary);
      
    for(ix=0;ix<Nx;ix++)
      for(iz=0;iz<Nz;iz++)
       swq11.write((char*)&v[ix][iz],sizeof(v[ix][iz]));
    swq11.close();
*/

     return 0;


}


