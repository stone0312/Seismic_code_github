#define Nx 400
#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "iostream.h"
#include "fstream.h"
#define pai 3.14159265

using namespace std;

int main()
{
  int ix;
/*
  for(ix=200;ix<300;ix++)
     cout<<ix<<"    "<<80*sin(pai/2-pai/3/100*(ix-200))<<endl;
  return 0;
*/
  float *dep=new float [Nx];
  for(ix=0;ix<Nx;ix++)
      dep[ix]=0;
  float *dep1=new float [200];  
/*
    for(ix=0;ix<200;ix++)
       dep[ix]=80*5; 
    for(ix=200;ix<300;ix++)
       dep[ix]=5*80*sin(pai/2-pai/3/100*(ix-200));
    for(ix=300;ix<Nx;ix++)
       dep[ix]=5*40;

    for(ix=0;ix<200;ix++)
       dep1[ix]=dep[ix+150];
*/
    for(ix=0;ix<200;ix++)
       dep1[ix]=5.0*(-0.2*ix+80);

    for(ix=0;ix<200;ix++)
       cout<<ix<<"  "<<dep1[ix]<<endl;

    ofstream swq;
    swq.open("/data2/swq/zj_slope_break_data/2d_mod/water_layer_depth_dz5m_p02.dat" ,ios::binary);
      
    for(ix=0;ix<200;ix++)
       swq.write((char*)&dep1[ix],sizeof(dep1[ix]));
      
    swq.close();
     return 0;


}


