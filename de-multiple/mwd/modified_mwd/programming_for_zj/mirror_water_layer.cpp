#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "iostream.h"
#include "fstream.h"
#define pai 3.14159265
//#define shot 0

using namespace std;

int main()
{
  char fn1[256];
  int Nx;
  float dx;
  float dz; 
  float theta;
  float fai;
  int shot;

  ifstream swq;
  swq.open("mirror_water_layer.par");
  swq>>fn1>>Nx>>dx>>dz;
  swq.close();

  int ix;

  float *dep=new float [Nx];
  for(ix=0;ix<Nx;ix++)
      dep[ix]=0;
   float *dep1=new float [Nx+1];
   float *mir_dep=new float [Nx+1];
   int *mir_receiver=new int [Nx];

  ifstream swq1;
  swq1.open(fn1,ios::binary);
  for(ix=0;ix<Nx;ix++)
    swq1.read((char*)&dep[ix],sizeof(dep[ix]));
  swq1.close();

  for(ix=0;ix<Nx;ix++)
     dep1[ix]=dep[ix];
  dep1[Nx]=dep[Nx-1];

/*
  for(ix=0;ix<Nx;ix++)
    cout<<ix<<"  "<<dep[ix]<<endl;
  return 0;
*/
      
  float k;
  float k1;
  float mir_k;

  int *receiver=new int [Nx];
  
   for(shot=0;shot<Nx;shot++)
   { 
    for(ix=0;ix<Nx;ix++)
       {
/*
         if(dep1[ix+1]==dep1[ix])
           {
              mir_dep[ix]=2*dep1[ix];
              receiver[ix]=2*ix-shot;
              
             cout<<shot<<"  "<<ix<<"   dep===="<<"   "<<dep[ix]<<"    "<<receiver[ix]<<"    "<<mir_dep[ix]<<endl;
           }
*/
//         else 
//           {
             k=(dep1[ix+1]-dep1[ix])/dx;
             theta=atan(k);

             k1=dep1[ix]/((ix-shot)*dx);
             fai=atan(k1);
 
             mir_k=tan(2*theta-fai);

             if(-(float(dep1[ix])-float(mir_k*ix*dx))/mir_k/dx>0.0)
              receiver[ix]=int(-(float(dep1[ix])-float(mir_k*ix*dx))/mir_k/dx+0.5);
             else
              receiver[ix]=int(-(float(dep1[ix])-float(mir_k*ix*dx))/mir_k/dx-0.5);

             mir_dep[ix]=fabs((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx))+sqrt(dep1[ix]*dep1[ix]+(ix-shot)*dx*(ix-shot)*dx))*sin(fai));

//             mir_dep[ix]=2*fabs((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx)))*sin(theta-fai));

          if(cos(fai)>=0)  
            mir_receiver[ix]=shot+int((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx))+sqrt(dep1[ix]*dep1[ix]+(ix-shot)*dx*(ix-shot)*dx))*cos(fai)/dx+0.5);
          else
            mir_receiver[ix]=shot-int((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx))+sqrt(dep1[ix]*dep1[ix]+(ix-shot)*dx*(ix-shot)*dx))*cos(fai)/dx-0.5);

             cout<<shot<<"  "<<ix<<"   dep===="<<"   "<<dep[ix]<<"    "<<receiver[ix]<<"    "<<mir_receiver[ix]<<"    "<<mir_dep[ix]<<"  fai======"<<fai*360/2/pai<<",   theta======"<<theta*360/2/pai<<endl;
//           cout<<shot<<"  "<<ix<<"   "<<"   "<<receiver[ix]<<"    "<<mir_receiver[ix]<<endl;
//          }
       
       }  
/*
   for(ix=0;ix<Nx;ix++)
     cout<<shot<<"  "<<ix<<"            "<<dep[ix]<<"            "<<receiver[ix]<<"          "<<mir_dep[ix]<<endl; 
*/
   }
     return 0;


}


