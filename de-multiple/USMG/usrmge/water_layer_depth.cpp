#define Nx 120
#include "math.h"
#include "iostream.h"
#include "fstream.h"

using namespace std;

int main()
{
  int i,j;

  int *dep=new int [Nx];
  for(i=0;i<Nx;i++)
      dep[i]=0;
/*   

  for(i=0;i<Nx;i++)
        dep[i]=5*(60+int(0.1*i+0.5));
*/   

    ofstream swq;
    swq.open("/data2/swq/data_for_mwd_test/water_layer_depth_80.dat" ,ios::binary);
      
    for(i=0;i<Nx;i++)
       swq.write((char*)&dep[i],sizeof(dep[i]));
      
    swq.close();
     return 0;


}


