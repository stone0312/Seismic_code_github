/* ============================================================================
 File:		  ricker.cpp
 Description: get zero offset data with ricker wave
 Author:      stone
 date:        2021.03.29
 Website:     https://github.com/stone0312/Seismic_code_github.git
 Version:     1.0
 Note:        

 ============================================================================*/
#include<fstream>
#include<iostream>
#include<math.h>

#define pi 3.1415926

using namespace std;

int main(int argc,char *argv[])
{
    int i,j;
    int Xnumber=512;
    int dx=4;
    int Tnumber=1024;
    float dt=0.002;
    int f0=35;
    float **trace;
    
	trace = new float *[Xnumber];

	
	for(i=0;i<Xnumber;i++)
         {
			trace[i] = new float [Tnumber];
		 }
	

	for(i=0;i<Xnumber;i++)
         {

            for(j=0;j<Tnumber;j++)
                {
                     if(i==255)
                     trace[i][j]=(1-2*(pi*f0*(j-150)*dt)*(pi*f0*(j-150)*dt))*exp(-(pi*f0*(j-150)*dt)*(pi*f0*(j-150)*dt));
					 else
                     trace[i][j]=0;

        }
    }

    ofstream ouF;
    ouF.open("./zero_offset_data.dat",ios::binary);
	for (i=0;i<Xnumber;i++)
	{

		for(j=0;j<Tnumber;j++)
		{

		ouF.write((char*)&trace[i][j],sizeof(trace[i][j]));
	
		}

	}
    ouF.close();

	for(i=0;i<Xnumber;i++)
	{
		delete [] trace[i];
	}
	delete [] trace;

	return 0;

}

