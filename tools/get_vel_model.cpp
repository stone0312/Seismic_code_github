/* ============================================================================
 File:		  get_vel_model.cpp
 Description: design the vel model 
 Author:      stone
 date:        2021.03.30
 Website:     https://github.com/stone0312/Seismic_code_github.git
 Version:     1.0
 Note:        

 ============================================================================*/
#include<fstream>
#include<iostream>
#include<math.h>

#define nx 497
#define nz 750

using namespace std;

float v[nx][nz] = {0.0};

int main(int argc, char *argv[])
{
	int i,j;
	
	for (i = 0; i < nx ; ++i) 
	{
		for (j = 0; j < nz; ++j) 
		{
			v[i][j]=1500+1.25*j;
		}
	}

	ofstream ouF;

	ouF.open("./vel_model.dat",ios :: binary);

	for (i = 0; i < nx ; ++i) 
	{
		for (j = 0; j < nz; ++j) 
		{
			ouF.write((char*)&v[i][j],sizeof(v[i][j]));
		}
	}

	ouF.close();

	return 0;

}
