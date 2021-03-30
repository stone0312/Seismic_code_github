/* ============================================================================
 File:		  vel_model_change.cpp
 Description: modify the vel model 
 Author:      stone
 date:        2021.03.30
 Website:     https://github.com/stone0312/Seismic_code_github.git
 Version:     1.0
 Note:        

 ============================================================================*/
#include<fstream>
#include<iostream>
#include<math.h>

//#define nx 750
//#define nz 497

using namespace std;


//float v[nx][nz]={0.0};


int main(int argc, char *argv[])
{
	int i,j;
	int nx=750;
	int nz=497;

	float **v;

	v= new float *[nx];
	
		for (i = 0; i < nx; ++i) 
		{
			v[i] = new float [nz];
		}

	ifstream inF;

	inF.open("./mar_v_dx12.5m_dz4m_750x497.dat",ios :: binary);

	for (i = 0; i < nx ;++i) 
	{
		for (j = 0; j < nz; ++j) 
		{
			inF.read((char*)&v[i][j],sizeof(v[i][j]));
		}
	}


	for (i = 0; i < nx ;++i) 
	{
		for (j = 0; j < nz; ++j) 
		{
			v[i][j]=v[i][j]/2;
		}
	}


	ofstream ouF;

	ouF.open("./vel_change.dat",ios :: binary);

	for (i = 0; i < nx ;++i) 
	{
		for (j = 0; j < nz; ++j) 
		{
			ouF.write((char*)&v[i][j],sizeof(v[i][j]));
		}
	}
	
	for (i = 0; i < nx; ++i) 
	{
		delete  [] v[i];
	}
	delete [] v;

	inF.close();
	ouF.close();

	return 0;

}
