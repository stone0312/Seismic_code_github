/* ============================================================================
 File:		  convloution.cpp
 Description: get forward result by convloution 
 Author:      stone
 date:        2021.08.19
 Website:     https://github.com/stone0312/Seismic_code_github.git
 Version:     1.0
 Note:        

 ============================================================================*/
#include<fstream>
#include<iostream>
#include<math.h>


#define nx 301
#define nz 501
#define PI 3.1415926
using namespace std;


//void conv(float x[],int m,float h[],int n,float y[],int l);

int main(int argc, char *argv[])
{
	
	// design r
	int i,j;
	int t = 100;
	int x = 10;
	float **r;

	r = new float *[10];

	for (i=0;i<x;i++)
	{
		r[i]= new float [t];
	}

	for (i=0;i<x;i++)
	{
		for(j=0;j<t;j++)
		{if(i==5)
				r[i][j]=1;
			else r[i][j]=0;
		}
	}
	
	// design ricker
	
	float freq=35;
    float dt=0.004;
	int number=30;
	float *ricker;
	float nt;

	ricker = new float [30];
    
	for(i=-number/2;i<=number/2;i++)
    {
		nt=dt*i;    
        ricker[i]=(1-2*PI*PI*freq*freq*nt*nt)
                    *exp(-PI*PI*freq*freq*nt*nt);
	}	
	
	float **result;

	result = new float *[10];

	for (i=0;i<x;i++)
	{
		result[i]= new float [t+number-1];
	}

	for (i=0;i<x;i++)
	{
		for(j=0;j<t+number-1;j++)
				result[i][j]=0;
	}

	cout << "this ok 1" << endl;
    
	// Do convolution r(n)*ricker(n)
//	conv(ricker,number,r,t,result,number+t-1);

	// output the result
	ofstream ouF;
    ouF.open("./data.dat",ios::binary);
	for (i=0;i<x;i++)
	{

		for(j=0;j<t+number-1;j++)
		{

		ouF.write((char*)&result[i][j],sizeof(result[i][j]));
	
		}

	}
    ouF.close();

	delete [] ricker;
	for (int i=0; i<x; i++) 
		{delete[] result[i];}
	delete [] result;
	for (int i=0; i<x; i++) 
		{delete[] r[i];}
	delete [] r;
	return 0;
}

// The function of convolution.

void conv(float x[],int m,float h[],int n,float y[],int l)
{ 
	int k,i;
	
	for(k=0;k<l;k++)
    { 
		y[k]=0.0;
		  for(i=0;i<m;i++)
			{ 
			if(k-i>=0&&k-i<=n-1)
			 y[k]=y[k]+x[i]*h[k-i];
			}
    }
}


