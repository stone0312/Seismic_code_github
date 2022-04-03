/* ============================================================================
 File:		  convloution_by_armadillo.cpp
 Description: use armadillo data to get forward result by convloution
 Author:      stone
 date:        2021.08.26
 Website:     https://github.com/stone0312/Seismic_code_github.git
 Version:     1.0
 Note:        

 ============================================================================*/
#include<fstream>
#include<iostream>
#include<math.h>
#include<armadillo>

#define PI 3.1415926

using namespace std;
using namespace arma;


int main(int argc, char *argv[])
{
	
	// design r
	int i,j;
	int t = 251;
	int x = 81;
	
	mat r;
	r.zeros(t,x);
	for (j=0; j<x; j++)
	{
		r(50,j)=1;
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
        ricker[i+number/2]=(1-2*PI*PI*freq*freq*nt*nt)
                    *exp(-PI*PI*freq*freq*nt*nt);
	}	
	
	vec ricker2;
	ricker2.zeros(30);
	for(i=0;i<30;i++)
	{
		ricker2(i)=ricker[i];
	}
	

	// Do convolution 
	
	//mat b = conv2(r,ricker2);
	mat b;
	b.zeros(t+number-1,x);
	for (j=0; j<x; j++)
	{
		b.col(j)=conv(r.col(j),ricker2);
	}
	
	mat c;
	c.zeros(t,x);
	c = b.submat(span(number/2,t+number/2-1),span(0,x-1));

	// output the result
	float **result;

	result = new float *[t];

	for (i=0;i<t;i++)
	{
		result[i]= new float [x];
	}

	for (i=0;i<t;i++)
	{
		for(j=0; j<x; j++)
		{	
			result[i][j]=c(i,j);
		}
	}
	
	ofstream ouF;
    ouF.open("./forward_result_by_conv_n1=251.dat",ios::binary);
	for (j=0;j<x;j++)
	{

		for(i=0;i<t;i++)
		{

		ouF.write((char*)&result[i][j],sizeof(result[i][j]));
	
		}

	}
    ouF.close();


	delete [] ricker;

	for (int i=0; i<x; i++) 
		{delete[] result[i];}
	return 0;
}



