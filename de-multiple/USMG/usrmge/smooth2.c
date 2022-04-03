#include "iostream.h"
#include "fstream.h"
#include "math.h"
#include "complex.h"
#include "alloc.c"
#define pai 3.14159265

void tripd(float *d, float *e, float *b, int n);
int smooth2 (float *vel, float r1,  float *tripd_d, float *tripd_e, float *tripd_f, int n1);
int smooth2 (float *vel, float r1,  float *tripd_d, float *tripd_e, float *tripd_f, int n1);
int     smooth2fb( float *vel, float r1, int nx);

int     smooth2fb( float *vel, float r1, int nx)
{

        /* scale the smoothing parameter */
        r1 = r1*r1*0.25;

        float   *dd;
        dd=alloc1float(nx);
        float   *ee;
        ee=alloc1float(nx);
        float   *ff;
        ff=alloc1float(nx);

        int     status  = 0 ;
        status  = smooth2(vel, r1, dd, ee, ff,nx);

        free(dd);
        free(ee);
        free(ff);

        return status;
}

void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Input:
d[]	the diagonal of A 
e[]	the superdiagonal of A
b[]	the rhs of Ax=b

Output:
b[]	b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
*****************************************************************************/
{
	int k; 
	float temp;
	
	/* decomposition */
	for(k=1; k<n; ++k){
           temp = e[k-1];
           e[k-1] = temp/d[k-1];
           d[k] -= temp*e[k-1];
	}

	/* substitution	*/
        for(k=1; k<n; ++k)  b[k] -= e[k-1]*b[k-1];
	
        b[n-1] /=d[n-1];
        for(k=n-1; k>0; --k)  b[k-1] = b[k-1]/d[k-1] - e[k-1]*b[k]; 
	
 }

 int smooth2 (float *vel, float r1, float *tripd_d, float *tripd_e, float *tripd_f, int n1)
 {
 	int ix;

        for(ix=0; ix<n1-1; ++ix)
	  {
	     tripd_d[ix] = 1.0+r1*2.0;
	     tripd_e[ix] = -r1;
	     tripd_f[ix] = vel[ix];
	  }
	tripd_d[0] -= r1;
	tripd_d[n1-1] = 1.0+r1;
	tripd_f[n1-1] = vel[n1-1];
	tripd(tripd_d,tripd_e,tripd_f,n1);
	for(ix=0; ix<n1; ++ix)
	   vel[ix] = tripd_f[ix];

	return 0;
 }

int main()
{
   int lt,ltt;
   lt=4000;
   ltt=lt+400;

   float a[ltt];

   float b[lt];
   
   int it;
  
   ifstream swq1;
   swq1.open("deghost_amp_1st_trace_100m_lamda100.dat",ios::binary);

   ofstream swq2;
   swq2.open("smooth20_deghost_amp_1st_trace_100m_lamda100.dat",ios::binary);
 
   for(it=0;it<ltt;it++)
      a[it]=0.0;

   for(it=200;it<lt+200;it++) 
      swq1.read((char*)&a[it],sizeof(a[it]));

   smooth2fb( a, 20.0, ltt);

   for(it=0;it<lt/2+1;it++)
      b[it]=a[it+150]; 
   for(it=lt/2+1;it<lt;it++)
      b[it]=b[lt-it]; 

   for(it=0;it<lt;it++)    
      swq2.write((char*)&b[it],sizeof(b[it]));

   swq1.close();
   swq2.close();

   return 0;
}

































































