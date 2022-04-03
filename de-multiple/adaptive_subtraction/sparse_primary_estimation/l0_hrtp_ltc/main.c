/*C's stadard*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

/* c99 stadard */
#include<complex.h>

/*   DEFINE SOMETHING   */
#ifndef FILE_NAME_LENGTH
#define FILE_NAME_LENGTH 256
#endif

#ifndef PAI
#define PAI (3.141592653589793)
#endif

#ifndef MAX
#define MAX(a,b) (a)>(b)?(a):(b)
#endif

#ifndef MIN
#define MIN(a,b) (a)>(b)?(b):(a)
#endif

/*  Feng-Bo's lib.  */
#include "alloc.c"

/* HR_RT's subroutine */
#include "./subroutine/norm.c"
#include "./subroutine/minn.c"
#include "./subroutine/maxn.c"
#include "./subroutine/sis_radon_liutc_v_beam_psdm.c"

int main()
{
	int i;
	int nx, np, ns;
  float *pcycle;
  float pmax, dx, dp;
  float *p, *cor_in;
  float **data, **taup;

/* The initial value */ 
  nx= 99;
  np= 61;
  ns= 1000;
  pmax= sin(80.0*PAI/180.0);
  dx= 15.0;
  dp= 2*pmax/(np-1);

  printf("dp=%f, p=%f\n", dp, pmax);

/* allocate the memory */
	p= alloc1float(np);
	pcycle= alloc1float(np);
	cor_in= alloc1float(nx);
	data= alloc2float(ns, nx);
	taup= alloc2float(ns, np);
  zero1float(p, np);
  zero1float(pcycle, np);
  zero1float(cor_in, nx);
  zero2float(data, ns, nx);
  zero2float(taup, ns, np);

/* Read the data */
  FILE *fp_read= NULL;
  fp_read= fopen("./data_20x.dat", "rb");
	if(fp_read==NULL)
  printf("the file is failed to read");
	fread(&data[0][0], sizeof(float), nx*ns, fp_read);
  fclose(fp_read);

/* Calculation of x coordation */
  for(i=0; i< nx; i++)
	{
		cor_in[i]=i* dx;
	}

/* Calculation of scanning P */
  for(i=(int)(-(np-1)/2); i<= (int)(np-1)/2; i++)
	{
	  p[i+(np-1)/2]=i* dp;
	}

/* Do the high resolution Radon Transform */
  sis_RT_HR_Time_2D(ns, nx, 0.001, cor_in, data, np, p, taup, pcycle);

/* Output the Radon Transform result */
  FILE *fp_write_tp= NULL;
  fp_write_tp= fopen("./taup.dat", "wb");
	if(fp_write_tp==NULL)
  printf("the file is failed to write");
	fwrite(&taup[0][0], sizeof(float), ns*np, fp_write_tp);
  fclose(fp_write_tp);

  FILE *fp_write= NULL;
  fp_write= fopen("./pcycle.dat", "wb");
	if(fp_write==NULL)
  printf("the file is failed to write");
	fwrite(&pcycle[0], sizeof(int), np, fp_write);
  fclose(fp_write);

/* Release the Memory */
	free1float(pcycle);
	free1float(p);
	free1float(cor_in);
	free2float(data);
	free2float(taup);

	return 1;
}
