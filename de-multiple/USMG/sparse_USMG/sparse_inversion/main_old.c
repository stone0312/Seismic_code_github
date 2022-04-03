#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*      C99 standard.   */
#include<complex.h>

#include <sys/time.h>
#include <unistd.h>

#ifndef FILE_NAME_MAX_LENGTH
#define FILE_NAME_MAX_LENGTH 1024
#endif

#ifndef PI
#define PI (3.141592653589793)
#endif

#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

/*      Feng-Bo's lib.  */
//#include "fbsegy.h"
#include "fballoc.h"
#include "fballoc.c"
#include "fbrw.h"
#include "fbrw.c"

#include "prox_sub.c"
//#include "tv_sub.c"
//#include "tv_sub_v2.c"
#include "tv_sub_v2.omp_old.c"

int     readpar(char *parfn, char *vtruefile, char *vobsfile, char *vinvfile, char *vresfile,
		int *nx, int *nz, int *nitermax, float *lambda, float *niu, float *epsmin,
		char *type );
void	normalize_data( int nx, int nz, float **data, float amp );

int	main( int argc , char *argv[] )
{
	if( argc != 2 )
	{
		printf(" Wrong Parameters!\n");
		return 1 ;
	}

	char	parfile[FILE_NAME_MAX_LENGTH] ;
	strcpy(parfile, argv[1]);

	char	vtruefile[FILE_NAME_MAX_LENGTH]="" ;
	char	vobsfile[FILE_NAME_MAX_LENGTH]="" ;
	char	vinvfile[FILE_NAME_MAX_LENGTH]="" ;
	char	vresfile[FILE_NAME_MAX_LENGTH]="" ;
	char	type[FILE_NAME_MAX_LENGTH]="" ;

	int	nx, nz;
	float	lambda, niu;
	int	nitermax;
	float	epsmin;

	// init.
	readpar( parfile, vtruefile, vobsfile, vinvfile, vresfile,
		&nx, &nz, &nitermax, &lambda, &niu, &epsmin, type);

	float	**vtrue	= NULL;
	float	**vobs	= NULL;
	float	**vinv	= NULL;
	float	**vres1	= NULL;
	float	**vres2	= NULL;
	vtrue	= alloc2float(nz, nx);
	vobs	= alloc2float(nz, nx);
	vinv	= alloc2float(nz, nx);
	vres1	= alloc2float(nz, nx);
	vres2	= alloc2float(nz, nx);

	// read data.
	read_2d_float_rb(vtrue, nx, nz, vtruefile);
	read_2d_float_rb(vobs,  nx, nz, vobsfile);

	int	norm_flag = 1;
	float	amp	= 1.0;
	if ( 1 == norm_flag )
	{
	normalize_data(nx, nz, vtrue, amp);
	normalize_data(nx, nz, vobs,  amp);
	}

	// sparse inversion.
//v2.
	l1norm_denoise_tv_reg( vobs, vinv, lambda, niu, type, nx, nz, nitermax, epsmin);
	
	for ( int ix = 0 ; ix < nx ; ++ ix )
	for ( int iz = 0 ; iz < nz ; ++ iz )
	{
		vres1[ix][iz]	= vobs[ix][iz] - vtrue[ix][iz] ;
		vres2[ix][iz]	= vinv[ix][iz] - vtrue[ix][iz] ;
	}


	// output.
	write_2d_float_wb(vinv,  nx, nz, vinvfile);
	write_2d_float_wb(vres2, nx, nz, vresfile);

	// free pointers.
	free2float(vtrue);
	free2float(vobs);
	free2float(vinv);
	free2float(vres1);
	free2float(vres2);

	return 0;
}

int     readpar(char *parfn, char *vtruefile, char *vobsfile, char *vinvfile, char *vresfile,
		int *nx, int *nz, int *nitermax, float *lambda, float *niu, float *epsmin,
		char *type )
{
	FILE	*fp=NULL;
	char	tmp[1024]="";
	if((fp=fopen(parfn,"r"))==NULL)
	{
		printf("Can not open the parameter file: %s to read!\n", parfn);
		return 1;
	}

	fscanf(fp, "%s%s", tmp, vtruefile);
	fscanf(fp, "%s%s", tmp, vobsfile);
	fscanf(fp, "%s%s", tmp, vinvfile);
	fscanf(fp, "%s%s", tmp, vresfile);
	fscanf(fp, "%s%d", tmp, nx);
	fscanf(fp, "%s%d", tmp, nz);
	fscanf(fp, "%s%d", tmp, nitermax);
	fscanf(fp, "%s%f", tmp, lambda);
	fscanf(fp, "%s%f", tmp, niu);
	fscanf(fp, "%s%f", tmp, epsmin);
	fscanf(fp, "%s%s", tmp, type);

	fclose(fp);

	return 0;
}

void	normalize_data( int nx, int nz, float **data, float amp )
{

	float	max = data[0][0];
	for ( int ix = 0 ; ix < nx ; ++ ix )
	for ( int iz = 0 ; iz < nz ; ++ iz )
	{
		if (fabs(data[ix][iz]) > max)
			max = fabs(data[ix][iz]);
	}

	for ( int ix = 0 ; ix < nx ; ++ ix )
	for ( int iz = 0 ; iz < nz ; ++ iz )
		data[ix][iz]	= data[ix][iz] * amp / max;
}
