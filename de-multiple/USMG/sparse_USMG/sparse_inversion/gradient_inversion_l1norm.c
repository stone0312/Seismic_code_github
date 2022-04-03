#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>
#include <unistd.h>

// C99 standard.
#include<complex.h>

#define	FILE_NAME_MAX_LENGTH	1000

#define	USE_OPENMP

// OpenMP
#ifndef USE_OPENMP
	#define omp_get_thread_num() 0
#else
	#include<omp.h>
#endif

// MPI
//#include "mpi.h"

#include"FB_DEFINE.h"

#include"fballoc.h"
#include"fballoc.c"
#include"fbrw.h"
#include"fbrw.c"

#include"prox_sub.c"
#include"tv_sub.c"

/*================================================================================================*/
int	calculate_grid_num( int *Ix1, int *Ix2, int *Iz1, int *Iz2,
		int *Nx1, int *Nz1,
		float x0, float z0, float dx, float dz,
		float x1, float x2, float z1, float z2 );

int	load_tsk( float **tsk, int nx, int nz, int npair, int ngrid1,
		float x0, float z0, float dx, float dz,
		float x1, float x2, float z1, float z2,
		char *fn );

int	precondition_gradient_solver( float **k, float *b, float *x,
		int nx, int nz, int npair, char *fnGrad, char *fnHessian );

int	readpar( char *parfn, int *nx, int *nz, int *npair, float *dx, float *dz,
		int *nitermax, float *lambda, float *niu, float *epsmin,
		char *fntres, char *fnwavepath, char *fnHessian, char *fnGrad,  char *fnVelInv );

int	main( int argc , char *argv[] )
{
	char	parfile[FILE_NAME_MAX_LENGTH] ;
	strcpy(parfile,    argv[1]);

	int	nx, nz, ngrid, npair;
	float	dx, dz;
        float   epsmin  = 1E-5;         // the default relative travel-time error.
        int     nitermax= 200;          // the default maximum iteration number.
	float	lambda, niu;		// the damp factor.

        char    fntres[FILE_NAME_MAX_LENGTH]="";
        char    fnwavepath[FILE_NAME_MAX_LENGTH]="";
        char    fnHessian[FILE_NAME_MAX_LENGTH]="";
        char    fnGrad[FILE_NAME_MAX_LENGTH]="";
        char    fnVelInv[FILE_NAME_MAX_LENGTH]="";

	readpar( parfile, &nx, &nz, &npair, &dx, &dz,
		&nitermax, &lambda, &niu, &epsmin,
		fntres, fnwavepath,
		fnHessian, fnGrad, fnVelInv );

	int	nx1, nz1, ngrid1;
	float	x0, z0, x1, x2, z1, z2;
	x0	= 0;
	z0	= 0;
	x1	= 1000;
	x2	= 5000;
	z1	= 0;
	z2	= 2000;
	dx	= 40;
	dz	= 40;

	int	tmp;
	ngrid1	= calculate_grid_num( &tmp, &tmp, &tmp, &tmp,  &nx1, &nz1,
		x0, z0, dx, dz, x1, x2, z1, z2 );

	float	**tsk	= alloc2float(ngrid1, npair);
	load_tsk(tsk, nx, nz, npair, ngrid1,
		x0, z0, dx, dz, x1, x2, z1, z2, fnwavepath );

	float	**vpert	= alloc2float(nz1, nx1);
	float	**vbg	= alloc2float(nz1, nx1);
	float	**vall	= alloc2float(nz1, nx1);
	float	**vinv	= alloc2float(nz1, nx1);
	float	*tres	= calloc(npair, sizeof(float));

	printf("ok2.\n");

       	read_1d_float_rb(tres, npair, fntres);

	precondition_gradient_solver(tsk, tres, vpert[0], nx1, nz1, npair, fnGrad, fnHessian);

	write_2d_float_wb(vpert, nx1, nz1, fnVelInv);

	float	steplength;

	steplength	= 1E-3;
	for ( int ix = 0 ; ix < nx1 ; ++ ix )
	for ( int iz = 0 ; iz < nz1 ; ++ iz )
	{
		vbg[ix][iz]	= 4000;
		vall[ix][iz]	= vbg[ix][iz] + steplength*vpert[ix][iz] ;
	}

	//denoise.
	float	*obj	= calloc(nitermax, sizeof(float));
	l1norm_denoise_tv_reg( vall, vinv, lambda, niu, nx1, nz1, nitermax, epsmin, obj);

	write_2d_float_wb(vall, nx1, nz1, "vel_raw.dat");
	write_2d_float_wb(vinv, nx1, nz1, "vel_denoise.dat");

	// free.
	free2float(tsk);
	free2float(vpert);
	free2float(vbg);
	free2float(vall);
	free2float(vinv);
	free(tres);

	return 0 ;
}

int	readpar( char *parfn, int *nx, int *nz, int *npair, float *dx, float *dz,
		int *nitermax, float *lambda, float *niu, float *epsmin,
		char *fntres, char *fnwavepath, char *fnHessian, char *fnGrad,  char *fnVelInv )
{
	FILE	*fp=NULL;
	char	tmp[1024]="";
	if((fp=fopen(parfn,"r"))==NULL)
	{
		printf("Can not open paramter file: %s to read!\n", parfn);
		return 1;
	}

	// Line 01-03.
	fscanf(fp, "%s%d", tmp, nx);
	fscanf(fp, "%s%d", tmp, nz);
	fscanf(fp, "%s%d", tmp, npair);
	fscanf(fp, "%s%f", tmp, dx);
	fscanf(fp, "%s%f", tmp, dz);

	fscanf(fp, "%s%d", tmp, nitermax);
	fscanf(fp, "%s%f", tmp, lambda);
	fscanf(fp, "%s%f", tmp, niu);
	fscanf(fp, "%s%f", tmp, epsmin);

	fscanf(fp, "%s%s", tmp, fntres);
	fscanf(fp, "%s%s", tmp, fnwavepath);
	fscanf(fp, "%s%s", tmp, fnHessian);
	fscanf(fp, "%s%s", tmp, fnGrad);
	fscanf(fp, "%s%s", tmp, fnVelInv);

	fclose(fp);

	return 0;
}

//	solve A * x = b, using the diag-hessian as precondition.
int	precondition_gradient_solver( float **k, float *b, float *x,
		int nx, int nz, int npair, char *fnGrad, char *fnHessian )
{
	int	ngrid	= nx * nz ;
	float	*grad	= calloc(ngrid, sizeof(float));
	float	*dhess	= calloc(ngrid, sizeof(float));

	// calculate the gradient: AT * b and the diagonal hessian.
#pragma omp parallel for
        for( int igrid = 0 ; igrid < ngrid; ++ igrid )
	{
	        for( int iray  = 0 ; iray  < npair ; ++ iray  )
		{
			float	temp	= k[iray][igrid];
			grad[igrid]	+= temp * b[iray];
			dhess[igrid]	+= temp * temp ;
		}
	}

	float	max	= -1;
	float	eps	= 1E-2;
        for( int igrid = 0 ; igrid < ngrid; ++ igrid )
	{
		if( dhess[igrid] > max )
			max	= dhess[igrid];
	}
	float	thres	= max *eps;

        for( int igrid = 0 ; igrid < ngrid; ++ igrid )
	{
		if ( dhess[igrid] > thres )
			x[igrid]	= grad[igrid]/dhess[igrid];
		else
			x[igrid]	= 0;
	}

	write_1d_float_wb(grad,  ngrid, fnGrad);
	write_1d_float_wb(dhess, ngrid, fnHessian);

	return 0;
}

int	load_tsk( float **tsk, int nx, int nz, int npair, int ngrid1,
		float x0, float z0, float dx, float dz,
		float x1, float x2, float z1, float z2,
		char *fn )
{
	float	***tsk1	= alloc3float(nz, nx, npair);
       	read_3d_float_rb(tsk1, npair, nx, nz, fn);

	int	ix1, ix2, iz1, iz2;
	int	nx1, nz1;

	calculate_grid_num(&ix1, &ix2, &iz1, &iz2, &nx1, &nz1,
		x0, z0, dx, dz, x1, x2, z1, z2);

	for ( int ir = 0 ; ir < npair ; ++ ir )
	{
		int	igrid	= 0;
		for ( int ix = ix1 ; ix <= ix2 ; ++ ix )
		for ( int iz = iz1 ; iz <= iz2 ; ++ iz )
		{
			tsk[ir][igrid] = tsk1[ir][ix][iz];
			++ igrid;
		}
		if( ngrid1 != igrid )
			printf("Error!\n igrid=%d  ngrid1=%d \n", igrid, ngrid1);
	}

	free3float(tsk1);

	return 0;
}

int	calculate_grid_num( int *Ix1, int *Ix2, int *Iz1, int *Iz2,
		int *Nx1, int *Nz1,
		float x0, float z0, float dx, float dz,
		float x1, float x2, float z1, float z2 )
{
	int	ix1, ix2, iz1, iz2;
	int	nx1, nz1;

	ix1	= (int)(x1-x0+0.5)/dx;
	ix2	= (int)(x2-x0+0.5)/dx;
	iz1	= (int)(z1-z0+0.5)/dz;
	iz2	= (int)(z2-z0+0.5)/dz;
	nx1	= ix2 - ix1 + 1;
	nz1	= iz2 - iz1 + 1;

	*Ix1	= ix1;
	*Ix2	= ix2;
	*Iz1	= iz1;
	*Iz2	= iz2;
	*Nx1	= nx1;
	*Nz1	= nz1;

	printf(" ix1 = %d ix2 = %d nx1 = %d \n", ix1, ix2, nx1);
	printf(" iz1 = %d iz2 = %d nz1 = %d \n", iz1, iz2, nz1);

	return	nx1*nz1;
}
