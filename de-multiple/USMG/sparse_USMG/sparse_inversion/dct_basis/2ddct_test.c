/*======================================================================*
*	TEST of DCT function.						*
*=======================================================================*
*	Author:	Feng Bo.						*
*	Date:	2013.07.11.						*
*=======================================================================*
*======================================================================*/

//Sys.
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<math.h>

#include<complex.h>
#include<fftw3.h>

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
#include "fbsegy.h"
#include "fballoc.h"
#include "fballoc.c"
#include "fbrw.h"
#include "fbrw.c"
#include "prox_sub.c"
#include "tv_sub.c"


void	test_dct2d_denoise( float **in, float **out, int nx, int nz,
		float lambda, float niu, int nitermax, float epsmin );

int	l1_denoise_dct2d( float **Xobs, float **X, float lambda,
		int nx, int nz, int nitermax, float epsmin );

int	l1_denoise_dct2d_test( float **Xobs, float **X, float lambda, float niu,
		int nx, int nz, int nitermax, float epsmin );

int     readpar(char *parfn, char *vtruefile, char *vspecfile, char *vinvfile,
		int *nx, int *nz, int *nitermax, float *lambda, float *niu, float *epsmin );

int main( int argc , char *argv[] )
{
	if( argc != 2 )
	{
		printf(" Wrong Parameters!\n");
		return 1 ;

	}

	char	parfile[FILE_NAME_MAX_LENGTH]="";
	strcpy(parfile, argv[1]);

	char	vtruefile[FILE_NAME_MAX_LENGTH]="";
	char	vspecfile[FILE_NAME_MAX_LENGTH]="";
	char	vinvfile[FILE_NAME_MAX_LENGTH]="";
	int	nx, nz;
	int	nitermax;
	float	lambda, niu, epsmin;

	readpar( parfile, vtruefile, vspecfile, vinvfile,
		&nx, &nz, &nitermax, &lambda, &niu, &epsmin );

	float	**vel	= alloc2float(nz, nx);
	float	**velDCT= alloc2float(nz, nx);

	read_2d_float_rb(vel, nx, nz, vtruefile);

	//test_dct2d_denoise(vel, velDCT, nx, nz, nitermax, lambda, niu, epsmin);
	l1_denoise_dct2d_test(vel, velDCT, lambda, niu, nx, nz, nitermax, epsmin);

	write_2d_float_wb( velDCT, nx, nz, vinvfile);
	
	free2float(vel);
	free2float(velDCT);

	return 0 ;
}

void	test_dct2d_denoise( float **in, float **out, int nx, int nz,
		float lambda, float niu, int nitermax, float epsmin )
{
	printf("begin test_dct2d().\n");
	int	NX	= 2*(nx-1);
	int	NZ	= 2*(nz-1);
	int	NXZ	= NX*NZ;

	fftwf_plan	pf, pb;

	pf	= fftwf_plan_r2r_2d( nx, nz, in[0], out[0], FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE );
	pb	= fftwf_plan_r2r_2d( nx, nz, out[0], in[0], FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE );

	// Forward DCT.
	fftwf_execute(pf);

	write_2d_float_wb( out, nx, nz, "vel.FDCT.dat");

	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
	{
		if ( ix <5 && iz <5 )
		;
		else
		out[ix][iz]	= 0;
	}

/*
	// Denoise.
	float	*obj	= calloc(nitermax, sizeof(float));
	//l1norm_denoise_tv_reg( in, out, lambda, niu, nx, nz, nitermax, epsmin, obj );
	l1norm_inversion_denoise( in[0], out[0], lambda, nx*nz, nitermax, epsmin );
	free(obj);
*/

	// Inverse DCT.
	fftwf_execute(pb);

	// normalize.
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
		in[ix][iz]	/= NXZ;

	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(pb);

	write_2d_float_wb( in, nx, nz, "vel.IDCT.dat");
}

// min{ || X - Xobs || ^ 2 + 2 * lambda * | B * X |1 }
int	l1_denoise_dct2d( float **Xobs, float **X, float lambda,
		int nx, int nz, int nitermax, float epsmin )
{
	int	Mystatus = -1;
	fftwf_plan	pf, pb;

	float	**xtemp	= alloc2float(nz, nx);
	float	**Bxtemp= alloc2float(nz, nx);
	float	**xold	= alloc2float(nz, nx);
	float	**resd	= alloc2float(nz, nx);
	float	**pk	= alloc2float(nz, nx);
	float	**BTpk	= alloc2float(nz, nx);

	memset((void*)&xtemp[0][0],  0, sizeof(float)*nx*nz);
	memset((void*)&xold[0][0],   0, sizeof(float)*nx*nz);
	memset((void*)&resd[0][0],   0, sizeof(float)*nx*nz);
	memset((void*)&pk[0][0],     0, sizeof(float)*nx*nz);
	memset((void*)&BTpk[0][0],   0, sizeof(float)*nx*nz);

	// Initialize variables.

	int	nxz	= nx * nz;
	float	beta	= 1.0/(lambda);

	pf	= fftwf_plan_r2r_2d( nx, nz, xtemp[0], Bxtemp[0], FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE );
	pb	= fftwf_plan_r2r_2d( nx, nz, pk[0],    BTpk[0],   FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE );

	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			xold[ix][iz] = xtemp[ix][iz];

		// calculate BT*p[k] using Inverse DCT, from p[k] to BTpk.
		fftwf_execute(pb);

		// denoise.
		// x[k+1] = Xobs - lambda * BT * p[k];
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			xtemp[ix][iz]	= Xobs[ix][iz] - lambda * BTpk[ix][iz]/nxz ;

		// calculate B*x[k+1] using Forward DCT, from x[k+1] to B*x[k+1].
		fftwf_execute(pf);

		// ptemp = p[k] + beta * B * x[k+1];
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			pk[ix][iz] = pk[ix][iz] + beta * Bxtemp[ix][iz];

		// calculate prox.
		sprox_l1norm_vector_dual( pk[0], pk[0], nx*nz, beta );

		// 4. calculate the residuals.
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			resd[ix][iz] = xtemp[ix][iz] - xold[ix][iz];

		double	sum1	= inner_product_vector_1d_f( resd[0],  resd[0],  nx*nz ) ;
		double	sum2	= inner_product_vector_1d_f( xtemp[0], xtemp[0], nx*nz ) ;
		double	re	= sum1/sum2;

		printf("res[%d]=%8.6e\n", iter, re);

		if ( re < epsmin )
		{
			Mystatus = 0;
			break;
		}
	}

	// output.
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
		X[ix][iz] = xtemp[ix][iz];
	
	free2float(xtemp);
	free2float(Bxtemp);
	free2float(xold);
	free2float(resd);
	free2float(pk);
	free2float(BTpk);

	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(pb);

	write_2d_float_wb( X, nx, nz, "vel.IDCT.dat");

	return	Mystatus;
}


int     readpar(char *parfn, char *vtruefile, char *vspecfile, char *vinvfile,
		int *nx, int *nz, int *nitermax, float *lambda, float *niu, float *epsmin )
{
	FILE	*fp=NULL;
	char	tmp[1024]="";
	if((fp=fopen(parfn,"r"))==NULL)
	{
		printf("Can not open the parameter file: %s to read!\n", parfn);
		return 1;
	}

	fscanf(fp, "%s%s", tmp, vtruefile);
	fscanf(fp, "%s%s", tmp, vspecfile);
	fscanf(fp, "%s%s", tmp, vinvfile);
	fscanf(fp, "%s%d", tmp, nx);
	fscanf(fp, "%s%d", tmp, nz);
	fscanf(fp, "%s%d", tmp, nitermax);
	fscanf(fp, "%s%f", tmp, lambda);
	fscanf(fp, "%s%f", tmp, niu);
	fscanf(fp, "%s%f", tmp, epsmin);

	fclose(fp);

	return 0;
}


/*
float   ricker_t0(float t, float fm, float t0)
{
        //  Ricker FB.
        float x=pow(PI*fm*(t-t0),2);
        return exp(-x)*(1-2*x);
}


void	rickerfunc( float *data, int ns, float dt, float fm, float t0 )
{
	for ( int it = 0 ; it < ns ; ++ it )
		data[it]	= ricker_t0(dt*(it+1), fm, t0);
}


int main( int argc , char *argv[] )
{
	int	ns;
	float	dt, fm, t0;
	float	*trace;
	float	*trace_dct;
	float	*trace_dft;

	ns	= 1000;
	dt	= 0.001;
	fm	= 10.0;
	t0	= 0.3;

	trace		= calloc(ns, sizeof(float));
	trace_dct	= calloc(ns, sizeof(float));
	trace_dft	= calloc(ns, sizeof(float));

	rickerfunc(trace, ns, dt, fm, t0);

	write_1d_float_wb(trace, ns, "trace.dat");

	test_dct1d( trace, trace_dct, ns );
	//test_dft1d( trace, trace_dft, ns );

	free(trace);
	free(trace_dct);
	free(trace_dft);

	return 0 ;
}

*/
// min{ || X - Xobs || ^ 2 + 2 * lambda * | B * X |1 }
int	l1_denoise_dct2d_test( float **Xobs, float **X, float lambda, float niu,
		int nx, int nz, int nitermax, float epsmin )
{
	int	Mystatus = -1;
	fftwf_plan	pf, pb;

	float	**xtemp	= alloc2float(nz, nx);
	float	**Bxtemp= alloc2float(nz, nx);
	float	**xold	= alloc2float(nz, nx);
	float	**resd	= alloc2float(nz, nx);
	float	**pk	= alloc2float(nz, nx);
	float	**BTpk	= alloc2float(nz, nx);

	memset((void*)&xtemp[0][0],  0, sizeof(float)*nx*nz);
	memset((void*)&xold[0][0],   0, sizeof(float)*nx*nz);
	memset((void*)&resd[0][0],   0, sizeof(float)*nx*nz);
	memset((void*)&pk[0][0],     0, sizeof(float)*nx*nz);
	memset((void*)&BTpk[0][0],   0, sizeof(float)*nx*nz);

	// Initialize variables.

	int	nxz	= nx * nz;
	float	beta	= 1.0/(lambda);
	niu	*= nxz;

	pf	= fftwf_plan_r2r_2d( nx, nz, xtemp[0], Bxtemp[0], FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE );
	pb	= fftwf_plan_r2r_2d( nx, nz, pk[0],    BTpk[0],   FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE );

	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			xold[ix][iz] = xtemp[ix][iz];

		// calculate BT*p[k] using Inverse DCT, from p[k] to BTpk.
		fftwf_execute(pb);

		// denoise.
		// x[k+1] = Xobs - lambda * BT * p[k];
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			xtemp[ix][iz]	= Xobs[ix][iz] - lambda * BTpk[ix][iz]/nxz ;

		// calculate B*x[k+1] using Forward DCT, from x[k+1] to B*x[k+1].
		fftwf_execute(pf);

		// ptemp = p[k] + beta * B * x[k+1];
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			pk[ix][iz] = pk[ix][iz] + 1.0 * Bxtemp[ix][iz];

		// calculate prox.
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
		{
			float	tmp1	= pk[ix][iz];
			float	tmp2;
			tmp2	= MAX( fabs(tmp1)-niu/lambda, 0 ) * SGN(tmp1);
			pk[ix][iz]	-= tmp2;
		}

/*
		// calculate prox.
		sprox_l1norm_vector_dual( pk[0], pk[0], nx*nz, beta );
*/

		// 4. calculate the residuals.
		for( int ix = 0 ; ix < nx ; ++ ix )
		for( int iz = 0 ; iz < nz ; ++ iz )
			resd[ix][iz] = xtemp[ix][iz] - xold[ix][iz];

		double	sum1	= inner_product_vector_1d_f( resd[0],  resd[0],  nx*nz ) ;
		double	sum2	= inner_product_vector_1d_f( xtemp[0], xtemp[0], nx*nz ) ;
		double	re	= sum1/sum2;

		printf("res[%d]=%8.6e\n", iter, re);

		if ( re < epsmin )
		{
			Mystatus = 0;
			break;
		}
	}

	// output.
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
		X[ix][iz] = xtemp[ix][iz];
	
	free2float(xtemp);
	free2float(Bxtemp);
	free2float(xold);
	free2float(resd);
	free2float(pk);
	free2float(BTpk);

	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(pb);

	//write_2d_float_wb( X, nx, nz, "vel.IDCT.dat");

	return	Mystatus;
}
