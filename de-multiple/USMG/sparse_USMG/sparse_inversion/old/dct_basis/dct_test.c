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

#define nn 8
#define NN ( 2*((nn)-1) )

int     readpar(char *parfn, char *vtruefile, char *vspecfile, char *vinvfile,
		int *nx, int *nz, int *nitermax, float *lambda, float *niu, float *epsmin );

void	test_dct2d( float **in, float **out, int nx, int nz )
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
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
		out[ix][iz]	= fabs(out[ix][iz]);
	write_2d_float_wb( out, nx, nz, "fzx.DCT.dat");
	return;

	float	percx	= 0.03;
	float	percz	= 0.1;
	int	nx_truncate	= (int)(nx*percx);
	int	nz_truncate	= (int)(nz*percz);
	nx_truncate	= 30;
	nz_truncate	= 15;
	printf("truncate the small coefficient.\t");
	printf("nx_truncate=%d nz_truncate=%d\n", nx_truncate, nz_truncate);

	double	total_energy=0;
	double	ratio_energy=0.98;
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
		total_energy	+= out[ix][iz]*out[ix][iz];

	int	x_idx=0, z_idx=0;
	for( int iz = 1 ; iz < nz ; ++ iz )
	{
		int	ix	= iz*nx/nz;
		double	part_energy=0;
		for( int jx = 0 ; jx < ix ; ++ jx )
		for( int jz = 0 ; jz < iz ; ++ jz )
			part_energy	+= out[jx][jz]*out[jx][jz];

		if ( part_energy >= ratio_energy*total_energy)
		{
			x_idx	= ix;
			z_idx	= iz;
			break;
		}
	}
	printf("x_idx=%d z_idx=%d\n", x_idx, z_idx);
	nx_truncate	= x_idx;
	nz_truncate	= z_idx;

	// truncate the small coefficient.
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
	{

		if( ix <= nx_truncate && iz <= nz_truncate )
			out[ix][iz]	*= 1.0;
		else
			out[ix][iz]	= 0;
	}

	// Inverse DCT.
	fftwf_execute(pb);

	// normalize.
	for( int ix = 0 ; ix < nx ; ++ ix )
	for( int iz = 0 ; iz < nz ; ++ iz )
		in[ix][iz]	/= NXZ;

	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(pb);

	write_2d_float_wb( in, nx, nz, "fzx.back.dat");
}

void	test_dct2d_fastdim_DCT( float **in, float **out, int nx, int nz )
{
	int	NX	= 2*(nx-1);
	int	NZ	= 2*(nz-1);
	int	NXZ	= NX*NZ;

	int	nz_truncate	= 20;
	float	*buf_in	= calloc(nz, sizeof(float));
	float	*buf_out= calloc(nz, sizeof(float));

	fftwf_plan	pf, pb;
	pf	= fftwf_plan_r2r_1d( nz, buf_in, buf_out, FFTW_REDFT00, FFTW_ESTIMATE );
	pb	= fftwf_plan_r2r_1d( nz, buf_out, buf_in, FFTW_REDFT00, FFTW_ESTIMATE );

	for( int ix = 0 ; ix < nx ; ++ ix )
	{
		for( int iz = 0 ; iz < nz ; ++ iz )
			buf_in[iz]	= in[ix][iz];

		fftwf_execute(pf);

		for( int iz = 0 ; iz < nz ; ++ iz )
			out[ix][iz]	= buf_out[iz] ;

		// truncate the small coefficient.
		for( int iz = nz_truncate ; iz < nz ; ++ iz )
			buf_out[iz] = 0 ;

		fftwf_execute(pb);

		for( int iz = 0 ; iz < nz ; ++ iz )
			in[ix][iz]	= buf_in[iz]/NZ ;

	}
	write_2d_float_wb( out, nx, nz, "fzx.DCT.dat");
	write_2d_float_wb( in,  nx, nz, "fzx.back.dat");

	fftwf_destroy_plan(pf);
	fftwf_destroy_plan(pb);

	free(buf_in);
	free(buf_out);

}

void	test_dct1d( float *in, float *out, int n )
{
	//write_1d_float_wb( in, n, "fx.dat");

	fftwf_plan	p;

	p	= fftwf_plan_r2r_1d( n, in, out, FFTW_REDFT11, FFTW_ESTIMATE );

	fftwf_execute(p);

	for( int i = 0 ; i < n ; ++ i )
		out[i]	= fabs(out[i]);
	write_1d_float_wb( out, n, "fw.dct.dat");

/*
	for( int i = 0 ; i < n ; ++ i )
		printf("in[%d]=%f out[i]=%f\n", i, in[i], out[i]/sqrt(NN));

*/

	p	= fftwf_plan_r2r_1d( n, out, in, FFTW_REDFT11, FFTW_ESTIMATE );

	fftwf_execute(p);

/*
	for( int i = 0 ; i < n ; ++ i )
		printf("in[%d]=%f out[i]=%f\n", i, in[i]/NN, out[i]);
*/

	fftwf_destroy_plan(p);

	write_1d_float_wb( in, n, "ft.idct.dat");
}

void	test_dft1d( float *in1, float *out1, int n )
{
	//write_1d_float_wb( in, n, "fx.dat");

	float complex	*in	= calloc(n, sizeof(float complex));
	float complex	*out	= calloc(n, sizeof(float complex));
	fftwf_plan	p;

	for( int i = 0 ; i < n ; ++ i )
		in[i]	= in1[i] + 0*I;

	p	= fftwf_plan_dft_1d( n, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

	fftwf_execute(p);

	for( int i = 0 ; i < n ; ++ i )
		out1[i]	= cabsf(out[i]);
	write_1d_float_wb( out1, n, "fw.dft.dat");


/*
	p	= fftwf_plan_dft_1d( n, out, in, FFTW_BACKWARD, FFTW_ESTIMATE );

	fftwf_execute(p);

	for( int i = 0 ; i < n ; ++ i )
		printf("in[%d]=%f\n", i, crealf(in[i])/n);
*/

	fftwf_destroy_plan(p);
	free(in);
	free(out);

}

void	test_dft_r2r_1d( float *in, float complex *out, int n )
{

	fftwf_plan	p;

	p	= fftwf_plan_dft_r2c_1d( n, in, out, FFTW_ESTIMATE );

	fftwf_execute(p);

	for( int i = 0 ; i < n ; ++ i )
		printf("in[%d]=%f out[i]=%f,%f,%f\n",
			i, in[i], crealf(out[i])/sqrt(n), cimagf(out[i])/sqrt(n), cabsf(out[i])/sqrt(n));

/*
*/
	for( int i = 2 ; i < n ; ++ i )
		out[i]	= 0;

	p	= fftwf_plan_dft_c2r_1d( n, out, in, FFTW_ESTIMATE );

	fftwf_execute(p);

	for( int i = 0 ; i < n ; ++ i )
		printf("in[%d]=%f\n", i, (in[i])/n);

	fftwf_destroy_plan(p);

}

/*
int main( int argc , char *argv[] )
{
	float	*fx	= calloc(nn, sizeof(float));
	float	*fu	= calloc(nn, sizeof(float));
	float complex	*fxc	= calloc(nn, sizeof(float complex));
	float complex	*fuc	= calloc(nn, sizeof(float complex));

	for( int i = 0 ; i < nn ; ++ i )
	{
		fx[i]	= 8*(i+1);
		fxc[i]	= 8*(i+1) + 0*I;
	}

	//test_dct1d( fx, fu, nn );
	//test_dft1d( fxc, fuc, nn );
	test_dft_r2r_1d( fx, fuc, nn );

	free(fx);
	free(fu);
	free(fxc);
	free(fuc);

	return 0 ;
}

*/


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
/*
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

	test_dct2d(vel, velDCT, nx, nz);

	free2float(vel);
	free2float(velDCT);

	return 0 ;
}

*/
