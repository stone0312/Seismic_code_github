using namespace std;
#include <iostream>
#include <fstream>
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

#include "fballoc.h"
#include "fballoc.c"

#include "prox_sub.c"
#include "tv_sub_v2.omp.c"



int     readpar(char *parfn, char *vtruefile, char *vobsfile, char *multiplefile, char *vinvfile, char *vresfile, char *faifile,char *evofile, char *evosubfile,int *nx, int *nz, float * dx, float * dz, int *nitermax, float *lambda, float *niu, float *epsmin,char *type ,float *step);

void	normalize_data( int nx, int nz, float **data, float amp );

int orientation_smooth_diffusion(float **res, float **mul , float **evo, float **evo_old, float **mulx, float **mult, float **tensor, float ** evox, float **evot, float **dm1, float **dm2, float **dm1x, float **dm2t,  float **div, float step, int nx, int lt, float dx, float dt, int nitermax, float epsmin)
{
   int ix, it,ix1,it1,iter;
   float sum1,sum2;

//calculted the first deritive of x and z
   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++)
        {
	  mulx[ix][it]=0.0;
          mult[ix][it]=0.0;
        }
     
   for(it=0;it<lt;it++)
      for(ix=0;ix<nx-1;ix++)
         mulx[ix][it]=(mul[ix][it]-mul[ix+1][it])/dx;     

   for(it=0;it<lt;it++)
         mulx[nx-1][it]=-(mul[nx-2][it]-mul[nx-1][it])/dx;     

   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt-1;it++)
         mult[ix][it]=(mul[ix][it]-mul[ix][it+1])/dt;     

   for(ix=0;ix<nx;ix++)
         mult[ix][lt]=-(mul[ix][lt-2]-mul[ix][lt-1])/dt;     

//start the evolution

 for(ix=0;ix<nx;ix++)
   for(it=0;it<lt;it++)
      evo[ix][it]=res[ix][it];

 for(iter=0;iter<nitermax;iter++)
  {

   for(ix=0;ix<nx;ix++)
    for(it=0;it<lt;it++)
      evo_old[ix][it]=evo[ix][it];

   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++)
       {
//calculate the structor tensor
           tensor[0][0]=pow(mulx[ix][it],2); 
           tensor[0][1]=mulx[ix][it]*mulx[ix][it]; 
           tensor[1][0]=mulx[ix][it]*mulx[ix][it]; 
           tensor[1][1]=pow(mult[ix][it],2); 

//calculate the deritive of the iter-th iteration---evo
           if(ix==(nx-1))
             evox[ix][it]=evo[ix-1][it]-evo[ix][it]; 
           else
             evox[ix][it]=evo[ix][it]-evo[ix+1][it]; 

           if(it==(lt-1))
             evot[ix][it]=evo[ix][it-1]-evo[ix][it]; 
           else
             evot[ix][it]=evo[ix][it]-evo[ix][it+1]; 

           dm1[ix][it]=tensor[0][0]*evox[ix][it]+tensor[0][1]*evot[ix][it];
           dm2[ix][it]=tensor[1][0]*evox[ix][it]+tensor[1][1]*evot[ix][it];
 
       } 

//calculate the diversity
    for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++) 
       {
	  dm1x[ix][it]=0.0;
          dm2t[ix][it]=0.0;
       }

   for(it=0;it<lt;it++)
      for(ix=0;ix<nx-1;ix++)
         dm1x[ix][it]=(dm1[ix][it]-dm1[ix+1][it])/dx;     
   for(it=0;it<lt;it++)
         dm1x[nx-1][it]=-(dm1[nx-2][it]-dm1[nx-1][it])/dx;     

   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt-1;it++)
         dm2t[ix][it]=(dm2[ix][it]-dm2[ix][it+1])/dt;     
   for(ix=0;ix<nx;ix++)
         dm2t[ix][lt]=-(dm2[ix][lt-2]-dm2[ix][lt-1])/dt;     

    for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++) 
        div[ix][it]=dm1x[ix][it]+dm2t[ix][it];

    for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++)
        evo[ix][it]+=step*div[ix][it];

    sum1=0.0;
    sum2=0.0;
    for(ix=0;ix<nx;ix++)
      for(it=0;it<lt;it++) 
        {
           sum1+=(evo[ix][it]-evo_old[ix][it])*(evo[ix][it]-evo_old[ix][it]);
           sum2+=evo[ix][it]*evo[ix][it];
        }

      float err;
      err=sum1/(sum2+0.01);

      cout<<"Evolution Time is===="<<iter+1<<",Error is===="<<err<<endl;  

    if(sum1/sum2<epsmin)
      break;

   }

 return 0;

}


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
	char	multiplefile[FILE_NAME_MAX_LENGTH]="" ;
	char	vinvfile[FILE_NAME_MAX_LENGTH]="" ;
	char	vresfile[FILE_NAME_MAX_LENGTH]="" ;
	char	faifile[FILE_NAME_MAX_LENGTH]="" ;
	char	evofile[FILE_NAME_MAX_LENGTH]="" ;
	char	evosubfile[FILE_NAME_MAX_LENGTH]="" ;
	char	type[FILE_NAME_MAX_LENGTH]="" ;

	int	nx, nz;
        float   dx, dz;
	float	lambda, niu;
	int	nitermax;
	float	epsmin,step;
	float	*obj;

	// init.
	readpar( parfile, vtruefile, vobsfile, multiplefile, vinvfile, vresfile, faifile,evofile, evosubfile,  
		&nx, &nz, &dx, &dz, &nitermax, &lambda, &niu, &epsmin, type, &step);


	float	**vtrue	= NULL;
	float	**vobs	= NULL;
	float	**vinv	= NULL;
	float	**vres1	= NULL;
	float	**vres2	= NULL;
	float	**fai	= NULL;
	vtrue	= alloc2float(nz, nx);
	vobs	= alloc2float(nz, nx);
	vinv	= alloc2float(nz, nx);
	vres1	= alloc2float(nz, nx);
	vres2	= alloc2float(nz, nx);
	fai	= alloc2float(nz, nx);
        obj	= alloc1float(nitermax);

	float	**p	= NULL;
	float	**d	= NULL;
        float   **D1    = NULL;
	float	**mul	= NULL;
	float	**evo	= NULL;
	float	**evo_old	= NULL;
	float	**mulx	= NULL;
	float	**mult	= NULL;
	float	**tensor	= NULL;
	float	**evox	= NULL;
	float	**evot	= NULL;
	float	**dm1	= NULL;
	float	**dm2	= NULL;
	float	**dm1x	= NULL;
	float	**dm2t	= NULL;
	float	**div	= NULL;

	p	= alloc2float(nz, nx);
	d	= alloc2float(nz, nx);
        D1      = alloc2float(nz, nx);
	mul	= alloc2float(nz, nx);
	evo	= alloc2float(nz, nx);
	evo_old	= alloc2float(nz, nx);
	mulx	= alloc2float(nz, nx);
	mult	= alloc2float(nz, nx);
	tensor	= alloc2float(2, 2);
	evox	= alloc2float(nz, nx);
	evot	= alloc2float(nz, nx);
	dm1	= alloc2float(nz, nx);
	dm2	= alloc2float(nz, nx);
	dm1x	= alloc2float(nz, nx);
	dm2t	= alloc2float(nz, nx);
	div	= alloc2float(nz, nx);
        
//        printf("2222\n");

	// read data.
//	read_2d_float_rb(vtrue, nx, nz, vtruefile);
//	read_2d_float_rb(vobs,  nx, nz, vobsfile);

        ifstream swq1;
        swq1.open(vtruefile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq1.read((char*)&vtrue[ix][iz],sizeof(vtrue[ix][iz]));

        ifstream swq2;
        swq2.open(vobsfile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq2.read((char*)&vobs[ix][iz],sizeof(vobs[ix][iz]));
 
        ifstream swq22;
        swq22.open(faifile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq22.read((char*)&fai[ix][iz],sizeof(fai[ix][iz]));

        ifstream swq222;
        swq222.open(multiplefile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq222.read((char*)&mul[ix][iz],sizeof(mul[ix][iz]));


        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
               {
                  p[ix][iz]=0.0;
                  d[ix][iz]=0.0;
               }

        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              {
                 if(fai[ix][iz]==0.0)
                   p[ix][iz]=vobs[ix][iz];
              }
          

        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              {
                 if(fai[ix][iz]!=0.0)
                   d[ix][iz]=vobs[ix][iz];
              }


//        printf("2222\n");


	int	norm_flag = 1;
	float	amp	= 1.0;
	if ( 1 == norm_flag )
	{
	  normalize_data(nx, nz, vtrue, amp);
	  normalize_data(nx, nz, vobs,  amp);
	}

	// sparse inversion.
	l1norm_denoise_tv_reg( d, vinv, lambda, niu, type, nx, nz, nitermax, epsmin, obj,fai,D1);
/*	
        ofstream swq33;
        swq33.open("/data1/swq/structure_tensor/test1.dat",ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq33.write((char*)&D1[ix][iz],sizeof(D1[ix][iz]));

        return 0;
*/


//        printf("2222\n");


	for ( int ix = 0 ; ix < nx ; ++ ix )
	for ( int iz = 0 ; iz < nz ; ++ iz )
	{
//		vres1[ix][iz]	= vobs[ix][iz] - vtrue[ix][iz] ;
//		vres2[ix][iz]	= vinv[ix][iz] - vtrue[ix][iz] ;
		vres1[ix][iz]	= d[ix][iz] - vtrue[ix][iz] ;
		vres2[ix][iz]	= d[ix][iz] - vinv[ix][iz] ;
	}


        for ( int ix = 0 ; ix < nx ; ++ ix )
        for ( int iz = 0 ; iz < nz ; ++ iz )
          {
            vinv[ix][iz] += p[ix][iz];
//            vres2[ix][iz] += p[ix][iz];
          }
	// output.
//	write_2d_float_wb(vinv,  nx, nz, vinvfile);
//	write_2d_float_wb(vres2, nx, nz, vresfile);
//	write_1d_float_wb(obj, nitermax, "./obj.dat");

        ofstream swq3;
        swq3.open(vinvfile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq3.write((char*)&vinv[ix][iz],sizeof(vinv[ix][iz]));

        ofstream swq4;
        swq4.open(vresfile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq4.write((char*)&vres2[ix][iz],sizeof(vres2[ix][iz]));

        ofstream swq5;
        swq5.open("./objective.dat",ios::binary);
        for ( int ix = 0 ; ix < nitermax ; ix ++ )
              swq5.write((char*)&obj[ix],sizeof(obj[ix]));


        orientation_smooth_diffusion(vres2, mul , evo, evo_old, mulx, mult, tensor,  evox, evot, dm1, dm2, dm1x, dm2t,  div,  step,  nx, nz,  dx, dz, nitermax,  epsmin);

        ofstream swq6;
        swq6.open(evofile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq6.write((char*)&evo[ix][iz],sizeof(evo[ix][iz]));

        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
               p[ix][iz]=vobs[ix][iz]-evo[ix][iz];

        ofstream swq7;
        swq7.open(evosubfile,ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq7.write((char*)&p[ix][iz],sizeof(p[ix][iz]));

	// free pointers.
	free2float(vtrue);
	free2float(vobs);
	free2float(vinv);
	free2float(vres1);
	free2float(vres2);
	free(obj);

        printf("Done!\n");

	return 0;
}

int     readpar(char *parfn, char *vtruefile, char *vobsfile, char *multiplefile, char *vinvfile, char *vresfile, char *faifile,char *evofile, char *evosubfile,int *nx, int *nz, float * dx, float * dz, int *nitermax, float *lambda, float *niu, float *epsmin,char *type ,float *step)
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
	fscanf(fp, "%s%s", tmp, multiplefile);
	fscanf(fp, "%s%s", tmp, vinvfile);
	fscanf(fp, "%s%s", tmp, vresfile);
	fscanf(fp, "%s%s", tmp, faifile);
	fscanf(fp, "%s%s", tmp, evofile);
	fscanf(fp, "%s%s", tmp, evosubfile);
	fscanf(fp, "%s%d", tmp, nx);
	fscanf(fp, "%s%d", tmp, nz);
	fscanf(fp, "%s%f", tmp, dx);
	fscanf(fp, "%s%f", tmp, dz);
	fscanf(fp, "%s%d", tmp, nitermax);
	fscanf(fp, "%s%f", tmp, lambda);
	fscanf(fp, "%s%f", tmp, niu);
	fscanf(fp, "%s%f", tmp, epsmin);
	fscanf(fp, "%s%s", tmp, type);
	fscanf(fp, "%s%f", tmp, step);

        cout<<type<<endl;

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
