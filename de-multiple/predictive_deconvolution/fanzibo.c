#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int gaus( float *a,float *b,int n );	//Using Gaussian-Seidel Iteration to compute the inverse of the matrix A.
int convolution(float *x,int nx,float *y,int ny,float *z);//compute the convolution of series x and y.
void ACF(float *x,float *y,int n );//compute the Autocorrelation Function of series x.
void Todiagonal( float *a , float *r , int n , int length , int lag , float baizao );
void Todiagonal( float *a , float *r , int n , int length , int lag , float baizao )
{
	int i,j,p;
	for( i=0;i<n;i++ )
	{
		for( j=0;j<n;j++ )
		{
			if( abs(lag+i-j)<length )
			{
				p=i*n+j;
				*(a+p)=r[abs(lag+i-j)];
			}
			else
			{
				p=i*n+j;
				*(a+p)=0.0;
			}
		}
	}
	for( i=0;i<n;i++ )
	{
		p=i*n+i;
		*(a+p)=(*(a+p))*(1+baizao);
	}
}
int convolution(float *x,int nx,float *y,int ny,float *z)
{
	int i,j;
	float tmp1[nx];
	float tmp2[ny];
	
	for(i=0;i<nx+ny-1;i++)
	{
		z[i]=0.0;
	}
	
	if(nx>ny)
	{
		for(i=0;i<ny;i++)
		{
			tmp2[i]=y[ny-1-i];
		}
		for(i=0;i<ny;i++)
		{
			for(j=0;j<=i;j++)
			{
				z[i]=z[i]+x[j]*tmp2[ny-1-i+j];
			}
		}
		for(i=ny;i<nx;i++)
		{
			for(j=0;j<ny;j++)
			{
				z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
			}
		}
		for(i=nx;i<nx+ny-1;i++)
		{
			for(j=0;j<ny-i+nx-1;j++)
			{
				z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
			}
		}
	}

	else if(nx==ny)
	{
		for(i=0;i<ny;i++)
		{
			tmp2[i]=y[ny-1-i];
		}
		for(i=0;i<ny;i++)
		{
			for(j=0;j<=i;j++)
			{
				z[i]=z[i]+x[j]*tmp2[ny-1-i+j];
			}
		}
		for(i=ny;i<2*ny-1;i++)
		{
			for(j=0;j<2*ny-i-1;j++)
			{
				z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
			}
		}
	}

	else
	{
		for(i=0;i<nx;i++)
		{
			tmp1[i]=x[nx-1-i];
		}
		for(i=0;i<nx;i++)
		{
			for(j=0;j<=i;j++)
			{
				z[i]=z[i]+y[j]*tmp1[nx-1-i+j];
			}
		}
		for(i=nx;i<ny;i++)
		{
			for(j=0;j<nx;j++)
			{
				z[i]=z[i]+tmp1[j]*y[i-nx+1+j];
			}
		}
		for(i=ny;i<nx+ny-1;i++)
		{
			for(j=0;j<nx-i+ny-1;j++)
			{
				z[i]=z[i]+tmp1[j]*y[i-nx+1+j];
			}
		}
	}
	return 1;
}
void ACF(float *x,float *y,int n)
{
	int i,j;
	float sum;
	for(i=0;i<n;i++)
	{
		y[i]=0.0;
	}

	for(i=0;i<n;i++)
	{
		sum=0;
		for(j=i;j<n;j++)
		{
			sum=x[j]*x[j-i]+sum;
		}
		y[i]=sum;
	}
	
}
int gaus( float *a,float *b,int n )
{
	int *js,l=1,k,i,j,is,p,q;
	float d,t;

	js=malloc( n*sizeof(int) );
	for( k=0;k<=n-2;k++ )
	{
		d=0.0;
		for( i=k;i<=n-1;i++ )
		for (j=k;j<=n-1;j++)
		{
			t=fabs(a[i*n+j]);	
			if(t>d)
			{
				d=t;
				js[k]=j;
				is=i;
			}
		}
		if(d+1.0==1.0)
			l=0;
		else
		{
			if (js[k]!=k)
				for(i=0;i<=n-1;i++)
				{
					p=i*n+k;
					q=i*n+js[k];
					t=a[p];
					a[p]=a[q];
					a[q]=t;

				}
	
			if(is!=k)
			{
				for(j=k;j<=n-1;j++)
				{
					p=k*n+j;
					q=is*n+j;
					t=a[p];
					a[p]=a[q];
					a[q]=t;

				}

				t=b[k];
				b[k]=b[is];
				b[is]=t;

			}
		}

		if(l==0)
		{
			free(js);
			printf("fail\n");
			return(0);

		}
		d=a[k*n+k];
		for( j=k+1;j<=n-1;j++ )
		{
			p=k*n+j;
			a[p]=a[p]/d;
		}
		b[k]=b[k]/d;
		for( i=k+1;i<=n-1;i++ )
		{
			for( j=k+1;j<=n-1;j++ )
			{
				p=i*n+j;
				a[p]=a[p]-a[i*n+k]*a[k*n+j];
			}
			b[i]=b[i]-a[i*n+k]*b[k];
		}
	}
	d=a[(n-1)*n+n-1];
	if(fabs(d)+1.0==1.0)
	{
		free(js);
		printf("fail\n");
		return(0);
	}
	b[n-1]=b[n-1]/d;
	for( i=n-2;i>=0;i-- )
	{
		t=0.0;
		for( j=i+1;j<=n-1;j++ )
			t=t+a[i*n+j]*b[j];
		b[i]=b[i]-t;
	}
	js[n-1]=n-1;
	for( k=n-1;k>=0;k-- )
	if( js[k]!=k )
	{
		t=b[k];
		b[k]=b[js[k]];
		b[js[k]]=t;
	}
	free(js);
	return(1);
}

int main( void )
{
	float *f,*acf_f,**a,epsilon,*inverse,*solu_conv;
	int i,j,leng_wav,leng_inv,flag;
	FILE *fp1,*fp2,*fp3;

/*********Please input the parameter of the program which include length of wavelet and**********
**********inverse wavelet respectively, and the epsilon respect the white noise should added****** 
*********in the digonal of the matrix A.***********************************************************/

	printf( "Please input the length of wavelet!" );
	scanf( "%d" , &leng_wav );
	printf( "Please input the length of the inverse wavelet!" );
	scanf( "%d" , &leng_inv );
	printf( "Please input the white noise factor! It could be 0.1, for example!" );
	scanf( "%f" , &epsilon );

	f=( float * )malloc( leng_wav*sizeof(float) );
	acf_f=( float *)malloc( leng_wav*sizeof(float) );
	inverse=( float *)malloc( leng_inv*sizeof(float) );
	solu_conv=( float * )malloc( (leng_inv+leng_wav-1)*sizeof(float) );
	a=( float ** )malloc( leng_inv*sizeof(float *) );
	for( i=0;i<leng_inv;i++ )
		a[i]=( float * )malloc( leng_inv*sizeof(float) );

	for( i=0;i<leng_wav;i++ )
	{
		f[i]=0.0;
		acf_f[i]=0.0;
	}
	for( i=0;i<leng_inv;i++ )
	{
		if( i==0 )
			inverse[i]=1.0;
		else
			inverse[i]=0.0;
		for( j=0;j<leng_inv;j++ )
			a[i][j]=0.0;
	}
// ********Read the data of wavelet  ***************** //
	fp1=fopen( "wavelet.dat" , "rb" );

	for( i=0;i<leng_wav;i++ )
		fread( &f[i] , sizeof(float) , 1 , fp1 );
	fclose( fp1 );

// *********compute the ACF of the wavelet and compute the inverse of the matrix A ** //
	ACF( f , acf_f , leng_wav );
	Todiagonal( &a[0][0] , acf_f , leng_inv , leng_wav , 0 , epsilon );
	flag=gaus( &a[0][0] , inverse , leng_inv );
	fp2=fopen( "inverse.dat" , "wb" );
	for( i=0;i<leng_inv;i++ )
		fwrite( &inverse[i] , sizeof(float) , 1 , fp2 );
	fclose( fp2 );
// ***********sure the inverse of the wavelet is right ***************************** //
	flag=convolution( inverse , leng_inv , f , leng_wav , solu_conv );
	fp3=fopen( "solu_conv.dat" , "wb" );
	for( i=0;i<leng_inv+leng_wav-1;i++ )
		fwrite( &solu_conv[i] , sizeof(float) , 1 , fp3 );
	fclose( fp3 );
// *******************Free the Allays ********************************************* //
	free( f );
	free( acf_f );
	free( inverse );
	free( solu_conv );
	for( i=0;i<leng_inv;i++ )
		free(a[i]);
	free( a );

	return 1;
}

