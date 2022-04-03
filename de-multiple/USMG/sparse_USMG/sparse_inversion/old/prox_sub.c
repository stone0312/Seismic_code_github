
//prototypes.
double	inner_product_vector_1d_f( float *a, float *b, int n );
void	sprox_l1norm_vector( float *in, float *out, int n, float gama );
void	sprox_l1norm_vector_dual( float *in, float *out, int n, float gama );
int	l1norm_inversion_denoise( float *Xobs, float *X, float lambda,
		int n, int nitermax, float epsmin );

double	inner_product_vector_1d_f( float *a, float *b, int n )
{
	double	sum	= 0 ;
	for ( int i = 0 ; i < n ; ++ i )
		sum	+= a[i] * b[i];
	return	sum;
}


/*
// initial s means single-precision, eg. float in C.
void	sprox_l1norm_vector_dual( float *in, float *out, int n, float gama )
{
	float	*tmp1	= calloc(n, sizeof(float));
	float	*tmp2	= calloc(n, sizeof(float));
	for ( int i = 0 ; i < n ; ++ i )
		tmp1[i]	= in[i]/gama;

	sprox_l1norm_vector(tmp1, tmp2, n, gama);

	for ( int i = 0 ; i < n ; ++ i )
		out[i]	= in[i] - gama * tmp2[i] ;

	free(tmp1);
	free(tmp2);
}
*/
//gama > 0.
// initial s means single-precision, eg. float in C.
void	sprox_l1norm_vector( float *in, float *out, int n, float gama )
{
	for ( int i = 0 ; i < n ; ++ i )
	{
		float   xi      = in[i];
		if ( xi > gama )
			out[i]  = xi - gama;
		else if ( fabs(xi) < gama )
			out[i]  = 0;
		else
			out[i]  = xi + gama;
	}
}


// initial s means single-precision, eg. float in C.
void	sprox_l1norm_vector_dual( float *in, float *out, int n, float gama )
{
	float	*tmp1	= calloc(n, sizeof(float));
	float	*tmp2	= calloc(n, sizeof(float));
	for ( int i = 0 ; i < n ; ++ i )
		tmp1[i]	= in[i]/gama;

	sprox_l1norm_vector(tmp1, tmp2, n, gama);

	for ( int i = 0 ; i < n ; ++ i )
		out[i]	= in[i] - gama * tmp2[i] ;

	free(tmp1);
	free(tmp2);
}

// min{ || X - Xobs || ^ 2 + 2 * lambda * | B * X |1 }
int	l1norm_inversion_denoise( float *Xobs, float *X, float lambda,
		int n, int nitermax, float epsmin )
{
	int	Mystatus = -1;
	float	*xtemp	= calloc( n, sizeof(float) );
	float	*resd	= calloc( n, sizeof(float) );
	float	*xold	= calloc( n, sizeof(float) );
	float	*xk	= calloc( n, sizeof(float) );
	float	*pk	= calloc( n, sizeof(float) );
	float	*ptemp	= calloc( n, sizeof(float) );
	float	beta;

	// Initialize variables.
	beta		= 1.0/lambda * 0.5;

	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{
		for ( int i = 0 ; i < n ; ++ i )
			xold[i] = xk[i];

		// 1. denoise.
		// x[k+1] = Xobs - lambda * p[k];
		for ( int i = 0 ; i < n ; ++ i )
			xk[i] = Xobs[i] - lambda * pk[i];

		// 2. ptemp = p[k] + beta * x[k+1];
		for ( int i = 0 ; i < n ; ++ i )
			ptemp[i] = pk[i] + beta * xk[i];

		// 3. update p[k+1];
		sprox_l1norm_vector_dual( ptemp, pk, n, beta );

		// 4. calculate the residuals.
		for ( int i = 0 ; i < n ; ++ i )
			resd[i]	= xk[i] - xold[i];
		double	sum1	= inner_product_vector_1d_f( resd, resd, n ) ;
		double	sum2	= inner_product_vector_1d_f( xk, xk, n ) ;
		double	re	= sum1/sum2;

		printf("res[%d]=%8.6lf\n", iter, re);

		if ( re < epsmin )
		{
			Mystatus = 0;
			break;
		}
	}

	for ( int i = 0 ; i < n ; ++ i )
		X[i] = xk[i];
	
	free(xtemp);
	free(resd);
	free(xold);
	free(xk);
	free(pk);
	free(ptemp);

	return	Mystatus;
}

/*
int	sparse_inversion_l1_denoise( float *dobs, float *ans, float lambda,
		int n, int nitermax, float epsmin )
{
	int	Mystatus = 0;
	float	steplength;
	float	*xtemp	= calloc( n, sizeof(float) );
	float	*xk	= calloc( n, sizeof(float) );
	float	*pk	= calloc( n, sizeof(float) );
	float	*ptemp	= calloc( n, sizeof(float) );
	float	beta;

	// Initialize variables.
	steplength	= 0.5;
	beta		= 1.0/lambda;
	for ( int i = 0 ; i < n ; ++ i )
		ans[i]	= 0;

	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{
		// 1. get xtemp, using the gradient and steplength.
		for ( int i = 0 ; i < n ; ++ i )
		{
			xtemp[i] = ans[i] - steplength*( ans[i] - dobs[i] );
		}

		// 2. denoise.
		// x[k+1] = xtemp - lambda * p[k];
		for ( int i = 0 ; i < n ; ++ i )
			xk[i] = xtemp[i] - lambda * pk[i];

		// ptemp = p[k] + beta * x[k+1];
		for ( int i = 0 ; i < n ; ++ i )
			ptemp[i] = pk[i] + beta * xk[i];

		// update p[k+1];
		sprox_l1norm_vector_dual( ptemp, pk, n, beta );

		// 3. calculate the residuals.
	}

}
*/
