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

