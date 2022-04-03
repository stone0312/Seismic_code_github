
void	cal_z_1st_deri( float **in, float **out, int nx, int nz );
void	cal_x_1st_deri( float **in, float **out, int nx, int nz );
void	Ltrans( float **X, float **Px1, float **Pz1, float **fai, float **Px, float **Pz, int nx, int nz );
void	Lforward( float **X, float **Px1, float **Pz1, float **fai, float **Px, float **Pz, int nx, int nz );


void    cal_z_1st_deri1( float **in, float **out, int nx, int nz );
void    cal_x_1st_deri1( float **in, float **out, int nx, int nz );
void    Ltrans1( float **X,  float **Px, float **Pz, int nx, int nz );
void    Lforward1( float **X, float **Px, float **Pz, int nx, int nz );

void    cal_z_1st_deri1( float **in, float **out, int nx, int nz )
{
#pragma omp parallel for
        for ( int ix = 0 ; ix < nx ;   ++ ix )
        for ( int iz = 0 ; iz < nz-1 ; ++ iz )
                out[ix][iz]     = in[ix][iz] - in[ix][iz+1];
}

void    cal_x_1st_deri1( float **in, float **out, int nx, int nz )
{
#pragma omp parallel for
        for ( int iz = 0 ; iz < nz ;   ++ iz )
        for ( int ix = 0 ; ix < nx-1 ; ++ ix )
                out[ix][iz]     = in[ix][iz] - in[ix+1][iz];
}

void    Ltrans1( float **X, float **Px, float **Pz, int nx, int nz )
{
        cal_z_1st_deri1( X, Pz, nx, nz );

        cal_x_1st_deri1( X, Px, nx, nz );
}


void    Lforward1( float **X, float **Px, float **Pz, int nx, int nz )
{
#pragma omp parallel for
        for ( int iz = 0 ; iz < nz ; ++ iz )
        for ( int ix = 0 ; ix < nx ; ++ ix )
                X[ix][iz] = 0 ;

#pragma omp parallel for
        for ( int iz = 0 ; iz < nz-1 ; ++ iz )
        for ( int ix = 0 ; ix < nx   ; ++ ix )
                X[ix][iz]       = Pz[ix][iz] ;

#pragma omp parallel for
        for ( int iz = 0 ; iz < nz   ; ++ iz )
        for ( int ix = 0 ; ix < nx-1 ; ++ ix )
                X[ix][iz]       += Px[ix][iz] ;

#pragma omp parallel for
        for ( int iz = 0 ; iz < nz-1 ; ++ iz )
        for ( int ix = 0 ; ix < nx   ; ++ ix )
                X[ix][iz+1]     -= Pz[ix][iz] ;

#pragma omp parallel for
        for ( int iz = 0 ; iz < nz   ; ++ iz )
        for ( int ix = 0 ; ix < nx-1 ; ++ ix )
                X[ix+1][iz]     -= Px[ix][iz] ;

}



//	calculate the 1st-derivative along Z-axis.
//	in[nx][nz];
//	out[nx][nz-1];

void	cal_z_1st_deri( float **in, float **out, int nx, int nz )
{
#pragma omp parallel for
	for ( int ix = 0 ; ix < nx ;   ++ ix )
	for ( int iz = 0 ; iz < nz-1 ; ++ iz )
		out[ix][iz]	= in[ix][iz] - in[ix][iz+1];
        
//         printf ("2222\n");
}

//	calculate the 1st-derivative along X-axis.
//	in[nx][nz];
//	out[nx-1][nz];
void	cal_x_1st_deri( float **in, float **out, int nx, int nz )
{
#pragma omp parallel for
	for ( int iz = 0 ; iz < nz ;   ++ iz )
	for ( int ix = 0 ; ix < nx-1 ; ++ ix )
		out[ix][iz]	= in[ix][iz] - in[ix+1][iz];


//        printf ("2222\n");

}



//	X[nx][nz];
//	Pz[nx][nz-1];
//	Px[nx-1][nz];
void	Ltrans( float **X, float **Px1, float **Pz1, float **fai, float **Px, float **Pz,  int nx, int nz )
{
	cal_z_1st_deri( X, Pz1, nx, nz );

	cal_x_1st_deri( X, Px1, nx, nz );

        for(int ix=0; ix<nx-1; ix++)
           for(int iz=0; iz<nz; iz++)
             Px[ix][iz]=0.0;

        for(int ix=0; ix<nx; ix++)
           for(int iz=0; iz<nz-1; iz++)
             Pz[ix][iz]=0.0;

        for(int ix=0; ix<nx-1; ix++) 
           for(int iz=0; iz<nz-1; iz++)
             { 
                Px[ix][iz]=Px1[ix][iz]*cos(fai[ix][iz]/180.0*PI)-Pz1[ix][iz]*sin(fai[ix][iz]/180.0*PI);
                Pz[ix][iz]=Px1[ix][iz]*sin(fai[ix][iz]/180.0*PI)+Pz1[ix][iz]*cos(fai[ix][iz]/180.0*PI);
             }

//        printf ("2222\n");
}

void	Lforward( float **X, float **Px1, float **Pz1, float **fai, float **Px, float **Pz, int nx, int nz )
{
         for(int ix=0; ix<nx-1; ix++)
           for(int iz=0; iz<nz; iz++)
             Px[ix][iz]=0.0;

        for(int ix=0; ix<nx; ix++)
           for(int iz=0; iz<nz-1; iz++)
             Pz[ix][iz]=0.0;

	for ( int ix = 0 ; ix < nx-1 ; ++ ix )
	  for ( int iz = 0 ; iz < nz-1 ; ++ iz )
            {
               Px[ix][iz]=Px1[ix][iz]*cos(fai[ix][iz]/180.0*PI)+Pz1[ix][iz]*sin(fai[ix][iz]/180.0*PI);
               Pz[ix][iz]=-Px1[ix][iz]*sin(fai[ix][iz]/180.0*PI)+Pz1[ix][iz]*cos(fai[ix][iz]/180.0*PI);
            }    

#pragma omp parallel for
	for ( int iz = 0 ; iz < nz ; ++ iz )
	for ( int ix = 0 ; ix < nx ; ++ ix )
		X[ix][iz] = 0 ;
	
	// matlab codes.
	//X(1:m-1,:)=P{1};
	//X(:,1:n-1)=X(:,1:n-1)+P{2};
	//X(2:m,:)=X(2:m,:)-P{1};
	//X(:,2:n)=X(:,2:n)-P{2};

#pragma omp parallel for
	for ( int iz = 0 ; iz < nz-1 ; ++ iz )
	for ( int ix = 0 ; ix < nx   ; ++ ix )
		X[ix][iz]	= Pz[ix][iz] ;

#pragma omp parallel for
	for ( int iz = 0 ; iz < nz   ; ++ iz )
	for ( int ix = 0 ; ix < nx-1 ; ++ ix )
		X[ix][iz]	+= Px[ix][iz] ;
	
#pragma omp parallel for
	for ( int iz = 0 ; iz < nz-1 ; ++ iz )
	for ( int ix = 0 ; ix < nx   ; ++ ix )
		X[ix][iz+1]	-= Pz[ix][iz] ;

#pragma omp parallel for
	for ( int iz = 0 ; iz < nz   ; ++ iz )
	for ( int ix = 0 ; ix < nx-1 ; ++ ix )
		X[ix+1][iz]	-= Px[ix][iz] ;
	
}

// min{ || X - Xobs || ^ 2 + 2 * lambda * | B * X |1 }
//	Xobs[nx][nz];
int	l1norm_denoise_tv_reg( float **Xobs, float **X, float lambda, float niu,
		char *type,
		int nx, int nz, int nitermax, float epsmin, float *obj , float **fai, float **D1)
{
	int	Mystatus = -1;

	float	**D	= NULL;
	float	**A	= NULL;
	float	**Dold	= NULL;
	float	**BTpk	= NULL;
	float	**Px	= NULL;
	float	**Pz	= NULL;
	float	**Qx	= NULL;
	float	**Qz	= NULL;
	float	**Px1	= NULL;
	float	**Pz1	= NULL;
	float	**Qx1	= NULL;
	float	**Qz1	= NULL;

	// Initialize variables.
	D	= alloc2float(nz, nx);
	A	= alloc2float(nz, nx);
	Dold	= alloc2float(nz, nx);
	BTpk	= alloc2float(nz, nx);
	Px	= alloc2float(nz, nx-1);
	Pz	= alloc2float(nz-1, nx);
	Qx	= alloc2float(nz, nx-1);
	Qz	= alloc2float(nz-1, nx);
	Px1	= alloc2float(nz, nx-1);
	Pz1	= alloc2float(nz-1, nx);
	Qx1	= alloc2float(nz, nx-1);
	Qz1	= alloc2float(nz-1, nx);

	memset((void*)&D[0][0],  0, sizeof(float)*nx*nz);
	memset((void*)&Px[0][0], 0, sizeof(float)*(nx-1)*nz);
	memset((void*)&Pz[0][0], 0, sizeof(float)*(nz-1)*nx);

	fprintf(stderr, "\t***********************************\n");
	fprintf(stderr, "\tSolving with Proximity**\n");
	fprintf(stderr, "\t***********************************\n");
	fprintf(stderr, "\t#iteration  relative-difference\n");
	fprintf(stderr, "\t-----------------------------------\n");



/*        
	Ltrans1( Xobs,  Qx, Qz, nx, nz );
        Lforward1( D1, Qx, Qz, nx, nz );
        return  Mystatus; 
*/       
/*
	Ltrans( Xobs,  Qx1, Qz1, fai, Qx, Qz,  nx, nz );
        Lforward( D1, Qx, Qz, fai, Px, Pz, nx, nz );
        return  Mystatus; 
*/


	float	eps1	= 1E-9;
	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{

		// Dold <= D.
#pragma omp parallel for
		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
			Dold[ix][iz]	= D[ix][iz];

		Lforward( BTpk,Px1, Pz1, fai,  Px, Pz, nx, nz );

/*
        ofstream swq333;
        swq333.open("./test222.dat",ios::binary);
        for ( int ix = 0 ; ix < nx ; ix ++ )
           for ( int iz = 0 ; iz < nz ; iz ++ )
              swq333.write((char*)&BTpk[ix][iz],sizeof(BTpk[ix][iz]));
        return 0;
*/


#pragma omp parallel for
		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
			D[ix][iz]	= Xobs[ix][iz] - lambda*BTpk[ix][iz];

		Ltrans( D, Qx1, Qz1,fai, Qx, Qz, nx, nz );

#pragma omp parallel for
		for ( int ix = 0 ; ix < nx-1 ; ++ ix )
		for ( int iz = 0 ; iz < nz   ; ++ iz )
			Px1[ix][iz]	+= Qx[ix][iz] ;

#pragma omp parallel for
		for ( int ix = 0 ; ix < nx   ; ++ ix )
		for ( int iz = 0 ; iz < nz-1 ; ++ iz )
			Pz1[ix][iz]	+= Qz[ix][iz] ;

		if ( 0 == strcmp("iso", type) )
		{
			// calculate A = sqrt(px^2+pz^2);
#pragma omp parallel for
			for ( int ix = 0 ; ix < nx ; ++ ix )
			for ( int iz = 0 ; iz < nz ; ++ iz )
			{
				if ( ix < nx-1 && iz < nz-1 )
				{
					float	px	= Px[ix][iz];
					float	pz	= Pz[ix][iz];
					A[ix][iz]	= sqrt(px*px + pz*pz);
				}
				else if (ix == nx-1 && iz < nz-1)
					A[ix][iz]	= fabs(Pz[ix][iz]);
				else if (iz == nz-1 && ix < nx-1)
					A[ix][iz]	= fabs(Px[ix][iz]);
				else if (iz == nz-1 && ix == nx-1)
					A[ix][iz]	= 0;
			}

#pragma omp parallel for
			for ( int ix = 0 ; ix < nx-1 ; ++ ix )
			for ( int iz = 0 ; iz < nz   ; ++ iz )
			{
				float	px	= Px[ix][iz];
				float	aa	= A[ix][iz];
				//float	tmp	= MAX( aa-niu/lambda, 0 ) * px / aa;
				//Px[ix][iz]	= Px[ix][iz] - tmp ;
				float	tmp;
				if ( fabs(aa) < eps1 )
					tmp	= 0;
				else
					tmp	= MAX( aa-niu/lambda, 0 ) * px / aa;
				Px[ix][iz]	= tmp ;
			}

#pragma omp parallel for
			for ( int ix = 0 ; ix < nx   ; ++ ix )
			for ( int iz = 0 ; iz < nz-1 ; ++ iz )
			{
				float	pz	= Pz[ix][iz];
				float	aa	= A[ix][iz];
				//float	tmp	= MAX( aa-niu/lambda, 0 ) * pz / aa;
				//Pz[ix][iz]	= Pz[ix][iz] - tmp ;
				float	tmp;
				if ( fabs(aa) < eps1 )
					tmp	= 0;
				else
					tmp	= MAX( aa-niu/lambda, 0 ) * pz / aa;
				Pz[ix][iz]	= tmp ;
			}
		}
		else if ( 0 == strcmp("l1", type) )
		{
#pragma omp parallel for
			for ( int ix = 0 ; ix < nx-1 ; ++ ix )
			for ( int iz = 0 ; iz < nz   ; ++ iz )
			{
				float	px	= Px1[ix][iz];
				float	tmp	= MAX( fabs(px)-niu/lambda, 0 ) * SGN(px) ;
				Px1[ix][iz]	= Px1[ix][iz] - tmp ;
//				Px1[ix][iz]	= tmp ;
			}

/*
#pragma omp parallel for
			for ( int ix = 0 ; ix < nx   ; ++ ix )
			for ( int iz = 0 ; iz < nz-1 ; ++ iz )
			{
				float	pz	= Pz1[ix][iz];
				float	tmp	= MAX( fabs(pz)-niu/lambda, 0 ) * SGN(pz) ;
				Pz1[ix][iz]	= Pz1[ix][iz] - tmp ;
//				Pz1[ix][iz]	= tmp ;
			}
*/

		}
		else
		{
			fprintf(stderr, "tv type error.\n");
			fprintf(stderr, "type should be (l1) or (iso).\n");
			return -1;
		}

		// calculate the residuals.
		double	sum1	= 0;
		double	sum2	= 0;
		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
		{
			float	xnow	= D[ix][iz];
			float	xold	= Dold[ix][iz];
			sum1	+= (xnow-xold) * (xnow-xold) ;
			sum2	+= xnow * xnow ;
		}
		if ( 0 == sum2 )
			break;
		double	re	= sum1/sum2;
		obj[iter]	= log10(re);

		if ( re < epsmin )
		{
			Mystatus = 0;
			break;
		}

		fprintf(stderr,"res[%d]=%8.6e\n", iter, re);

	}

	for ( int ix = 0 ; ix < nx ; ++ ix )
	for ( int iz = 0 ; iz < nz ; ++ iz )
		X[ix][iz]	= D[ix][iz];
	
	free2float(Px);
	free2float(Pz);
	free2float(Qx);
	free2float(Qz);
	free2float(D);
	free2float(A);
	free2float(Dold);
	free2float(BTpk);

	return	Mystatus;
}

/*
// min{ || X - Xobs || ^ 2 + 2 * lambda * | B * X |1 }
//	Xobs[nx][nz];
int	l1norm_denoise_tv_reg_using_prox( float **Xobs, float **X, float lambda, float niu,
		int nx, int nz, int nitermax, float epsmin, float *obj )
{
	int	Mystatus = -1;
	float	beta;

	float	**D	= NULL;
	float	**Dold	= NULL;
	float	**BTpk	= NULL;
	float	**Px	= NULL;
	float	**Pz	= NULL;
	float	**Px1	= NULL;
	float	**Pz1	= NULL;
	float	**Qx	= NULL;
	float	**Qz	= NULL;

	// Initialize variables.
	D	= alloc2float(nz, nx);
	Dold	= alloc2float(nz, nx);
	BTpk	= alloc2float(nz, nx);
	Px	= alloc2float(nz, nx-1);
	Pz	= alloc2float(nz-1, nx);
	Px1	= alloc2float(nz, nx-1);
	Pz1	= alloc2float(nz-1, nx);
	Qx	= alloc2float(nz, nx-1);
	Qz	= alloc2float(nz-1, nx);

	memset((void*)&D[0][0],  0, sizeof(float)*nx*nz);
	memset((void*)&Px[0][0], 0, sizeof(float)*(nx-1)*nz);
	memset((void*)&Pz[0][0], 0, sizeof(float)*(nz-1)*nx);

	fprintf(stdout, "\t***********************************\n");
	fprintf(stdout, "\tSolving with Proximity**\n");
	fprintf(stdout, "\t***********************************\n");
	fprintf(stdout, "\t#iteration  relative-difference\n");
	fprintf(stdout, "\t-----------------------------------\n");

	beta	= niu/lambda;
	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{

		// Dold <= D.
		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
			Dold[ix][iz]	= D[ix][iz];

		Lforward( BTpk, Px, Pz, nx, nz );

		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
			D[ix][iz]	= Xobs[ix][iz] - lambda*beta*BTpk[ix][iz];

		Ltrans( D, Qx, Qz, nx, nz );

		for ( int ix = 0 ; ix < nx-1 ; ++ ix )
		for ( int iz = 0 ; iz < nz   ; ++ iz )
			Px[ix][iz]	+= Qx[ix][iz] ;

		for ( int ix = 0 ; ix < nx   ; ++ ix )
		for ( int iz = 0 ; iz < nz-1 ; ++ iz )
			Pz[ix][iz]	+= Qz[ix][iz] ;

		// denoise.
		sprox_l1norm_vector( &Px[0][0], &Px1[0][0], (nx-1)*nz, 1.0/beta);
		sprox_l1norm_vector( &Pz[0][0], &Pz1[0][0], (nz-1)*nx, 1.0/beta);
		
		// update.
		for ( int ix = 0 ; ix < nx-1 ; ++ ix )
		for ( int iz = 0 ; iz < nz   ; ++ iz )
		{
			Px[ix][iz]	= Px[ix][iz] - Px1[ix][iz] ;
		}
		for ( int ix = 0 ; ix < nx   ; ++ ix )
		for ( int iz = 0 ; iz < nz-1 ; ++ iz )
		{
			Pz[ix][iz]	= Pz[ix][iz] - Pz1[ix][iz] ;
		}

		// calculate the residuals.
		double	sum1	= 0;
		double	sum2	= 0;
		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
		{
			float	xnow	= D[ix][iz];
			float	xold	= Dold[ix][iz];
			sum1	+= (xnow-xold) * (xnow-xold) ;
			sum2	+= xnow * xnow ;
		}
		double	re	= sum1/sum2;
		obj[iter]	= log10(re);

		if ( re < epsmin )
		{
			Mystatus = 0;
			break;
		}

		printf("res[%d]=%8.6e\n", iter, re);

	}

	for ( int ix = 0 ; ix < nx ; ++ ix )
	for ( int iz = 0 ; iz < nz ; ++ iz )
		X[ix][iz]	= D[ix][iz];
	
	free2float(Px);
	free2float(Pz);
	free2float(Px1);
	free2float(Pz1);
	free2float(Qx);
	free2float(Qz);
	free2float(D);
	free2float(Dold);
	free2float(BTpk);

	return	Mystatus;
}
*/
