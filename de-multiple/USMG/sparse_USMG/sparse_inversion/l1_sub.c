
//void	Ltrans( float **X, float **Px, float **Pz, int nx, int nz );
//void	Lforward( float **X, float **Px, float **Pz, int nx, int nz );

void Ltrans(float **u,float **b,float **c,float *p,float *x,int *interc,float **sum1,float **sem,int ntrace,int lt,int np,int N,float theata_max,float theata_min,float dt,float dx,float v,float threhold,int win_time);
void  Lforward(float **u,float *ufinal,int lt,int np);


int Ltrans(float **u,float **b,float **c,float *p,float *x,int *interc,float **sum1,float **sem,int ntrace,int lt,int np,int N,float theata_max,float theata_min,float dt,float dx,float v,float threhold,int win_time)
{
   int i,j,k,l,ll,ii,jj,kk;

         float **sum3=new float *[np];
         for(i=0;i<np;i++)
                 sum3[i]=new float [lt];

         float **sum4=new float *[np];
         for(i=0;i<np;i++)
                 sum4[i]=new float [lt];

         float *tmp=new float [win_time*2+1];
         float *tmp1=new float [win_time*2+1];

         for(i=0;i<win_time*2+1;i++)
           {
                          tmp[i]=0.0;
                                tmp1[i]=0.0;
           }
         for(i=0;i<np;i++)
                  for(j=0;j<lt;j++)
                        {
                                sum3[i][j]=0.0;
                                sum4[i][j]=0.0;
                        }

   for(i=0;i<np;i++)
      {
            for(j=0;j<lt;j++)
                {
                    for(k=0;k<N+1;k++)
                       interc[k]=int(j+1000*p[i]*x[k]/dt+0.5);

                    for(ii=0;ii<N+1;ii++)
                       {
                           kk=interc[ii];
                           if(kk<0)
                                 sum1[i][j]+=0.0;
                           else if(kk>lt)
                                 sum1[i][j]+=0.0;
                           else
                                 sum1[i][j]+=u[ii][kk];
                        }
                }
            for(j=win_time;j<lt-win_time;j++)
                {
                   k=j-win_time;
                   for(l=k;l<k+win_time*2+1;l++)
                      {
                         for(ll=0;ll<N+1;ll++)
                            interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);

                         for(ii=0;ii<N+1;ii++)
                             {
                                kk=interc[ii];
                                if(kk<0)
                                   {
                                     tmp[l-k]+=0.0;
                                     tmp1[l-k]+=0.0;
                                   }
                                else if(kk>lt)
                                   {
                                        tmp[l-k]+=0.0;
                                        tmp1[l-k]+=0.0;
                                   }
                                else
                                   {
                                        tmp[l-k]+=u[ii][kk];
                                        tmp1[l-k]+=pow(u[ii][kk],2);
                                   }
                             }
                       }
                   for(jj=0;jj<win_time*2+1;jj++)
                      {
                        sum3[i][j]+=tmp[jj];
                        sum4[i][j]+=tmp1[jj];
                      }
                   for(jj=0;jj<win_time*2+1;jj++)
                      {
                        tmp[jj]=0.0;
                        tmp1[jj]=0.0;
                      }
               }
          for(j=0;j<win_time;j++)
              {
                for(l=0;l<win_time*2+1;l++)
                  {
                    for(ll=0;ll<N+1;ll++)
                      interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                    for(ii=0;ii<N+1;ii++)
                      {
                         kk=interc[ii];
                         if(kk<0)
                           {
                             tmp[l]+=0.0;
                             tmp1[l]+=0.0;
                           }
                         else if(kk>lt)
                           {
                              tmp[l]+=0.0;
                              tmp1[l]+=0.0;
                           }
                        else
                           {
                              tmp[l]+=u[ii][kk];
                              tmp1[l]+=pow(u[ii][kk],2);
                           }
                      }
                   }
                  for(jj=0;jj<win_time*2+1;jj++)
                   {
                     sum3[i][j]+=tmp[jj];
                     sum4[i][j]+=tmp1[jj];
                   }
                 for(jj=0;jj<win_time*2+1;jj++)
                   {
                     tmp[jj]=0.0;
                     tmp1[jj]=0.0;
                   }
              }
           for(j=lt-win_time;j<lt;j++)
              {
                 for(l=lt-2*win_time-1;l<lt;l++)
                    {
                        for(ll=0;ll<N+1;ll++)
                          interc[ll]=int(j+1000*p[i]*x[ll]/dt+0.5);
                        for(ii=0;ii<N+1;ii++)
                           {
                              kk=interc[ii];
                              if(kk<0)
                                {
                                  tmp[l-lt+2*win_time+1]+=0.0;
                                  tmp1[l-lt+2*win_time+1]+=0.0;
                                }
                              else if(kk>lt)
                                {
                                   tmp[l-lt+2*win_time+1]+=0.0;
                                   tmp1[l-lt+2*win_time+1]+=0.0;
                                }
                              else
                                {
                                   tmp[l-lt+2*win_time+1]+=u[ii][kk];
                                   tmp1[l-lt+2*win_time+1]+=pow(u[ii][kk],2);
                                }
                            }
                     }
                 for(jj=0;jj<win_time*2+1;jj++)
                     {
                        sum3[i][j]+=tmp[jj];
                        sum4[i][j]+=tmp1[jj];
                     }
                for(jj=0;jj<win_time*2+1;jj++)
                     {
                        tmp[jj]=0.0;
                        tmp1[jj]=0.0;
                     }
              }
      }

   for(i=0;i<np;i++)
      for(j=0;j<lt;j++)
        sem[i][j]=pow(sum3[i][j],2)/(sum4[i][j]+0.01);
   for(i=0;i<np;i++)
      for(j=0;j<lt;j++)
        sem[i][j]/=((N+1)*(2*win_time+1));
     for(i=0;i<np;i++)
        for(j=0;j<lt;j++)
           c[i][j]=sem[i][j];

         for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          {
              if(c[i][j]<threhold)
                 sem[i][j]=0.0;
              else
                 sem[i][j]=1.0;
          }

   for(i=0;i<np;i++)
       for(j=0;j<lt;j++)
          b[i][j]=sum1[i][j]*sem[i][j];

         return 0;
}

int Lforward(float **u,float *ufinal,int lt,int np)
{
    int i,j;

    for(i=0;i<lt;i++)
      ufinal[i]=0.0;

    for(i=0;i<lt;i++)
       {
           for(j=0;j<np;j++)
               ufinal[i]+=u[j][i];
       }
    return 0;
}

// min{ || X - Xobs || ^ 2 + 2 * lambda * | B * X |1 }
//	Xobs[nx][nz];
int	l1norm_denoise_tv_reg( float **Xobs, float **X, float lambda, float niu, char *type, int nx, int nz, int nitermax, float epsmin, float *obj , 
                               float **c, float *p, float *x, float *interc, float **sum,float ** sem, int  ntrace, int lt, int np, int N, 
                               float theata_max, float theata_min, float dt, float dx, float v, float threhold, int win_time)
{
	int	Mystatus = -1;

	float	**D	= NULL;
	float	**A	= NULL;
	float	**Dold	= NULL;
	float	**BTpk	= NULL;
	float	**Px	= NULL;
	float	**Qx	= NULL;

	// Initialize variables.
	D	= alloc2float(nz, nx);
	A	= alloc2float(nz, nx);
	Dold	= alloc2float(nz, nx);
	BTpk	= alloc2float(nz, nx);
	Px	= alloc2float(nz, np);
	Qx	= alloc2float(nz, np);

	float	*x	= NULL;
	float	*p	= NULL;
	float	*interc	= NULL;
	float	**sum	= NULL;
	float	**sem	= NULL;
	float	**c	= NULL;
       
        x=alloc1float(N+1);
        interc=alloc1float(N+1);
        p=alloc1float(np);
        sum=alloc2float(nz,np);
        sem=alloc2float(nz,np);
        c=alloc2float(nz,np);
        
	memset((void*)&sum[0][0], 0, sizeof(float)*np*nz);
	memset((void*)&sem[0][0], 0, sizeof(float)*np*nz);
	memset((void*)&c[0][0], 0, sizeof(float)*np*nz);
	memset((void*)&D[0][0],  0, sizeof(float)*nx*nz);
	memset((void*)&Px[0][0], 0, sizeof(float)*np*nz);

	fprintf(stdout, "\t***********************************\n");
	fprintf(stdout, "\tSolving with Proximity**\n");
	fprintf(stdout, "\t***********************************\n");
	fprintf(stdout, "\t#iteration  relative-difference\n");
	fprintf(stdout, "\t-----------------------------------\n");

	float	eps1	= 1E-9;
	for ( int iter = 1 ; iter <= nitermax ; ++ iter )
	{

		// Dold <= D.
		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
			Dold[ix][iz]	= D[ix][iz]; //x-t domain raw data


//inverse transform tp-->xt
                Lforward( BTpk, Px ,nz , np);

		for ( int ix = 0 ; ix < nx ; ++ ix )
		for ( int iz = 0 ; iz < nz ; ++ iz )
		     D[ix][iz]	= Xobs[ix][iz] - lambda*BTpk[ix][iz];

//forward transform xt-->tp
                Ltrans(D, Qx, c, p, x, interc, sum, sem, ntrace, lt, np, N, theata_max, theata_min, dt, dx, v, threhold, win_time);

		for ( int ix = 0 ; ix < nx   ; ++ ix )
		for ( int iz = 0 ; iz < nz   ; ++ iz )
			Px[ix][iz]	+= Qx[ix][iz] ;

			for ( int ix = 0 ; ix < nx ; ++ ix )
			for ( int iz = 0 ; iz < nz ; ++ iz )
			{
				float	px	= Px[ix][iz];
				float	tmp	= MAX( fabs(px)-niu/lambda, 0 ) * SGN(px) ;
				//Px[ix][iz]	= Px[ix][iz] - tmp ;
				Px[ix][iz]	= tmp ;
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

		printf("res[%d]=%8.6e\n", iter, re);

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
	free2float(c);
	free2float(sum);
	free2float(sem);
	free1float(p);
	free1float(interc);
	free1float(x);

	return	Mystatus;
}

