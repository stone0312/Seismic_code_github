/*============================ The Tau-p Transform In The Frequence Domain ===================*/
int sis_RT_fre_2D(int nf, float df, int nx, float *cor_in, float complex **datainfx,
		          int np, float *p, float complex **dataoutfp)
{
	  int i, j, k;
  	float fre;
	  float cor_min, cor_mid, cor_max;
    float *cor;
  	float complex phase, **dataout_fp;

		/* allocate the memory */
	  cor= alloc1float(nx);
    zero1float(cor, nx);
    dataout_fp= alloc2floatcomplex(np, nf);

    cor_min= cor_in[0];
    cor_max= cor_in[0];

    for(i= 1; i< nx; i++)
	  {
	   	 if(cor_min> cor_in[i])cor_min= cor_in[i];
       if(cor_max< cor_in[i])cor_max= cor_in[i];
	  }

    cor_mid= 0.5*(cor_max + cor_min);

	  for(i= 0; i< nx; i++)
	  {
		   cor[i]= cor_in[i]- cor_mid;
	  }

    for(i= 0; i< nf; i++)
	  {
       fre= 2* PAI* df* i;
       for(j= 0; j< np; j++)
			 {
          dataout_fp[i][j]= 0;
          for(k= 0; k< nx; k++)
				  {
             phase= cos(fre* p[j]* cor[k])+
                  I* sin(fre* p[j]* cor[k]);
//           printf("fre=%f, cor=%f, p=%f, phase=%f\n",fre, cor[k], p[j], crealf(phase));
            dataout_fp[i][j]+= phase* datainfx[i][k];
          }
         dataoutfp[i][j]= dataout_fp[i][j]/ (float)nx;
      }
    }
    free1float(cor);
    free2floatcomplex(dataout_fp);

    return 1;
}

/*============== The Tau-p Transform In The Frequence Domain And Some Filter ===============*/
int sis_RT_HR_Time_2D(int ns, int nx, float dt, float *cor_in, float **data,
    int np, float *p, float **taup, float *pcycle)
{

    int i, j, k, l;
    double tmp;
    float df, nf_max;
    float cor_min, cor_max, cor_mid;
    float *cor, *coe, *sig_t;
    float **tr_in, **filt_max;
    float **tr_pt_filt;
    float complex *sig_f, **data_fx, **data_fp;

    //use 2-D array to access input data and out data
    df= 1.0/ns;
    nf_max=0.5*ns + 1;

    //alloc memory
    cor= alloc1float(nx);
    coe= alloc1float(ns);
    sig_t= alloc1float(ns);
    sig_f= alloc1floatcomplex(nf_max);

    tr_in= alloc2float(ns, nx);
    filt_max= alloc2float(ns, np);
    tr_pt_filt= alloc2float(ns, np);

    data_fx= alloc2floatcomplex(nx, nf_max);
    data_fp= alloc2floatcomplex(np, nf_max);

    for(i=0; i< nx; i++)
    {
       for(j= 0; j< ns; j++)
       {
          tr_in[i][j]= data[i][j];
       }
    }

    //renew the coordinate of input data
    cor_min= minn(cor_in, nx);
    cor_max= maxn(cor_in, nx);
    cor_mid= 0.5*(cor_max + cor_min);
    //printf("cor_min= %f, cor_max= %f, cor_mid=%f\n", cor_min, cor_max, cor_mid);

    for(i= 0; i< nx; i++)
		{
        cor[i]=cor_in[i] - cor_mid;
    }

    //create fftwf plan
    fftwf_plan t2f, f2t;
    t2f= fftwf_plan_dft_r2c_1d(ns, sig_t, (fftwf_complex*)sig_f, FFTW_ESTIMATE);
    f2t= fftwf_plan_dft_c2r_1d(ns, (fftwf_complex*)sig_f, sig_t, FFTW_ESTIMATE);

    //create a gaussian window function
    tmp=sqrtf(ns);
    for(i= 0; i< ns; i++)
    {
       coe[i]= exp(-pow(i/tmp, 2));
    }

    //change data to frequency-space domain
    memset(&data_fx[0][0], 0, nf_max* nx* sizeof(float complex));

    for(i= 0; i< nx; i++)
		{
        for(j= 0; j< ns; j++)sig_t[j] = tr_in[i][j]/ns;

        memset(sig_f, 0, nf_max*sizeof(float complex));
        fftwf_execute(t2f);

        for(j= 0; j< nf_max; j++)
				{
            data_fx[j][i]=sig_f[j];
        }
    }

/* get the information of data in frequency domain */
    double ener_total=0;
    for(i= 0; i< nf_max; i++)
		{
        for(j= 0; j< nx; j++)ener_total+= norm(data_fx[i][j]);
    }
    float const norm_coe= sqrtf(ener_total);

// normalizate the energy to 1
    for(i= 0; i< nf_max; i++)
    {
       for(j= 0; j< nx; j++)
       {
          data_fx[i][j]= data_fx[i][j]/ norm_coe;
       }
    }

    double const ener_out= (1- 1.0e-8);
    tmp= 0.5;
    ener_total= 0;

// get the half part energy and almost total energy num 
// nf_mid is the half  energy num
// nf     is the whole energy num

    for(i= 0; i< nf_max; i++)
		{
        for(j= 0; j< nx; j++)ener_total+= norm(data_fx[i][j]);
        if(ener_total>= tmp)break;
    }

    int const nf_mid= i;

    for(; i<nf_max; i++)
		{
       for(j= 0; j< nx; j++)ener_total+= norm(data_fx[i][j]);
       if(ener_total>= ener_out)break;
    }
    int const nf= i;

//  printf("nf_mid=%d,nf= %d\n", nf_mid, nf);

/*  begin to get the Radon spectrum of the input data */

    double filt_total= 0;
    double filt_now= 0;
    float max;
    float const dp= fabsf(p[1]- p[0]);
    i= 4.0/ (nf* df* fabsf(cor_max- cor_min)* dp);
    i= MAX(i, 5);
    int const p_len= MIN(i, (int)sqrtf(np));
    int t_len;
    int p1, p2;
    int t1, t2;

    // Change the data to F-P domain
    sis_RT_fre_2D(nf_max, df, nx, cor, data_fx, np, p, data_fp);

    // Change the data to P-T domain
    for(i= 0; i< np; i++)
    {
        memset(sig_f, 0, nf_max* sizeof(float complex));
        for(j= 0; j< nf_max; j++)sig_f[j]= data_fp[j][i];
        memset(sig_t, 0, ns*sizeof(float));
        fftwf_execute(f2t);

        memcpy(taup[i], sig_t, ns*sizeof(float));
    }

    // Form the filter with the P-T domain Radon Spectrum
    t_len=sqrtf(ns);
    t_len=MIN(t_len, (int)(4.0*ns/nf_mid));

    memset(tr_pt_filt[0], 0, np*ns*sizeof(float));
    printf("p_len=%d, t_len=%d\n", p_len, t_len);
    for(i= 0; i< np; i++)
    {
       for(j= 0; j< ns; j++)
       {
          for(k= 0; k< p_len; k++)
          {
              p1= (i+ k+ np)% np;
              p2= (i- k+ np)% np;
              for(l= 0; l< t_len; l++)
							{
                 t1=(j+ l+ ns)%ns;
                 t2=(j- l+ ns)%ns;

                 tr_pt_filt[i][j]+= taup[p1][t1]*
                                    taup[p1][t2]*
                                    taup[p2][t1]*
                                    taup[p2][t2];
               }
           }
           tr_pt_filt[i][j]=pow(tr_pt_filt[i][j],2);
        }
     }

/* Output the loop_time's filter and the taup result */
    filt_total=0;
    for(i=0; i<np; i++){
        for(j=0; j<ns; j++){
            filt_total+=tr_pt_filt[i][j];
        }
    }

    for(i=0; i<np; i++){
        for(j=0; j<ns; j++){
            tr_pt_filt[i][j]= tr_pt_filt[i][j]/filt_total;
        }
    }

/* Debug Information */
   FILE *fp_write= NULL;
	 fp_write= fopen("./tr_pt_filt.dat", "wb");
	 if(fp_write==NULL)
   printf("the file is failed to write");
   fwrite(&tr_pt_filt[0][0], sizeof(float), np*ns, fp_write);
   fclose(fp_write);

    
    filt_total= 0.05;

    filt_now= 0;
    for(i= 0; i< np; i++){
        max= 0.0;
        for(j= 0; j< ns; j++){
            max+= tr_pt_filt[i][j];
        }
        if(max> filt_total)pcycle[i]= 1;
    }

/*
   FILE *fp_write1= NULL;
	 fp_write1= fopen("./pcyclein.dat", "wb");
	 if(fp_write1==NULL)
   printf("the file is failed to write");
   fwrite(&pcycle[0], sizeof(float), np, fp_write1);
   fclose(fp_write1);
*/

/* ====================== Get the err trace ======================= */
    free1float(cor);
    free1float(coe);
    free1float(sig_t);
    free1floatcomplex(sig_f);
    
    free2float(tr_in);
    free2float(filt_max);
    free2float(tr_pt_filt);
    free2floatcomplex(data_fx);
    free2floatcomplex(data_fp);

//    fftwf_destroy_plan(t2f);
//    fftwf_destroy_plan(f2t);
    return 0;
}
