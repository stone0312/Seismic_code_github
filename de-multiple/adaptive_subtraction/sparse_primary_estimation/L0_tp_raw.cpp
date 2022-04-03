using namespace std;
#include <iostream>
#include <fstream>
#include <string.h>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"

#ifndef PAI
#define PAI (3.141592653589793)
#endif

#ifndef MAX
#define MAX(a,b) (a)>(b)?(a):(b)
#endif

#ifndef MIN
#define MIN(a,b) (a)>(b)?(b):(a)
#endif

float norm( complex<float> a)
{
        float num;

         num= a.real()*a.real()+ a.imag()* a.imag();

        return num;
}

float minn(float *a, int n)
{
        int i;
        float min;

        min= a[0];

        for(i=1; i< n; i++)
        {
                if(a[i]< min)min=a[i];
        }

        return min;
}

float maxn(float *a, int n)
{
        int i;
        float max;

        max= a[0];

        for(i=1; i< n; i++)
        {
                if(a[i]> max)max=a[i];
        }

        return max;
}
int sis_RT_fre_2D(int nf, float df, int nx, float *cor_in, complex<float> **datainfx,
                          int np, float *p, complex<float> **dataoutfp)
{
    int i, j, k;
    float fre;
    float cor_min, cor_mid, cor_max;
    float *cor;
    complex<float> phase;
    complex<float> **dataout_fp;

    cor= alloc1float(nx);
    dataout_fp= alloc2complex (nf, np);

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

//    cout<<"1"<<endl; 

    for(i= 0; i< nf; i++)
    {
       fre= 2* PAI* df* i;
       for(j= 0; j< np; j++)
       {
          dataout_fp[j][i]= 0;
          for(k= 0; k< nx; k++)
          {
             phase.real()= cos(fre* p[j]* cor[k]);
             phase.imag()=sin(fre* p[j]* cor[k]);
             dataout_fp[j][i]+= phase* datainfx[k][i];
          }
         dataoutfp[j][i]= dataout_fp[j][i]/ (float)nx;
       }
    }

//   cout<<"2"<<endl; 

    free1float(cor);
    free2complex(dataout_fp);

    return 1;
}

int sis_RT_HR_Time_2D(int ns, int nx, float dt, float *cor_in, float **data,
    int np, float *p, float **taup,float **tp_filter , float **tr_pt)
{

    int i, j, k, l;
    double tmp;
    float df, nf_max;
    float cor_min, cor_max, cor_mid;
    int **tr_pt_filt_record;
    float *cor, *coe, *sig_t;
    float **tr_pt_hr, **tr_in, **tr_out, **tr_pt_get, **filt_max;
    double **tr_pt_filt;
    complex<float> *sig_f, **data_fx, **data_fp;

    //use 2-D array to access input data and out data
    df= 1.0/ns;
    nf_max=0.5*ns + 1;

    //alloc memory
    cor= alloc1float(nx);
    coe= alloc1float(ns);
    sig_t= alloc1float(ns);
    sig_f= alloc1complex(nf_max);

    tr_pt_hr= alloc2float(ns, np);
    tr_pt_filt_record= alloc2int(ns, np);
    tr_in= alloc2float(ns, nx);
    tr_out= alloc2float(ns, nx);
    tr_pt_get= alloc2float(ns, np);
    filt_max= alloc2float(ns, np);
    tr_pt_filt= alloc2double(ns, np);

    data_fx= alloc2complex(nf_max, nx);
    data_fp= alloc2complex(nf_max, np);

    for(i=0; i< nx; i++)
    {
       for(j= 0; j< ns; j++)
       {
          tr_in[i][j]= data[i][j];
       }
    }

    for(i=0; i< np; i++)
    {
       for(j= 0; j< ns; j++)
       {
          tr_pt_get[i][j]= taup[i][j];
       }
    }

    //renew the coordinate of input data
    cor_min= minn(cor_in, nx);
    cor_max= maxn(cor_in, nx);
    cor_mid= 0.5*(cor_max + cor_min);
//    printf("cor_min= %f, cor_max= %f, cor_mid=%f\n", cor_min, cor_max, cor_mid);

    for(i= 0; i< nx; i++)
    {
        cor[i]=cor_in[i] - cor_mid;
    }

    for(i= 0; i< np; i++)
    {
        for(j= 0; j< ns; j++)
        {
            tr_pt_filt_record[i][j]= 1;
        }
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
    memset(&data_fx[0][0], 0, nf_max* nx* sizeof(complex<float>));

    for(i= 0; i< nx; i++)
    {
        for(j= 0; j< ns; j++)
        sig_t[j] = tr_in[i][j]/ns;

        memset(sig_f, 0, nf_max*sizeof(complex<float>));
        fftwf_execute(t2f);

        for(j= 0; j< nf_max; j++)
	{
            data_fx[i][j]=sig_f[j];
        }
    }

/* get the information of data in frequency domain */
    double ener_total=0;
    for(i= 0; i< nf_max; i++)
    {
        for(j= 0; j< nx; j++)
             ener_total+= norm(data_fx[j][i]);
    }
    float const norm_coe= sqrtf(ener_total);

/* normalizate the energy to 1 */
    for(i= 0; i< nf_max; i++)
    {
       for(j= 0; j< nx; j++)
       {
          data_fx[j][i]= data_fx[j][i]/ norm_coe;
       }
    }


    cout<<"3333333333"<<endl;


    double const ener_out= (1- 1.0e-8);
    tmp= 0.5;
    ener_total= 0;

/* get the half part energy and almost total energy num */



    for(i= 0; i< nf_max; i++)
    {
        for(j= 0; j< nx; j++)ener_total+= norm(data_fx[j][i]);
        if(ener_total>= tmp)break;
    }

    int const nf_mid= i;

    for(; i< nf_max; i++)
    {
       for(j= 0; j< nx; j++)
         ener_total+= norm(data_fx[j][i]);
       if(ener_total>= ener_out)break;
    }
    int const nf= i;

//  printf("nf_mid=%d,nf= %d\n", nf_mid, nf);

/*  begin to get the Radon spectrum of the input data */

    float energy_aver;
    double err[5];
    double err_pre= 0;
    double err_now;
    double filt_total= 0;
    double filt_pre= 0, filt_now= 0;
//	int nf_stat[10];
//	float nf_ave;
    double const ener_end= MAX(1.0-5*ener_out, 1.0e-30);
    double ener_now= ener_out;
    float const dp= fabsf(p[1]- p[0]);
    i= 4.0/ (nf* df* fabsf(cor_max- cor_min)* dp);
    i= MAX(i, 5);
    int const p_len= MIN(i, (int)sqrtf(np));
    int const agc_p_win= sqrtf(np);
    int const agc_t_win= sqrtf(ns);
    int const t_len_min= 3;
    int t_len;
    int p1, p2;
    int t1, t2;
    int t_index= 0, p_index= 0;
    int nf_now;
    float maxnum;
    float fre;
    complex<float> phase;

    memset(tr_pt_get[0], 0, np* ns* sizeof(float));
    memset(tr_out[0], 0, nx* ns* sizeof(float));

    FILE *fout1=fopen("filter", "wb+");

    int loop_time=0;
    while(ener_now>= ener_end)
     {

      loop_time++;

      // Change the data to F-P domain
      sis_RT_fre_2D(nf_max, df, nx, cor, data_fx, np, p, data_fp);

      // Change the data to P-T domain
      for(i= 0; i< np; i++)
      {
        memset(sig_f, 0, nf_max* sizeof(complex<float>));
        for(j= 0; j< nf_max; j++)sig_f[j]= data_fp[i][j];

        memset(sig_t, 0, ns*sizeof(float));
        fftwf_execute(f2t);

        memcpy(tr_pt_hr[i], sig_t, ns*sizeof(float));
      }
      if(loop_time==1)
      for(i= 0; i< np; i++)
      {
        for(j= 0; j< ns; j++)
        {
          tr_pt[i][j]= tr_pt_hr[i][j];
        }
      }
      fwrite(&tr_pt[0][0], sizeof(float), ns*np, fout1);
     // break;

        // Form the filter with the P-T domain Radon Spectrum
        t_len=sqrtf(ns);
        t_len=MIN(t_len, (int)(4.0*ns/nf_mid));
        memset(tr_pt_filt[0], 0, np*ns*sizeof(double));
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
                        tr_pt_filt[i][j]+= (double)tr_pt_hr[p1][t1]*
                                          (double)tr_pt_hr[p1][t2]*
                                          (double)tr_pt_hr[p2][t1]*
                                          (double)tr_pt_hr[p2][t2];
                    }
                }

                tr_pt_filt[i][j]=pow(tr_pt_filt[i][j],2)*tr_pt_filt_record[i][j];
            }
        }

/* Output the loop_time's filter and the taup result */
        if(loop_time == 1)
	{
/*
          FILE *fp_wr= NULL;
          fp_wr=fopen("/data1/swq/sparse_primary/raw_taup.dat", "wb");
          if(fp_wr==NULL)
          printf("the file is failed to open");
*/
//          fwrite(&tr_pt[0][0], sizeof(float), np*ns, fp_wr);
//          fclose(fp_wr);

          for(i= 0; i< np; i++)
          {
             for(j= 0; j< ns; j++)
             {
                tp_filter[i][j]= (float)tr_pt_filt[i][j];
//              fwrite(&maxnum, sizeof(float), 1, fout1);
             }
          }
//          fclose(fout1);
        }

/* Get the tr_pt_filter's max energy and the total energy of first loop_time ???*/

        if(loop_time== 1)
	 {
            for(i= 0; i< np; i++)
	     {
                for(j= 0; j< ns; j++)
		 {
                    filt_max[i][j]= -1;
                    for(k= -agc_p_win; k<= agc_p_win; k++)
		     {
                        p1= (i+ k);
                        if(p1< 0|| p1>= np)continue;
                        for(l= -agc_t_win; l<= agc_t_win; l++)
			{
                            t1= (j+ l);
                            if(t1< 0||t1>= ns)continue;
                            if(filt_max[i][j]< tr_pt_filt[p1][t1])
                               filt_max[i][j]= tr_pt_filt[p1][t1];
                        }
                    }
                }
             }

            filt_total= 0;
            for(i= 0; i< np; i++)
	     {
                for(j= 0; j< ns; j++)
		{
                   filt_total+= tr_pt_filt[i][j];
                }
            }
        // Define the energy end boundary 
            filt_pre= 1.1* filt_total/ (np* ns);
            filt_total= 1.0e-15* filt_total;

        }

        // Get the total loop_time's energy of filter
        filt_now= 0;
        for(i= 0; i< np; i++)
	{
            for(j= 0; j< ns; j++)
	    {
                filt_now+= tr_pt_filt[i][j];
            }
        }

        //judge weather the energy of filter is small enough to break
        if(filt_now< filt_total)break;

        //if the energy of filter is growing, break
        filt_now= filt_now/ (np*ns);
        if(filt_now> 1.1* filt_pre)break;
        
				//!!! Here change the filt_pre !!!
        filt_pre= filt_now;

        energy_aver= filt_now;

        for(i= 0; i< np; i++)
	{
            for(j= 0; j< ns; j++)
	    {

                tr_pt_filt[i][j]= pow(tr_pt_filt[i][j], 2)/
                    ((energy_aver+ filt_max[i][j]));

//				if(tr_pt_filt[i][j]<=energy_aver)tr_pt_filt[i][j]=0;
//				else{
//					tr_pt_filt[i][j]=powf(tr_pt_filt[i][j],2)/
//									 ((energy_aver + filt_max[i][j]));
//				}

            }
        }


        // Find the P-T component with max energy
        maxnum= -1;
        filt_now= 0;
        for(i= 0; i< np; i++)
	{
            for(j= 0; j< ns; j++)
	    {
                filt_now+= tr_pt_filt[i][j];

                if(maxnum< (tr_pt_filt[i][j]))
	        {
                    maxnum= tr_pt_filt[i][j];
                    p_index= i;
                    t_index= j;
                }
            }
        }
        if(maxnum== 0)break;

        // get the selected wavelet and its spectrum
        for(i= 0; i< ns; i++)
	{
            t1= (i+ ns- t_index)%ns;
            sig_t[i]= coe[t1]*tr_pt_hr[p_index][i]/ns;
        }
        memset(sig_f, 0, nf_max*sizeof(complex<float>));
        fftwf_execute(t2f);

        // get the main frequency of this wavelet
        tmp= 0;
        for(i= 0; i< nf; i++)
	{
           sig_t[i]= norm(sig_f[i]);
           tmp+= sig_t[i];
        }

        tmp= 0.6* tmp;
        for(i= 1; i< nf; i++)
	{
           sig_t[0]+= sig_t[i];
           if(sig_t[0]> tmp)break;
        }
        nf_now= i;

        t_len= 2.5*(float)ns/(float)nf_now;

//      printf("loop time:%-4d\tt_len:%-5d\tp:%-10f\tt:%-5d\t\n", loop_time, t_len, p[p_index], t_index);

        tr_pt_filt_record[p_index][t_index]= 0;

        memset(sig_t, 0, ns*sizeof(float));
        for(i= t_index-t_len, j= -t_len; i<= t_index+ t_len; i++, j++)
        {
            t1= (i + ns)%ns;
            tmp= exp(-4.0*pow(j/t_len, 2));

            //add the wavelet to result only and if only t_len>t_len_min
            if(t_len> t_len_min)tr_pt_get[p_index][t1]+= tr_pt_hr[p_index][t1]* tmp;
           //tr_pt_get[p_index][t1]+=tr_pt_hr[p_index][t1]*tmp;

            //record the signal to subtract
            sig_t[t1]=tr_pt_hr[p_index][t1]*tmp;
        }

        // renew the data in F-x domain
        memset(sig_f, 0, nf_max*sizeof(complex<float>));
        for(i=0; i<ns; i++)sig_t[i] = sig_t[i]/ns;
        fftwf_execute(t2f);

        for(i= 0; i< nf; i++)
	{
            fre= 2*PAI*df*i;
            for(j= 0; j< nx; j++)
	    {
                phase.real()= cos(-fre*p[p_index]*cor[j]);
                phase.imag()=sin(-fre*p[p_index]*cor[j]);
                //data_fx[j][i]= data_fx[j][i] - phase*sig_f[i];
            }
        }

        ener_now= 0;
        for(i= 0; i< nf; i++)
	{
            for(j= 0; j< nx; j++)
               ener_now+= norm(data_fx[j][i]);
        }
        // printf("%-10f\n", 100*ener_now/ener_total);

        if(loop_time==1)
				{
            for(i= 0; i< 5; i++)err[i]= 1.5*ener_now;
            err_pre= 1.5* ener_now;
        }
        else err[loop_time% 5]= ener_now;

        err_now=0;
        for(i= 0; i< 5; i++)err_now+= err[i];
        err_now= err_now* 0.2;
//  		if((err_now > (1.0 - 1.0e-3)*err_pre))break;
//      if((err_now > err_pre && loop_time>20))break;
        if((err_now> 0.995* err_pre && loop_time> 10))break;
        err_pre= err_now;

    }
    printf("loop_time=%d\n", loop_time);

    for(i=0; i< np; i++)
    {
			for(j= 0; j< ns; j++)
			{
        taup[i][j]= tr_pt_get[i][j]* norm_coe;
		  }
    }

/* ====================== Get the err trace ======================= */
    free1float(cor);
    free1float(coe);
    free1float(sig_t);
    free1complex(sig_f);
    
    free2int(tr_pt_filt_record);
    free2float(tr_in);
    free2float(tr_out);
    free2float(filt_max);
    free2float(tr_pt_get);
    free2double(tr_pt_filt);
    free2complex(data_fx);
    free2complex(data_fp);

//    fftwf_destroy_plan(t2f);
//    fftwf_destroy_plan(f2t);
    return 0;
}

int main()
{
   int is,ir,ix,ip,it,iwin;
   int nshot,nx, np, ns,nx_lh, nx1,nwin,ns1;
   float pmax, dx, dp,dt;
   float *p, *cor_in,**tp_filter,**tp_filter1;
   float **data1, **data, **data2,  **taup,**taup1,**tp_raw,**tp_raw1;

    char input[256],output1[256],output2[256],output3[256];

   ifstream swq;
   swq.open("L0_tp.par");
   swq>>input>>output1>>output2>>output3>>nshot>>nx>>np>>ns>>dt>>dx>>nx_lh>>nwin;
   swq.close();

   ns1=ns/nwin;

   /* The initial value */
   pmax= sin(80.0*PAI/180.0);
   dp= 2*pmax/(np-1);

   nx1=2*nx_lh+1;

   /* allocate the memory */
   p= alloc1float(np);
   cor_in= alloc1float(nx);
   data1= alloc2float(ns, nx+2*nx_lh);
   data= alloc2float(ns, 2*nx_lh+1);
   data2= alloc2float(ns1, 2*nx_lh+1);
   taup= alloc2float(ns, np);
   tp_filter= alloc2float(ns, np);
   tp_raw= alloc2float(ns, np);
   taup1= alloc2float(ns1, np);
   tp_filter1= alloc2float(ns1, np);
   tp_raw1= alloc2float(ns1, np);

   for(ix=nx;ix<nx+2*nx_lh;ix++)
     for(it=0;it<ns;it++)
        data1[ix][it]=0.0;

  ifstream swq1;
  swq1.open(input,ios::binary);
  if(!swq1)
   {
     cout<<"Cannot Open "<<input<<endl;
     return 0;
   }

  ofstream swq2;
  swq2.open(output1,ios::binary);
  if(!swq2)
   {
     cout<<"Cannot Open "<<output1<<endl;
     return 0;
   }

  ofstream swq3;
  swq3.open(output2,ios::binary);
  if(!swq3)
   {
     cout<<"Cannot Open "<<output2<<endl;
     return 0;
   }

  ofstream swq4;
  swq4.open(output3,ios::binary);
  if(!swq4)
   {
     cout<<"Cannot Open "<<output3<<endl;
     return 0;
   }

  /* Calculation of x coordation */
  for(ix=0; ix< nx; ix++)
     cor_in[ix]=ix* dx;

/* Calculation of scanning P */
  for(ip=(int)(-(np-1)/2); ip<= (int)(np-1)/2; ip++)
     p[ip+(np-1)/2]=ip* dp;

 for(is= 0; is< nshot; is++)
 {
  cout<<is+1<<" Shot Begins ..."<<endl;   

  for(ix=nx_lh;ix<nx+nx_lh;ix++)
    for(it=0;it<ns;it++)
      swq1.read((char*)&data1[ix][it],sizeof(data1[ix][it]));
 

  for(ir= 0; ir< nx; ir++)   
   {
     for(ip=0;ip<np;ip++)
          for(it=0;it<ns1;it++)
           {
              taup1[ip][it]=0.0;
              tp_filter1[ip][it]=0.0;
              tp_raw1[ip][it]=0.0;
           } 

     for(ix= 0; ix< 2*nx_lh+ 1; ix++)
       for(it= 0; it<ns;it++)
          data[ix][it]=0.0;

    for(ix=0;ix<2*nx_lh+1;ix++) 
       for(it=0;it<ns;it++)
         data[ix][it]=data1[ir+ix][it];   

    for(iwin=0;iwin<nwin;iwin++) 
     {
        for(ix=0;ix<2*nx_lh+1;ix++)
          for(it=0;it<ns1;it++)
            data2[ix][it]=data[ix][iwin*ns1+it];   

         cout<<"2222"<<endl;

         /* Do the high resolution Radon Transform */
         sis_RT_HR_Time_2D(ns1, nx1, dt, cor_in, data2, np, p, taup1, tp_filter1, tp_raw1);

         cout<<"2222"<<endl;

         for(ip=0;ip<np;ip++)
          for(it=0;it<ns1;it++)
           {
              taup[ip][iwin*ns1+it]=taup1[ip][it];
              tp_filter[ip][iwin*ns1+it]=tp_filter1[ip][it];
              tp_raw[ip][iwin*ns1+it]=tp_raw1[ip][it];
           }
        
         cout<<"2222"<<endl;

     }
/* Output the Radon Transform result */

     for(ip=0;ip<np;ip++)
      for(it=0;it<ns;it++)
       { 
          swq2.write((char*)&tp_raw[ip][it],sizeof(tp_raw[ip][it]));
          swq3.write((char*)&taup[ip][it],sizeof(taup[ip][it]));
          swq4.write((char*)&tp_filter[ip][it],sizeof(tp_filter[ip][it]));
       }

     cout<<is+1<<" Shot, "<<ir+1<<" Trace Done ..."<<endl;   

   }
   cout<<is+1<<" Shot Done!"<<endl;   
  }

  cout<<" ALL DONE!"<<endl;   

  /* Release the Memory */
  free1float(p);
  free1float(cor_in);
  free2float(data);
  free2float(taup);

  return 1;
}
