//-------------phase-shift------------//
//-------------2020.06.09-------------//
//-------------by sy------------------//


#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <fftw3.h>
#define pi 3.1415926

float   main()
{
    int i,j,k;
    int Xn = 512;
    int Tn = 1024;
    int Zn = 1024;
    int dx = 4;
    int dz = 4;
    float dt = 0.002;
    float v = 2000;

    float *trace_in,*result_out,*temp_real,*temp_ima;
    trace_in = (float*)malloc(sizeof(float)*Xn*Tn);
    result_out = (float*)malloc(sizeof(float)*Xn*Zn);
    temp_real = (float*)malloc(sizeof(float)*Xn*Zn);
    temp_ima = (float*)malloc(sizeof(float)*Xn*Zn);

    FILE *fi,*fo;
    fi = fopen("/home/shiyu/stolt/zero_offset_data1.dat","rb");
    fread(trace_in,sizeof(float),Xn*Tn,fi);
    fclose(fi);

//------------------------design-fft-plan--------------------------------//

    fftwf_complex *in,*out;
    in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Xn*Tn);
    out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Xn*Tn);
    fftwf_plan p1,p2;

    for(i=0;i<Xn;i++)
    {
        for(j=0;j<Tn;j++)
        {
            in[i*Tn+j][0]=trace_in[i*Tn+j];
            in[i*Tn+j][1]=0.0;
        }
    }

    p1 = fftwf_plan_dft_2d(Xn,Tn,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftwf_execute(p1);

//-----------------------process-the-exp---------------------------------//
   for(i=0;i<Xn;i++)
   {
        for(j=0;j<Tn;j++)
        {
            temp_real[i][j] = result_out[i*Tn+j][0];
            temp_ima[i][j] = result_out[i*Tn+j][1];
        }
   }
    



























}
