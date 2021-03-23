//---------------fftw3.h_test-------------//
//---------------2020.4.20----------------//
//---------------by sy--------------------//
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"fftw3.h"
#define PI 3.1415926
int main()
{
    printf("-------------------------start----------------------\n");

    int len=8;
    int i;
    float *in=NULL;
    fftwf_complex *out=NULL;
    fftwf_plan p;
    in =(float*)fftwf_malloc(sizeof(float)*len);
    out=(fftwf_complex*)fftwf_malloc(sizeof(fftw_complex)*len);
    float dx=1.0/len;

    for(i=0;i<len;i++)
    {
        in[i]=sinf(2*PI*dx*i)+sinf(4*PI*dx*i);
        printf("%.2f",in[i]);

    }
    printf("\n\n");

    p=fftwf_plan_dft_r2c_1d(len,in,out,FFTW_ESTIMATE);
    fftwf_execute(p);

    for(i=0;i<len;i++)
    {
        float len=sqrt(out[i][0]*out[i][0]+out[i][0]*out[i][1]);
        printf("%.2f",len);
    }
    printf("\n");

    //fftwf_destory_plan(p);
    fftwf_free(in);
    fftwf_free(out);

    printf("----------------------end----------------------\n");
    return 0;
}
