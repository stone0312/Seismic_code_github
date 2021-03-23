//---------------fftw3_test-----------------//
//---------------2020.5.1-------------------//
//---------------by sy----------------------//

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#define PI 3.1415926

int main()
{
    int len=128;
//    int len=256;
    int i;
    double  *in=NULL;
    double  *amp;
    amp=(double*)malloc(sizeof(double)*len);
    fftw_complex *out=NULL;
    fftw_plan p;
    in= (double *)fftw_malloc(sizeof(double) * len);
    out= (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * len);
    double dx=1.0/len;
//    double dx=2.0/len;

    //enter real number
    for (i=0;i<len;i++)
    {
        in[i]=0.5*sin(2*PI*20*dx*i)+2*sin(2*PI*40*dx*i);//+sin(4*PI*dx*i);
    }

    //fft
    p = fftw_plan_dft_r2c_1d(len, in, out, FFTW_ESTIMATE);
    fftw_execute(p);
    //Amplitude
    for (i=0;i<len;i++)
    {
         amp[i] = sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);//real_amp=amp/N/2
    }

   for(i=0;i<len;i++)
    {
         printf("%.2f\n",amp[i]);
    }
    FILE *fp;
    fp=fopen("./Amplitude_test_contine","wd");
    fwrite(amp,sizeof(double),128,fp);
    fclose(fp);

    printf("\n");

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;
}

