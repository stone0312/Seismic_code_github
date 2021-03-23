//----------------fftw3-test-discrete---------------//
//----------------2020.05.31------------------------//
//----------------sy--------------------------------//


#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<string.h>
#include<complex.h>

int main()
{
    int N=3;
    int i;
    FILE *fp, *fp1;

//    fftw_complex *in, *out;
    fftw_complex **in, **out;

    in = (fftw_complex**) fftw_malloc(sizeof(fftw_complex *) *N);

    out = (fftw_complex**) fftw_malloc(sizeof(fftw_complex *) *N);


    fftw_plan p;

    for(i=0;i<N;i++)
    {
        in[i][0]=i+1;
        in[i][1]=0;
    }

    printf("----------------input-data-------------------\n");

    for(i=0;i<N;i++)
    {
        printf("%f %fi\n",in[i][0],in[i][1]);
    }

//    fp1=fopen("./input_data_discrete.dat","wd");
//    fwrite(in,sizeof(double),128,fp1);
//    fclose(fp1);

    memset(out,0,sizeof(out));//initialization

    printf("-----------------forward---------------------\n");

    p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_execute(p);

    for(i=0;i<N;i++)
    {
        printf("%f %fi\n",out[i][0],out[i][1]);
    }


    printf("------backward-without-amp-processing--------\n");


    p = fftw_plan_dft_1d(N,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);

    fftw_execute(p);

    for(i=0;i<N;i++)
    {
        printf("%f %fi\n",in[i][0],in[i][1]);
    }

    printf("------backward-with-amp-processing----------\n");

    for(i=0;i<N;i++)
    {
        printf("%f %fi\n",in[i][0]/N,in[i][1]/N);
    }

//    fp=fopen("./backward_with_Amp_processing.dat","wd");
//    fwrite(in,sizeof(double),128,fp);
//    fclose(fp);
    fftw_destory_plan(p);
    fftw_free(in);
    fftw_free(out);
    return 0;

}

