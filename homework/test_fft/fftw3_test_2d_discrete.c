//---------------fftw3_test_2d_discrete--------------//
//---------------2020.06.03--------------------------//
//---------------sy----------------------------------//


#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<string.h>
#include<complex.h>

#define N 3
#define ELEM(r,c)(r*N+c)
int main()

{
    int i,j;
    fftw_complex *in, *out;
    fftw_plan p1,p2;

    in  = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) *N *N);

    out = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) *N *N);

    p1 = fftw_plan_dft_2d(N,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

    p2 = fftw_plan_dft_2d(N,N,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);


    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
            in[ELEM(i,j)][0]=0;
            in[ELEM(i,j)][1]=0;
            }
    }

            in[ELEM(1,1)][0]=1;
            in[ELEM(2,1)][0]=2;


    fftw_execute(p1);

    showresult(in,out);

    fftw_execute(p2);

    printf("BACK:\n");

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%lf %lfi\t",in[ELEM(i,j)][0]/N/N,in[ELEM(i,j)][1]/N/N);
            }
    printf("\n");
    }

    fftw_free(in);
    fftw_free(out);

    return 1;
}

int showresult(fftw_complex*in,fftw_complex*out)
{
    int i,j;

    printf("In:\n");

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%lf %lfi\t",in[ELEM(i,j)][0],in[ELEM(i,j)][1]);
            }
    printf("\n");
    }

    printf("Out:\n");

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%lf %lfi\t",out[ELEM(i,j)][0],out[ELEM(i,j)][1]);
            }
    printf("\n");
    }


return 1;

}
