#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <fftw3.h>
#define pi 3.1415926
#define N 5
float main()
{

    float sinc_interpolation(float *y0,int n,float x0);
    int i,j,k;
    int Xn=512;
    int Tn=1024;
    int Zn=1024;
    int dx=4;
    int dz=4;
    float v=2000;
    float dt=0.002;//define Fs

    float *trace;
    trace=(float*)malloc(sizeof(float)*Xn*Tn);

    FILE *fp1;
    fp1=fopen("/home/shiyu/stolt/zero_offset_data1.dat","rb");
    fread(trace,sizeof(float),Xn*Tn,fp1);
    fclose(fp1);

    fftwf_complex *in1,*out1,*in2,*out2;
    in1=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Xn*Tn);
    out1=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Xn*Tn);

    in2=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Xn*Zn);
    out2=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Xn*Zn);

    fftwf_plan p1,p2;

    for(i=0;i<Xn;i++)
    {
        for(j=0;j<Tn;j++)
        {
            in1[i*Tn+j][0]=trace[i*Tn+j];
            in1[i*Tn+j][1]=0.0;
        }
    }

    p1=fftwf_plan_dft_2d(Xn,Tn,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE);
    fftwf_execute(p1);

    float dw;
    dw=(2.0*pi)/(dt*Tn);
//  dw=1/(dt*Tn);
    float dkx;
    dkx=(2.0*pi)/(dx*Xn);
//  dkx=1/(dx*Xn);
    float dkz;
    dkz=(2.0*pi)/(dz*Zn);
//  dkz=1/(dz*Zn);



    float kx,kz,w;
    float F_nyquist = (2.0*pi)/dt/2.0;//F_nyquist=2Fs
    float h,H,X;
    int h1,h2;


    float tempt_shi[1024]={0};
    float tempt_xu[1024]={0};
    float *Num_shi,*Num_xu;



    // loop for kx(error)
    for(i=0;i<Xn;i++)//array through rows by rouws
    {

        for(k=0;k<1024;k++)
           {
             tempt_shi[k]=out1[i*Zn+k][0];
             tempt_xu[k]=out1[i*Zn+k][1];
           }
          Num_shi=tempt_shi;
          Num_xu=tempt_xu;


        // loop for kx and kz
        for(j=0;j<Zn;j++)
        {
            in2[i*Zn+j][0]= 0;
            in2[i*Zn+j][1]= 0;

            if(i<=Xn/2-1)//kz
                kx=i*dkx;
            else        //-kz
                kx=(i-Xn)*dkx;

            if(j<=Zn/2-1)//kz
                kz=j*dkz;
            else        //-kz
                kz=(j-Zn)*dkz;

            if(kz>=0)//w
            {
                w=v*sqrt(kx*kx+kz*kz);

            }
            else//-w
            {
                w=-v*sqrt(kx*kx+kz*kz);
            }
            if (fabs(w)>F_nyquist)
            {continue;}

            h=w/dw;
            H=floor(h);
            if(h==0.0)
            {
                in2[0][0]=out1[0][0];
                in2[0][1]=out1[0][1];
            }
            else if(h>0.0)
            {
                    if(H!=h)
                     {
                       // h1=floor(h);
                      //  h2=ceil(h);
                          in2[i*Zn+j][0]=sinc_interpolation(Num_shi,1024,h)*v*kz/w;
                          in2[i*Zn+j][1]=sinc_interpolation(Num_xu,1024,h)*v*kz/w;
                     }
                    if(H==h)
                     {
                        in2[i*Zn+j][0]=out1[i*Zn+(int)(h)][0]*v*kz/w;
                        in2[i*Zn+j][1]=out1[i*Zn+(int)(h)][1]*v*kz/w;

                     }

            }
            else
            {
                  if(H!=h)
                   {
                     X=h+Zn;
                    // h1=floor(h)+Zn;
                   //  h2=ceil(h)+Zn;
                    in2[i*Zn+j][0]=sinc_interpolation(Num_shi,1024,X)*v*kz/w;
                    in2[i*Zn+j][1]=sinc_interpolation(Num_xu,1024,X)*v*kz/w;


                   }

                  if(H==h)
                   {
                      in2[i*Zn+j][0]=out1[i*Zn+(int)(h+Zn)][0]*v*kz/w;
                      in2[i*Zn+j][1]=out1[i*Zn+(int)(h+Zn)][1]*v*kz/w;

                   }

           }
       }
       memset(tempt_shi,0,sizeof(float)*1024);
       memset(tempt_xu,0,sizeof(float)*1024);

    }



    p2=fftwf_plan_dft_2d(Xn,Zn,in2,out2,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(p2);

    //xiugai//
    memset(trace,0,sizeof(float)*Xn*Tn);
    for (i=0;i<Xn;i++)
    {
       for(j=0;j<Tn;j++)
       {
         trace[i*Tn+j]=out2[i*Tn+j][0]/(Xn*Tn);
       }
    }



    FILE *fp2;
    fp2=fopen("Sinc2_interpolation_f-kjieguo_normal_frequency.dat","wb");
    fwrite(trace,sizeof(float),Xn*Tn,fp2);
    fclose(fp2);




    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);

    fftwf_free(in1);
    fftwf_free(out1);
    fftwf_free(in2);
    fftwf_free(out2);

}



 float sinc_interpolation(float *y0,int n,float x0)
 {
 //    const int N=5;
    float a_indx=0;
    float front[N]={0};
    float back[N]={0};
    int i0=(int)(x0);
    int i;
    for(i=0;i<N;i++)
    {
        if(i0+1+i<n)
          back[i]=y0[i0+1+i];
        else
          back[i]=y0[n-1];
        //   back[i]=0;
        if(i0-i>=0)
         front[i]=y0[i0-i];
        else
          front[i]=y0[0];
         //   front[i]=0;

     }
    float dist=0;
    float coef=0;
    float coef_sum=0;
    float w_cos=1;
    for(i=0;i<N;i++)
    {
       w_cos=(1.0+cos((pi*i)/(N+0.0001)))*0.5;
       dist=(x0-i0)+i;
       coef=(sin(pi*dist)/(pi*dist+0.0001))*w_cos;
       a_indx=a_indx+front[i]*coef;
       coef_sum=coef_sum+coef;


       dist=(i0+1-x0)+i;
       coef=(sin(pi*dist)/(pi*dist+0.0001))*w_cos;
       a_indx=a_indx+back[i]*coef;
       coef_sum=coef_sum+coef;

    }


    return(a_indx/(float)(coef_sum+0.000001));

 }
