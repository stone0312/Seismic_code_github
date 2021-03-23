//------------------2020.5.12--------------//
//--------------phase-shift----------------//
//-------------------sy--------------------//

#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include"fftw3.h"

void syps(float **p, float**o, int nx, int nz, int nt, float dx, float dt, float v);

int main()
{
    int nx,nz,nt;
    int i;
    float v;
    float dx,dt;
    float **p,**o;
    FILE *fi,*fo;

    //input par
    nx=50;
    nz=250;
    nt=250;
    dx=8;
    dt=0.001;
    v=3000;

    //opeaning memory
    p=(float **)malloc(sizeof(float *)*nx);//put
   
    for(i=0;i<=nx;i++)
    {
         p[i]=(float *)malloc(sizeof(float )*nt);
    }
    o=(float **)malloc(sizeof(float *)*nx);//out
   
    for(i=0;i<nx;i++)
    {
         o[i]=(float *)malloc(sizeof(float )*nz);
    }

    //read
    fi=fopen("/home/shiyu/stolt/zero_offset_data3.dat","rb");
    for(i=0;i<nx;i++)
    {
    fread(p[i],sizeof(float),nt,fi);
    }
    syps(p,o,nx,nz,nt,dx,dt,v);
    fo=fopen("phase-shift_out2.dat","wb");
    for(i=0;i<nx;i++)
    {
    fwrite(o[i],sizeof(float),nz,fo);
    }
    //free memory
    for(i=0;i<nx;i++)
    {
    free(p[i]);
    }
    for(i=0;i<nx;i++)
    {
    free(o[i]);
    }
    free(p);
    free(o);

    fclose(fi);
    fclose(fo);

    return 0;
}


void syps(float **p, float**o, int nx, int nz, int nt, float dx, float dt, float v)
{
    
    FILE *kk;
    int i,j,k;
    float **q;
    float dz,kz;
    float a,b;
    float kx,omg;
    int Nt;

    fftw_complex *in,*out;
    fftw_plan pn;

    dz=5;
    //dz=dt*v*2.5;
    Nt=80;
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *nt *nx);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *nt *nx);

    q=(float **)malloc(sizeof(float *)*nx);
    for(i=0;i<nx;i++)
    {
    q[i]=(float *)malloc(sizeof(float )*nt);
    }
    //2d-fft
    for(i=0;i<nx;i++)
    {
    for(j=0;j<nt;j++)
    {
    in[i*nt+j][0]=p[i][j];
    in[i*nt+j][1]=0.0;
    }
    }

    pn=fftw_plan_dft_2d(nx, nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pn);

    for(i=0;i<nx;i++)
    {
    for(j=0;j<nt;j++)
    {
        p[i][j]=out[i*nt+j][0];
        q[i][j]=out[i*nt+j][1];

    }
    }
    for(i=0;i<nx;i++)
    {
    for(j=0;j<nz;j++)
    {
        o[i][j]=0.0;
    }
    }



    for(k=0;k<nz;k++)
    {
    for(i=0;i<Nt;i++)
    {
    for(j=0;j<nx;j++)
    {
    in[j][0]=p[j][i];
    in[j][0]=q[j][i];
    }
    pn=fftw_plan_dft_1d(nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pn);
    for(j=0;j<nx;j++)
    {
        o[j][k]+=out[j][0];
    }

    }

    //cal_nest_p
    if(k==nz-1)
    {
    break;
    }

    for(i=1;i<Nt;i++)
    {
    for(j=0;j<nx;j++)
    {
    omg=i/(nt*dt);
    if(j<(nx/2+1))
        kx=j/(nx*dx);
    else
        kx=(nx-j)/(nx*dx);
    if((4*omg*omg/v/v-kx*kx)<0)
        kz=0.0;
    else
        kz=sqrt(4*omg*omg/v/v-kx*kx);
    a=p[j][i];
    b=q[j][i];
    p[j][i]=a*cos(kz*dz)+b*sin(kz*dz);
    q[j][i]=-a*sin(kz*dz)+b*cos(kz*dz);
    }
    }
    }
    for(i=0;i<nx;i++)
    {
        free(q[i]);
    }
    free(q);
//    fftw_destory_plan(pn);
    fftw_free(in);
    fftw_free(out);
}

