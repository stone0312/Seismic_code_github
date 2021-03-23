#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi 3.1415926
main()
{
    int i,j;
    int Xnumber=512;
    int dx=4;
    int Tnumber=1024;
    float dt=0.002;
    int f0=35;
    float *trace;
    trace=(float*)malloc(sizeof(float)*Xnumber*Tnumber);
    for(i=0;i<Xnumber;i++)
         {

            for(j=0;j<Tnumber;j++)
                {
                     if(i==255)
                     trace[i*1024+j]=(1-2*(pi*f0*(j-150)*dt)*(pi*f0*(j-150)*dt))*exp(-(pi*f0*(j-150)*dt)*(pi*f0*(j-150)*dt));
                     else
                     trace[i*1024+j]=0;

        }
    }





    FILE *fp;
    fp=fopen("./zero_offset_data.dat","wb");
    fwrite(trace,sizeof(float),512*1024,fp);
    fclose(fp);
    free(trace);

}


