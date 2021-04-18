/**
  * filename   :  Acoustic Forward Program
  * CopyRight  (c) 2019-
  *
  * name       :  ShengShen TJU
  * date       :  2019-11-29
  * describe   :  FD with regular grid and sponge boundary & 2nd in time & 16th in space
  **/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.1415926

int main()
{
    //*************************************************************************//
    //                          define par                                     //
    //*************************************************************************//
    float **P1,**P2,**Vp,**DEN,**Vx1,**Vx2,**Vz1,**Vz2,**record,**R;/*{{{*/
    float dt, dx, dz;
    int   i, ii, j, jj, k, Lw, NX, NZ, PML, Nx,Nz,Nt,Snapshot;
    float C1,C2,C3,C4,C5,C6,C7,C8;
    float LS1,LS2;/*}}}*/
    //*************************************************************************//
    //                          assignment par                                 //
    //*************************************************************************//
    NX  = 100; //{{{
    NZ  = 100; //
    PML = 20;  //
    Nx  = NX+2*PML;//X方向加上吸收边界后的大小
    Nz  = NZ+2*PML;//Z方向加上吸收边界后的大小
    Nt  = 2000;    //正演时间
    dx  = 10.0;    //x方向采样
    dz  = 10.0;    //z方向采样
    dt  = 0.0005;  //时间采样
    //有限差分系数
    C1 = 1.2340911;
    C2 = -1.0664985e-1;
    C3 = 2.3036367e-2;
    C4 = -5.3423856e-3;
    C5 = 1.0772712e-3;
    C6 = -1.6641888e-4;
    C7 = 1.7021711e-5;
    C8 = -8.5234642e-7;

    float absorb[NX+2*PML][NZ+2*PML];/*}}}*/
    //*************************************************************************//
    //                          open up the arraies                            //
    //*************************************************************************//
    P1     = new float* [Nx];/*{{{*/
    P2     = new float* [Nx];
    Vp     = new float* [Nx];
    Vx1    = new float* [Nx];
    Vz1    = new float* [Nx];
    Vx2    = new float* [Nx];
    Vz2    = new float* [Nx];
    DEN    = new float* [Nx];
    record = new float* [Nx];
    R      = new float* [Nx];
    for(i = 0;i<Nx;i++)
    {
        P1[i]     = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        P2[i]     = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        Vp[i]     = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        Vx1[i]    = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        Vz1[i]    = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        Vx2[i]    = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        Vz2[i]    = new float [Nz];
    }
    for(i = 0;i<Nx;i++)
    {
        DEN[i]    = new float [Nz];
    }
    for(i = 0;i<NX;i++)
    {
        record[i] = new float [Nt];
    }
    for(i = 0;i<Nx;i++)
    {
        R[i]      = new float [Nz];
    }
/*}}}*/
    //*************************************************************************//
    //                          initializing par                               //
    //*************************************************************************//
    for(i = 0; i < Nx; i++)/*{{{*/
    {
        for(j = 0; j < Nz;j++)
        {
            P1[i][j]  = 0.0;
            P2[i][j]  = 0.0;
            Vx1[i][j] = 0.0;
            Vx2[i][j] = 0.0;
            Vz1[i][j] = 0.0;
            Vz2[i][j] = 0.0;
            R[i][j]   = 0.0;
        }
    }
    for(i = 0; i < NX; i++)
    {
        for(j = 0; j < Nt;j++)
        {
            record[i][j] = 0.0;
        }
    }

    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Nz;j++)
        {
            absorb[i][j]=1.0;//吸收衰减函数初始值为1.0
        }
    }

    //*************************************************************************//
    //                          速度模型设计                                   //
    //*************************************************************************//
    for(i=0; i<Nx; i++)
        for(j=0; j<Nz;j++)
        {
            Vp[i][j]=2000.0;
            DEN[i][j]=2400.0;
        }

    //速度模型保存
    FILE *fp_model;
    if((fp_model = fopen ("../file/Model.dat", "wb"))!=NULL)
    {

        for (i=PML;i<Nx-PML;i++)
        {
            for (j=PML;j<Nz-PML;j++)
            {
                fwrite (&Vp[i][j] , sizeof(float), 1, fp_model);

            }
        }
        fclose (fp_model);
    }
/*}}}*/
    //*************************************************************************//
    //                          design_vel_model                               //
    //*************************************************************************//
    //                          design_sponge_boundary                         //
    for(i=PML;i<Nx-PML;i++)/*{{{*/
    {
        for(j=0;j<PML;j++)
        {
           absorb[i][j]=-0.38*(PML-j-1)*(PML-j-1)/(1.0*PML*PML)+1;
        }
    }
    //右边界//
    for(i=PML;i<Nx-PML;i++)
    {
        for(j=Nz-PML;j<Nz;j++)
        {
            absorb[i][j]=-0.38*(j-Nz+PML)*(j-Nz+PML)/(1.0*PML*PML)+1;
        }
    }
    //上边界//
    for(i=0;i<PML;i++)
    {
        for(j=0;j<Nz;j++)
        {
            absorb[i][j]=-0.38*(PML-i-1)*(PML-i-1)/(1.0*PML*PML)+1;
        }
    }
    //下边界//
    for(i=Nx-PML;i<Nx;i++)
    {
        for(j=0;j<Nx;j++)
        {
            absorb[i][j]=-0.38*(i-Nx+PML)*(i-Nx+PML)/(1.0*PML*PML)+1;
        }
    }
/*}}}*/
    //*************************************************************************//
    //                          design_sourc                                   //
    //*************************************************************************//
    float signal[Nt];/*{{{*/
    float freq=28;
    //雷克子波，freq为频率,signal为雷克子波
    float t, t1;
    for(i=0;i<Nt;i++)
    {
        t=dt*i;
        t1=1.0/freq;//双边雷克子波
        signal[i]=50*(1-2*PI*PI*freq*freq*(t-t1)*(t-t1))
                    *exp(-PI*PI*freq*freq*(t-t1)*(t-t1));
    }
/*}}}*/
    //*************************************************************************//
    //                          location of sourc                              //
    //*************************************************************************//
    int nx_location,nz_location;/*{{{*/
    nx_location=Nx/2;
    nz_location=Nz/2;
/*}}}*/
    //*************************************************************************//
    //                          remove files                                   //
    //*************************************************************************//
    remove("../file/Vxall.dat");/*{{{*/
    remove("../file/Vzall.dat");
    remove("../file/record.dat");
    remove("../file/Pall.dat");
/*}}}*/
    //*************************************************************************//
    //                          progress bar                                   //
    //*************************************************************************//
    int i_timebar = 0;/*{{{*/
    char timebar[25];
    const char *timebar_lable = "|/-\\";
    timebar[0] = 0;
/*}}}*/
    //*************************************************************************//
    //                          cal_wave_recursion                             //
    //*************************************************************************//
    for(k=0; k<Nt; k++) //begain cal time of wave_recursion{{{
    {
        //进度条显示
        i_timebar=k*25/Nt;
        printf("[%-25s][%d%%][%c]\r", timebar, (i_timebar+1)*4, timebar_lable[k%4]);
        fflush(stdout);
        timebar[i_timebar] = '#';
        timebar[i_timebar+1] = 0;

        //加入震源
        if(k<Nt)
        {
            P2[nx_location][nz_location]=P2[nx_location][nz_location]
                -dt*DEN[nx_location][nz_location]*Vp[nx_location][nz_location]*signal[k];
        }
        //时间递推计算

        //计算Vx
        for(i=8; i<Nx-8; i++)
        {
            for(j=8; j<Nz-8; j++)
            {
                Vx1[i][j]=Vx2[i][j]-(1.0/DEN[i][j])*(dt/dx)*(
                        C1*(P2[i+1][j]-P2[i][j])+
                        C2*(P2[i+2][j]-P2[i-1][j])+
                        C3*(P2[i+3][j]-P2[i-2][j])+
                        C4*(P2[i+4][j]-P2[i-3][j])+
                        C5*(P2[i+5][j]-P2[i-4][j])+
                        C6*(P2[i+6][j]-P2[i-5][j])+
                        C7*(P2[i+7][j]-P2[i-6][j])+
                        C8*(P2[i+8][j]-P2[i-7][j]));
            }
        }

        //计算Vz
        for(i=8; i<Nx-8; i++)
        {
            for(j=8; j<Nz-8; j++)
            {
                Vz1[i][j]=Vz2[i][j]-(1.0/DEN[i][j])*(dt/dz)*(
                        C1*(P2[i][j+1]-P2[i][j])+
                        C2*(P2[i][j+2]-P2[i][j-1])+
                        C3*(P2[i][j+3]-P2[i][j-2])+
                        C4*(P2[i][j+4]-P2[i][j-3])+
                        C5*(P2[i][j+5]-P2[i][j-4])+
                        C6*(P2[i][j+6]-P2[i][j-5])+
                        C7*(P2[i][j+7]-P2[i][j-6])+
                        C8*(P2[i][j+8]-P2[i][j-7]));
            }
        }

        //计算波场值P
        for(i=8; i<Nx-8; i++)
        {
            for(j=8;j<Nz-8;j++)
            {
                P1[i][j]=P2[i][j]-DEN[i][j]*Vp[i][j]*Vp[i][j]*dt*((1.0/dx)*(
                            C1*(Vx1[i][j]-Vx1[i-1][j])+
                            C2*(Vx1[i+1][j]-Vx1[i-2][j])+
                            C3*(Vx1[i+2][j]-Vx1[i-3][j])+
                            C4*(Vx1[i+3][j]-Vx1[i-4][j])+
                            C5*(Vx1[i+4][j]-Vx1[i-5][j])+
                            C6*(Vx1[i+5][j]-Vx1[i-6][j])+
                            C7*(Vx1[i+6][j]-Vx1[i-7][j])+
                            C8*(Vx1[i+7][j]-Vx1[i-8][j]))
                        +(1.0/dz)*
                        (
                         C1*(Vz1[i][j]-Vz1[i][j-1])+
                         C2*(Vz1[i][j+1]-Vz1[i][j-2])+
                         C3*(Vz1[i][j+2]-Vz1[i][j-3])+
                         C4*(Vz1[i][j+3]-Vz1[i][j-4])+
                         C5*(Vz1[i][j+4]-Vz1[i][j-5])+
                         C6*(Vz1[i][j+5]-Vz1[i][j-6])+
                         C7*(Vz1[i][j+6]-Vz1[i][j-7])+
                         C8*(Vz1[i][j+7]-Vz1[i][j-8])));

            }
        }

        //时间片波场值替换
        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Nz; j++)
            {
                P2[i][j]  = P1[i][j];
                Vx2[i][j] = Vx1[i][j];
                Vz2[i][j] = Vz1[i][j];
            }
        }

        //应用吸收边界条件
        //左边界//
        for(i = 0; i < Nx; i++)
        {
            for(j = 0; j < Nz; j++)
            {
                Vx1[i][j] = Vx1[i][j]*absorb[i][j];
                Vx2[i][j] = Vx2[i][j]*absorb[i][j];
                Vz1[i][j] = Vz1[i][j]*absorb[i][j];
                Vz2[i][j] = Vz2[i][j]*absorb[i][j];
                P2[i][j]  = P2[i][j]*absorb[i][j];
                P1[i][j]  = P1[i][j]*absorb[i][j];
            }
        }

        //输出Snapshot时刻的时间切片Vx
        FILE *fp_snap;
        Snapshot=500;  //选择输出的时间切片
        if(k==Snapshot)
        {
            fp_snap=fopen("../file/Snap.dat","wb");
            for(ii=PML;ii<Nx-PML;ii++)
                for(jj=PML;jj<Nz-PML;jj++)
                {
                    LS1=float(Vx1[ii][jj]);
                    fwrite(&LS1,sizeof(float),1,fp_snap);
                }
            fclose(fp_snap);
        }

        //输出所有时间的P
        if(k%10==0)
        {
            FILE *fpp;
            if((fpp = fopen ("../file/Pall.dat", "a+"))!=NULL)
            {
                for(i=PML;i<Nx-PML;i++)
                {
                    for(j=PML;j<Nz-PML;j++)
                    {
                        fwrite (&P1[i][j] , sizeof(float), 1, fpp);

                    }
                }
            }
            fclose (fpp);
        }

        //地表地震记录剖面record
        for(i=PML; i<Nx-PML; i++)
            record[i-PML][k]=P1[nx_location][i];


    }
    //波场时间递推计算结束


    //out_record
    FILE *fp_record;
    fp_record=fopen("../file/record.dat", "wb");
    for(i=0; i<NX; i++)
        for(k=0; k<Nt; k++)
        {
            LS2=float(record[i][k]);
            fwrite(&LS2,sizeof(float),1, fp_record);
        }
    fclose(fp_record);

    printf("\nEnd of Calculation!!\n");
/*}}}*/
    //*************************************************************************//
    //                          free memory                                    //
    //*************************************************************************//
    for(i=0;i<Nx;i++)/*{{{*/
        delete []	P1[i];
    delete []P1;
    for(i=0;i<Nx;i++)
        delete []	P2[i];
    delete []P2;
    for(i=0;i<Nx;i++)
        delete []	Vp[i];
    delete []Vp;
    for(i=0;i<Nx;i++)
        delete []	Vx1[i];
    delete []Vx1;
    for(i=0;i<Nx;i++)
        delete []	Vz1[i];
    delete []Vz1;
    for(i=0;i<Nx;i++)
        delete []	Vx2[i];
    delete []Vx2;
    for(i=0;i<Nx;i++)
        delete []	Vz2[i];
    delete []Vz2;
    for(i=0;i<Nx;i++)
        delete []	DEN[i];
    delete []DEN;
    for(i=0;i<NX;i++)
        delete []	record[i];
    delete []record;
    for(i=0;i<Nx;i++)
        delete []	R[i];
    delete []R;


    return 0;/*}}}*/
}
