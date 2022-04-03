#include "iostream.h"
#include "fstream.h"
#include "math.h"

int gaus( float *a,float *b,int n )
{
        int *js,l=1,k,i,j,is,p,q;
        float d,t;

        js=(int *)malloc( n*sizeof(int) );
        for( k=0;k<=n-2;k++ )
        {
                d=0.0;
                for( i=k;i<=n-1;i++ )
                for (j=k;j<=n-1;j++)
                {
                        t=fabs(a[i*n+j]);
                        if(t>d)
                        {
                                d=t;
                                js[k]=j;
                                is=i;
                        }
                }
                if(d+1.0==1.0)
                        l=0;
                else
                {
                        if (js[k]!=k)
                                for(i=0;i<=n-1;i++)
                                {
                                        p=i*n+k;
                                        q=i*n+js[k];
                                        t=a[p];
                                        a[p]=a[q];
                                        a[q]=t;

                                }

                        if(is!=k)
                        {
                                for(j=k;j<=n-1;j++)
                                {
                                        p=k*n+j;
                                        q=is*n+j;
                                        t=a[p];
                                        a[p]=a[q];
                                        a[q]=t;

                                }

                                t=b[k];
                                b[k]=b[is];
                                b[is]=t;

                        }
                }

                if(l==0)
                {
                        free(js);
                        printf("fail\n");
                        return(0);

                }
                d=a[k*n+k];
                for( j=k+1;j<=n-1;j++ )
                {
                        p=k*n+j;
                        a[p]=a[p]/d;
                }
                b[k]=b[k]/d;
                for( i=k+1;i<=n-1;i++ )
                {
                        for( j=k+1;j<=n-1;j++ )
                        {
                                p=i*n+j;
                                a[p]=a[p]-a[i*n+k]*a[k*n+j];
                 }
                        b[i]=b[i]-a[i*n+k]*b[k];
                }
        }
        d=a[(n-1)*n+n-1];
        if(fabs(d)+1.0==1.0)
        {
                free(js);
                printf("fail\n");
                return(0);
        }
        b[n-1]=b[n-1]/d;
        for( i=n-2;i>=0;i-- )
        {
                t=0.0;
                for( j=i+1;j<=n-1;j++ )
                        t=t+a[i*n+j]*b[j];
                b[i]=b[i]-t;
        }
        js[n-1]=n-1;
        for( k=n-1;k>=0;k-- )
        if( js[k]!=k )
        {
                t=b[k];
                b[k]=b[js[k]];
                b[js[k]]=t;
        }
        free(js);
        return(1);
}

int convolution(float *x,int nx,float *y,int ny,float *z)
{
        int i,j;
        float tmp1[nx];
        float tmp2[ny];

        for(i=0;i<nx+ny-1;i++)
        {
                z[i]=0.0;
        }

        if(nx>ny)
        {
                for(i=0;i<ny;i++)
                {
                        tmp2[i]=y[ny-1-i];
                }
                for(i=0;i<ny;i++)
                {
                        for(j=0;j<=i;j++)
                        {
                                z[i]=z[i]+x[j]*tmp2[ny-1-i+j];
                        }
                }
                for(i=ny;i<nx;i++)
                {
                        for(j=0;j<ny;j++)
                        {
                                z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
                        }
                }
                for(i=nx;i<nx+ny-1;i++)
                {
                        for(j=0;j<ny-i+nx-1;j++)
                        {
                                z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
                        }
                }
        }

        else if(nx==ny)
        {
                for(i=0;i<ny;i++)
                {
                        tmp2[i]=y[ny-1-i];
                }
                for(i=0;i<ny;i++)
                {
                        for(j=0;j<=i;j++)
                        {
                                z[i]=z[i]+x[j]*tmp2[ny-1-i+j];
                        }
                }
                for(i=ny;i<2*ny-1;i++)
                {
                        for(j=0;j<2*ny-i-1;j++)
                        {
                                z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
                        }
                }
        }

        else
        {
                for(i=0;i<nx;i++)
                {
                        tmp1[i]=x[nx-1-i];
                }
                for(i=0;i<nx;i++)
                {
                        for(j=0;j<=i;j++)
                        {
                                z[i]=z[i]+y[j]*tmp1[nx-1-i+j];
                        }
                }
                for(i=nx;i<ny;i++)
                {
                        for(j=0;j<nx;j++)
                        {
                                z[i]=z[i]+tmp1[j]*y[i-nx+1+j];
                        }
                }
                for(i=ny;i<nx+ny-1;i++)
                {
                        for(j=0;j<nx-i+ny-1;j++)
                        {
                                z[i]=z[i]+tmp1[j]*y[i-nx+1+j];
                        }
                }
        }
        return 1;
}


int main()
{
    char input[256],output[256],autocorr[256];
    int i,j,k,l;
    int ntrace,lt,length,lag;
    float lamda;

    ifstream swq1;
    swq1.open("predictive_deconvolution.par");
    if(!swq1)
        {
           cout<<"cannot open"<<endl;
           abort();
        }
   swq1>>input>>output>>autocorr>>ntrace>>lt>>length>>lag>>lamda;
   swq1.close();

   float **u=new float *[ntrace];//initial data
   for(i=0;i<ntrace;i++)
      u[i]=new float [lt];

   float **u_pre=new float *[ntrace];
   for(i=0;i<ntrace;i++)
      u_pre[i]=new float [lt];

   float **final=new float *[ntrace];
   for(i=0;i<ntrace;i++)
      final[i]=new float [lt+length-1];
   
   float **autoco=new float *[ntrace];
   for(i=0;i<ntrace;i++)
      autoco[i]=new float [lt];
  
   float *pre_oper=new float [length];
   float *lag_right=new float [length];  
   float *sintrace=new float [length]; 
   float *final_tmp=new float [lt+length-1];
   float **toep=new float *[length];
   for(i=0;i<length;i++)
       toep[i]=new float [length];
   float *toep1d=new float [lt*lt];

   for(i=0;i<length;i++)
       {
          pre_oper[i]=0.0;
          lag_right[i]=0.0;
          sintrace[i]=0.0;
       }
    
   for(i=0;i<length;i++)
      for(j=0;j<length;j++)
          toep[i][j]=0.0;
 
   for(i=0;i<lt*lt;i++)
      toep1d[i]=0.0;

   ifstream swq2;
   swq2.open(input,ios::binary);
   if(!swq2)
   {
        cout<<"cannot open"<<endl;
        abort();
   }

     for(i=0;i<ntrace;i++)
           for(j=0;j<lt;j++)
                swq2.read((char*)&u[i][j],sizeof(u[i][j]));
   swq2.close();
 
   ifstream swq3;
   swq3.open(autocorr,ios::binary);
   if(!swq3)
   {
        cout<<"cannot open"<<endl;
        abort();
   }

     for(i=0;i<ntrace;i++)
           for(j=0;j<lt;j++)
                swq3.read((char*)&autoco[i][j],sizeof(autoco[i][j]));
   swq3.close();

    ofstream swq4;
    swq4.open("multiple_predicted_by_deconvolution_lag280.dat",ios::binary);
    if(!swq4)
    {
        cout<<"cannot open"<<endl;
        abort();
    }

    ofstream swq5;
    swq5.open(output,ios::binary);
    if(!swq5)
    {
        cout<<"cannot open"<<endl;
        abort();
    }
   
   for(i=0;i<ntrace;i++)
      {
          for(j=0;j<length;j++)
              sintrace[j]=u[i][j];    

//form the right-hand side of Toeplitz equation;
          for(j=0;j<lt-lag;j++)
             lag_right[j]=autoco[i][j+lag];
          for(j=lt-lag;j<lt;j++)
             lag_right[j]=0.0;

//form Toeplitz matrix
          for(j=0;j<lt;j++)
             {
                for(k=j;k<lt;k++)
                    toep[j][k]=autoco[i][k-j];
                for(k=0;k<j;k++)
                    toep[j][k]=autoco[i][j-k];
             }             


          for(j=0;j<lt;j++)
             toep[j][j]*=(1+lamda);
/*
           ofstream swq33;
           swq33.open("test_toep.dat",ios::binary);
           if(!swq33)
            {
               cout<<"cannot open"<<endl;
               abort();
            }

           for(j=0;j<lt;j++)
               for(k=0;k<lt;k++)
                  swq33.write((char*)&toep[j][k],sizeof(toep[j][k]));
           swq33.close();
 
           return 0;
*/

           for(j=0;j<lt;j++)
                {
                  for(k=0;k<lt;k++)
                        {
                            l=j*lt+k;
                            toep1d[l]=toep[j][k];
                        }
                }


         gaus(toep1d,lag_right,lt);
         convolution(sintrace,lt,lag_right,length,final_tmp);

         for(j=0;j<lag;j++)
             u_pre[i][j]=0.0;
         for(j=lag;j<lt;j++)
             u_pre[i][j]=final_tmp[j-lag];
           
         for(j=0;j<lt;j++)
            swq4.write((char*)&u_pre[i][j],sizeof(u_pre[i][j]));

          for(j=0;j<lt;j++)
           final[i][j]=u[i][j]-u_pre[i][j];

         for(j=0;j<lt;j++)
            swq5.write((char*)&final[i][j],sizeof(final[i][j]));
       
         for(j=0;j<lt;j++)
            {
                 sintrace[j]=0.0;
                 lag_right[j]=0.0;
            }
         for(j=0;j<lt+length-1;j++)
                  final_tmp[j]=0.0;
         for(j=0;j<lt;j++)
             for(k=0;k<lt;k++)
                 toep[j][k]=0.0;        
         cout<<i<<"   trace predicted deconvolution done!"<<endl;
 
      }

    swq4.close();
    swq5.close();

    return 0;


}











































