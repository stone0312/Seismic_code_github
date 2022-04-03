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

int trmul( float *a , float *b , int m , int n , int k , float *c )
  { int i,j,l,u;
    for (i=0; i<=m-1; i++)
    for (j=0; j<=k-1; j++)
      { u=i*k+j; c[u]=0.0;
        for (l=0; l<=n-1; l++)
          c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
    return 0;
  }

  int grad( float *a , int n , float *b , float eps ,  float *x )
  { int i,k;
    float *p,*r,*s,*q,alpha,beta,d,e;
    int  trmul(float [],float [],int,int,int,float []);
    p=(float*)malloc(sizeof(float)*n);
    r=(float*)malloc(sizeof(float)*n);
    s=(float*)malloc(sizeof(float)*n);
    q=(float*)malloc(sizeof(float)*n);
    for (i=0; i<=n-1; i++)
      { x[i]=0.0; p[i]=b[i]; r[i]=b[i]; }
    i=0;
    while (i<=n-1)
      { trmul(a,p,n,n,1,s);
        d=0.0; e=0.0;
        for (k=0; k<=n-1; k++)
           { d=d+p[k]*b[k]; e=e+p[k]*s[k]; }
        alpha=d/e;
        for (k=0; k<=n-1; k++)
           x[k]=x[k]+alpha*p[k];
        trmul(a,x,n,n,1,q);
        d=0.0;
        for (k=0; k<=n-1; k++)
          { r[k]=b[k]-q[k]; d=d+r[k]*s[k]; }
        beta=d/e; d=0.0;
        for (k=0; k<=n-1; k++) d=d+r[k]*r[k];
        d=sqrt(d);
        if (d<eps)
          { free(p); free(r); free(s); free(q);return 0;}
        for (k=0; k<=n-1; k++)
           p[k]=r[k]-beta*p[k];
        i=i+1;
        cout<<"ITERATION TIME   "<<i+1<<"   "<<d<<endl;
      }
    free(p); free(r); free(s); free(q);
    return 0;
  }

int normL2(int length, float *vector,float *result)
{
        int k;float tmp;
        tmp=0;
        for (k =0; k< length; k++)
        {
                tmp+=vector[k]*vector[k];
        }
        *result=sqrt(tmp);
        return 0;
}
int L2_product(int length, float *x,float *y,float *result)
{
        int k; float tmp=0.0;
        for (k=0;k<length;k++)
                tmp+=x[k]*y[k];
        *result=tmp;
        return 0;
}

int LSCG_Under(int m,int n,float **A,float *b,float *x,float lamda,float epsilon,int maxiter)
{
        // CG method for regularization equation (AA')mu=b  and set x=A'mu in Underdetermine Case;
        // m : row  n: col and m<=n; 
        // 
        float *r,*q,*p,*z,*mu;
        int i,j,iter, n_yn;
        float nrmL2;
        float alpha, beta;
        float sum , tmp;
        float gamma,gamma1;
        int r0;
        r = (float *) malloc (sizeof(float)*m);
        q = (float *) malloc (sizeof(float)*n);
        p = (float *) malloc (sizeof(float)*m);
        z = (float *) malloc (sizeof(float)*m);
        mu = (float *) malloc (sizeof(float)*m);
        for( i = 0; i < m; i++)
                mu[i]=0.0;    // mu=0.0 as initial guess 
        for (j =0 ; j < n ;j++)
                x[j] = 0.0;
        for(i=0;i<m;i++)
        {
                r[i]=b[i];  // compute the residual
                z[i]=r[i];
                p[i]=z[i];
        }

        // compute the initial residual z ,r, and conjugate basis p
        //iteratively construct the conjugate basis pj
        n_yn = normL2(m,z,&gamma);
        if (0 != n_yn)
        {
                printf ("'norml2' Error.\n");
                return 1;
        }
        for(iter=0;iter<maxiter;iter++)
        {
                //compute the q=Ap;
                for(j=0;j<n;j++)
                {       sum=0.0;
                        for(i=0;i<m;i++)
                                sum+=A[i][j]*p[i];   // q=A'p;
                        q[j]=sum;

                }

                  // the length of z gamma=sqrt(z'*z);
                if (iter==0)
                        r0=gamma;

                L2_product(n,q,q,&tmp);
                if (0 != n_yn)
                {
                        printf ("'L2_product' Error.\n");
                        return 1;
                }   
                sum=0.0;
                for(i=0;i<m;i++)
                        sum+=lamda*p[i]*p[i];
                alpha=gamma*gamma/(tmp+sum);

                for(i=0;i<m;i++)
                        mu[i]+=alpha*p[i];   //mu = mu + alpha*p;
                for(j=0;j<n;j++)
                {
                        x[j]=x[j]+alpha*q[j];  //update x=x+alpha*q;
                }
                sum=0.0;
                for(i=0;i<m;i++)
                        sum+=lamda*p[i]*p[i];
                alpha=gamma*gamma/(tmp+sum);
                for(i=0;i<m;i++)
                        mu[i]+=alpha*p[i];   //mu = mu + alpha*p;
                for(j=0;j<n;j++)
                {
                        x[j]=x[j]+alpha*q[j];  //update x=x+alpha*q;
                }
                sum=0.0;
                for(i=0;i<m;i++)
                {       for (j =0; j< n ; j++)
                                sum+=A[i][j]*x[j];
                        r[i]=b[i]-sum;
                }
                for(i=0;i<m;i++)
                        z[i]=r[i]-lamda*mu[i];
                normL2(m,z,&gamma1);
                if (0 != n_yn)
                {
                        printf ("'normL2' Error.\n");
                        return 1;
                }
                normL2(m,r,&nrmL2);
                if (0 != n_yn)
                {
                        printf ("'normL2' Error.\n");
                        return 1;
                }
                beta=pow(gamma1/gamma,2);
                if(nrmL2<epsilon*r0)
                    break;
                //printf("residual=%f\n",nrmL2);
                for(i=0;i<m;i++)
                        p[i]=z[i]+beta*p[i];

                gamma=gamma1;
        printf ("Iteration Time  %d, %f\n", iter, nrmL2);
        }

        free(p);free(q);free(r);free(z);

        return 0;
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


 //        LSCG_Under(lt,lt,toep,lag_right,pre_oper,0.1,0.001,100);
 //        grad(toep1d,lt,lag_right,0.01,pre_oper);
         gaus(toep1d,lag_right,lt);
 //        convolution(sintrace,lt,pre_oper,length,final_tmp);
         convolution(sintrace,lt,lag_right,length,final_tmp);
         for(j=0;j<lt+length-1;j++)
            u_pre[i][j]=final_tmp[j];

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
 
      }
/*
    for(i=0;i<ntrace;i++)
       for(j=0;j<lt;j++)
           final[i][j]=u[i][j]-u_pre[i][j+length-1];
*/
    ofstream swq4;
    swq4.open("multiple_predicted_by_deconvolution",ios::binary);
    if(!swq4)
    {
        cout<<"cannot open"<<endl;
        abort();
    }

    for(i=0;i<ntrace;i++)
         for(j=0;j<lt+length-1;j++)
            swq4.write((char*)&u_pre[i][j],sizeof(u_pre[i][j]));
    swq4.close();
 /*   
    ofstream swq5;
    swq5.open(output,ios::binary);
    if(!swq5)
    {
        cout<<"cannot open"<<endl;
        abort();
    }

    for(i=0;i<ntrace;i++)
         for(j=0;j<lt;j++)
            swq5.write((char*)&final[i][j],sizeof(final[i][j]));
    swq5.close();
*/
    return 0;


}











































