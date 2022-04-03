#include "iostream.h"
#include "fstream.h"
#include "math.h"
void ACF(float *x,float *y,int n)
{
        int i,j;
        float sum;
        for(i=0;i<n;i++)
        {
                y[i]=0.0;
        }

        for(i=0;i<n;i++)
        {
                sum=0;
                for(j=i;j<n;j++)
                {
                        sum=x[j]*x[j-i]+sum;
                }
                y[i]=sum;
        }

}


int main()
{
   int i,j,k,l;
   char input[256],output[256];
   int ntrace,lt;
   
   ifstream swq1;
    swq1.open("autocorrelation.par");
    if(!swq1)
        {
           cout<<"cannot open"<<endl;
           abort();
        }
   swq1>>input>>output>>ntrace>>lt;
   swq1.close();

   float **u=new float *[ntrace];
   for(i=0;i<ntrace;i++)
      u[i]=new float [lt];
  
   float *tmp1=new float [lt];
   float *tmp2=new float [lt];
   float **acf=new float *[ntrace];
   for(i=0;i<ntrace;i++)
      acf[i]=new float [lt];

   for(i=0;i<lt;i++)
       {
          tmp1[i]=0.0;
          tmp2[i]=0.0;
       } 
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

   for(i=0;i<ntrace;i++)
      {
         for(j=0;j<lt;j++)
             tmp1[j]=u[i][j];
         ACF(tmp1,tmp2,lt);
         for(k=0;k<lt;k++)
             acf[i][k]=tmp2[k];
         for(l=0;l<lt;l++)
             {
                tmp1[l]=0.0;
                tmp2[l]=0.0;
             }       
   
      }      
  
    ofstream swq3;
    swq3.open(output,ios::binary);
    if(!swq3)
    {
        cout<<"cannot open"<<endl;
        abort();
    }
     
    for(i=0;i<ntrace;i++)
         for(j=0;j<lt;j++)
            swq3.write((char*)&acf[i][j],sizeof(acf[i][j]));
    swq3.close();

    return 0;
  
}




















