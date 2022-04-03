#include "iostream.h"
#include "fstream.h"
#include "fftw3.h"
#include "math.h"
#include "complex.h"

using namespace std;

int inner_product(complex<float> *a,complex<float> *b,complex<float> *c,int n) //vector c=<vector a,vector b>
{
   int i;
   for(i=0;i<n;i++)
      c[i]=a[i]*b[i];

   return 0;
} 

int main()
{
   int i,j,k,is,it;
   int nshot,nt,lt;
   char input[256],inputcp[256],output[256];
   
   ifstream swq;
   swq.open("convolution_srm.par");
   if(!swq)
     {
        cout<<"cannot open convolution_srm.par"<<endl;
        abort();
     }
   swq>>input>>inputcp>>output>>nshot>>nt>>lt;
   swq.close();

   cout<<"fna of original shot gather is===="<<input<<endl;
   cout<<"fna of output multiple model is===="<<output<<endl;
   cout<<"No. of shots is===="<<nshot<<endl;
   cout<<"No. of traces per shot is===="<<nt<<endl;
   cout<<"No. of temperal samples is===="<<lt<<endl;

   float **us=new float *[nt];
   for(i=0;i<nt;i++)
     us[i]=new float [lt];

   float **ur=new float *[nshot];
   for(i=0;i<nshot;i++)
     ur[i]=new float [lt];
 

   complex<float> **usf=new complex<float> *[nt];
   for(i=0;i<nt;i++)
     usf[i]=new complex<float> [lt];


   complex<float> **urf=new complex<float> *[nshot];
   for(i=0;i<nshot;i++)
     urf[i]=new complex<float> [lt];
 
 
   complex<float> *mf=new complex<float> [lt];
   for(i=0;i<lt;i++)
      {
         mf[i].real()=0.0;
         mf[i].imag()=0.0;
      }

   float *m=new float [lt];
    
   float mem=0.0;
   mem=(nt*lt+nshot*lt+nt*lt*2+nshot*lt*2+lt*2+lt)*4/(1024.0*1024.0);
   cout<<"Memory Needed is===="<<mem<<"MB"<<endl;

  ifstream swq1;
  swq1.open(input,ios::binary);
  if(!swq1)
    {
       cout<<"cannot open "<<input<<endl;
       abort();
    } 
  ifstream swq2;
  swq2.open(inputcp,ios::binary);
   if(!swq2)
    {
       cout<<"cannot open "<<inputcp<<endl;
       abort();
    }

  ofstream swq3;
  swq3.open(output,ios::binary);
  if(!swq3)
    {
       cout<<"cannot open "<<output<<endl;
       abort();
    }

  
  for(is=0;is<nshot;is++)
     {
        for(i=0;i<nt;i++)
          for(j=0;j<lt;j++)
             swq1.read((char*)&us[i][j],sizeof(us[i][j]));//read common shot gather 

//fft to frequency domain       
       fftwf_complex *in1,*out1;
       fftwf_plan p1;
       in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
       out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

       for(i=0;i<nt;i++)
         {
            if((in1==NULL)||(out1==NULL))
               cout<<"memory insufficient"<<endl;
            else
              {
                 for(j=0;j<lt;j++)
                    {
                       in1[j][0]=us[i][j];
                       in1[j][1]=0.0;
                    }
              }
           p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE);

           fftwf_execute(p1);

           for(k=0;k<lt;k++)
             {
                usf[i][k].real()=out1[k][0];
                usf[i][k].imag()=out1[k][1];
             }
         }
      fftwf_destroy_plan(p1);
 
        for(it=0;it<nt;it++)
           {
              swq2.seekg(it*lt*4,ios::beg);
              for(j=0;j<lt;j++)
                 swq2.read((char*)&ur[0][j],sizeof(ur[0][j]));

              for(i=1;i<nshot;i++)//read common receiver gather
                 {
                     swq2.seekg((nt-1)*lt*4,ios::cur);
                     for(j=0;j<lt;j++)
                        swq2.read((char*)&ur[i][j],sizeof(ur[i][j]));
                 }
              
//reordered data fft to frequency domain
              fftwf_complex *in2,*out2;
              fftwf_plan p2;
              in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
              out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

              for(i=0;i<nshot;i++)
                {
                   if((in2==NULL)||(out2==NULL))
                      cout<<"memory insufficient"<<endl;
                   else
                     {
                        for(j=0;j<lt;j++)
                           {
                              in2[j][0]=ur[i][j];
                              in2[j][1]=0.0;
                           }
                     }
                  p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);

                  fftwf_execute(p2);

                  for(k=0;k<lt;k++)
                    {
                       urf[i][k].real()=out2[k][0];
                       urf[i][k].imag()=out2[k][1];
                    }
                }
              fftwf_destroy_plan(p2);

//convolution (equivalent to multiplication in frequency domain) to predict multiple in shot-receiver pair (is,it).
              for(i=0;i<lt;i++)
                  {
                      for(j=0;j<nshot;j++)
                         mf[i]+=usf[j][i]*urf[j][i];
                  }

              fftwf_complex *in3,*out3;
              fftwf_plan p3;
              in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
              out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

                   if((in3==NULL)||(out3==NULL))
                      cout<<"memory insufficient"<<endl;
                   else
                     {
                        for(j=0;j<lt;j++)
                           {
                              in3[j][0]=mf[j].real();
                              in3[j][1]=mf[j].imag();
                           }
                     }
                  p3=fftwf_plan_dft_1d(lt,in3,out3,FFTW_BACKWARD,FFTW_ESTIMATE);

                  fftwf_execute(p3);

                  for(k=0;k<lt;k++)
                     m[k]=out3[k][0];
            
               fftwf_destroy_plan(p3);              

              for(i=0;i<lt;i++)
                  {
                      mf[i].real()=0.0;
                      mf[i].imag()=0.0;
                  }
              
              for(i=0;i<lt;i++)
                   swq3.write((char*)&m[i],sizeof(m[i]));

             for(i=0;i<nshot;i++)
                 for(j=0;j<lt;j++)
                      {
                         ur[i][j]=0.0;
                         urf[i][j].real()=0.0;
                         urf[i][j].imag()=0.0;
                      } 
           if(it%100==0)
              cout<<"(shot,receiver)====("<<is+1<<","<<it+1<<") SRM finished"<<endl;
           }
          cout<<is+1<<" shots prediction is done!"<<endl;
     } 
     cout<<"ALL DONE!"<<endl;
     swq1.close();
     swq2.close();
     swq3.close();

 return 0;

}











