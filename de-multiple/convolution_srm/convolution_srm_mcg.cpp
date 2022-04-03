#include <iostream>
#include <fstream>
#include <complex>

using namespace std;

#include "fftw3.h"
#include "alloc.c"
#include "math.h"


int main()
{
   int i,j,k,is,it,itmp,ir,ir1,ir2;
   int nshot,nt,lt;
   char input[256],inputcp[256],output[256],output_mcg[256];
   int ltt;   

   ifstream swq;
   swq.open("convolution_srm_mcg.par");
   if(!swq)
     {
        cout<<"cannot open convolution_srm_mcg.par"<<endl;
        abort();
     }
   swq>>input>>inputcp>>output>>output_mcg>>nshot>>nt>>lt;
   swq.close();

   cout<<"fna of original shot gather is===="<<input<<endl;
   cout<<"fna of output multiple model is===="<<output<<endl;
   cout<<"fna of output MCG is===="<<output_mcg<<endl;
   cout<<"No. of shots is===="<<nshot<<endl;
   cout<<"No. of traces per shot is===="<<nt<<endl;
   cout<<"No. of temperal samples is===="<<lt<<endl;

   ltt=2*lt;
   
   float **us=new float *[nt];
   for(i=0;i<nt;i++)
     us[i]=new float [ltt];

   float **ur=new float *[nshot];
   for(i=0;i<nshot;i++)
     ur[i]=new float [ltt];
 

   complex<float> **usf=new complex<float> *[nt];
   for(i=0;i<nt;i++)
     usf[i]=new complex<float> [ltt];


   complex<float> **urf=new complex<float> *[nshot];
   for(i=0;i<nshot;i++)
     urf[i]=new complex<float> [ltt];
 
 
   complex<float> *mf=new complex<float> [ltt];
   for(i=0;i<ltt;i++)
      {
         mf[i].real()=0.0;
         mf[i].imag()=0.0;
      }

   complex<float> **mcgf;
  mcgf=alloc2complex(ltt,nt);

  float **mcg;
  mcg=alloc2float(ltt,nt);

   float *m=new float [ltt];
    
   float mem=0.0;
   mem=(nt*ltt+nshot*ltt+nt*ltt*2+nshot*ltt*2+ltt*2+ltt*6*2)*4/(1024.0*1024.0);
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

  ofstream swq4;
  swq4.open(output_mcg,ios::binary);
  if(!swq4)
    {
       cout<<"cannot open "<<output_mcg<<endl;
       abort();
    }


       fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3;
       fftwf_plan p1, p2, p3;
       in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
       out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
       in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
       out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
       in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
       out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);    

       p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
       p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_FORWARD,FFTW_MEASURE);
       p3=fftwf_plan_dft_1d(ltt,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
  
//   for(is=0;is<nshot;is++)
   for(is=0;is<1;is++)
     {
        for(i=0;i<nt;i++)
           for(j=0;j<ltt;j++)
               us[i][j]=0.0;
  
        for(i=0;i<nt;i++)
          for(j=0;j<lt;j++)
             swq1.read((char*)&us[i][j],sizeof(us[i][j]));//read common shot gather 

//fft to frequency domain       
       for(i=0;i<nt;i++)
         {
            if((in1==NULL)||(out1==NULL))
               cout<<"memory insufficient"<<endl;
            else
              {
                 for(j=0;j<ltt;j++)
                    {
                       in1[j][0]=us[i][j];
                       in1[j][1]=0.0;
                    }
              }

           fftwf_execute(p1);

           for(k=0;k<ltt;k++)
             {
                usf[i][k].real()=out1[k][0];
                usf[i][k].imag()=out1[k][1];
             }
         }

        swq2.seekg(0,ios::beg);

        for(it=0;it<nt;it++)
           {
              for(i=0;i<nshot;i++)
                for(j=0;j<ltt;j++)
                   ur[i][j]=0.0;

              for(i=0;i<nshot;i++)//read common receiver gather
                 for(j=0;j<lt;j++)
                     swq2.read((char*)&ur[i][j],sizeof(ur[i][j]));
              
//reordered data fft to frequency domain
              for(i=0;i<nshot;i++)
                {
                   if((in2==NULL)||(out2==NULL))
                      cout<<"memory insufficient"<<endl;
                   else
                     {
                        for(j=0;j<ltt;j++)
                           {
                              in2[j][0]=ur[i][j];
                              in2[j][1]=0.0;
                           }
                     }

                  fftwf_execute(p2);

                  for(k=0;k<ltt;k++)
                    {
                       urf[i][k].real()=out2[k][0];
                       urf[i][k].imag()=out2[k][1];
                    }
                }

              for(ir1=0;ir1<nt;ir1++)
                 for(j=0;j<ltt;j++)
                   {
                      mcgf[ir1][j].real()=0.0;
                      mcgf[ir1][j].imag()=0.0;
                   }

//convolution (equivalent to multiplication in frequency domain) to predict multiple in shot-receiver pair (is,it).
              for(i=0;i<ltt;i++)
                  {
                    mf[i].real()=0.0;
                    mf[i].imag()=0.0;

                    for(ir1=int(0.4*it);ir1<int(0.6*it);ir1++)
                      {  
                        mcgf[ir1][i]=usf[ir1][i]*urf[ir1][i];
                        mf[i]+=usf[ir1][i]*urf[ir1][i];
                      }

                    mf[i].real()/=(it+1);
                    mf[i].imag()/=(it+1);
                  }

                   if((in3==NULL)||(out3==NULL))
                      cout<<"memory insufficient"<<endl;
                   else
                     {
                        for(j=0;j<ltt;j++)
                           {
                              in3[j][0]=mf[j].real();
                              in3[j][1]=mf[j].imag();
                           }
                     }

                  fftwf_execute(p3);

                  for(k=0;k<ltt;k++)
                     m[k]=out3[k][0]/ltt;

                    



                for(ir1=0;ir1<nt;ir1++) 
                 {   
                     {
                        for(j=0;j<ltt;j++)
                           {
                              in3[j][0]=mcgf[ir1][j].real();
                              in3[j][1]=mcgf[ir1][j].imag();
                           }
                     }

                   fftwf_execute(p3);

                   for(j=0;j<ltt;j++)
                     mcg[ir1][j]=out3[j][0]/ltt;
                 } 

              for(i=0;i<lt;i++)
                   swq3.write((char*)&m[i],sizeof(m[i]));
          
              for(ir1=0;ir1<nt;ir1++) 
                for(i=0;i<lt;i++)
                   swq4.write((char*)&mcg[ir1][i],sizeof(mcg[ir1][i]));   

 
            if(it%100==0)
              cout<<"(shot,receiver)====("<<is+1<<","<<it+1<<") SRM finished"<<endl;
           }
          cout<<is+1<<" shots prediction is done!"<<endl;
     } 
     cout<<"ALL DONE!"<<endl;
     swq1.close();
     swq2.close();
     swq3.close();
     swq4.close();

     fftwf_destroy_plan(p1);
     fftwf_destroy_plan(p2);
     fftwf_destroy_plan(p3);


//     free(in1);
//     free(out1);
//     free(in2);
//     free(out2);
//     free(in3);
//     free(out3);

 return 0;

}











