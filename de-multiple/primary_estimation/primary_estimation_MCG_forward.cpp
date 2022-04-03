using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

//subroutine a=b*c, a is m dimision vector, b is m*n matrix, c is n dimision vector
int complex_matrix_multiply_vector(complex<float> *a, complex<float>**b, complex<float>*c, int m, int n)
{
   int ix,iz;
   for(ix=0;ix<m;ix++)
      a[ix]=(0.0,0.0);

   for(ix=0;ix<m;ix++)
     {
       for(iz=0;iz<n;iz++)
         a[ix]+=b[ix][iz]*c[iz];
     }

   return 0;

}


int main()
{
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256];
  int ns, nr,lt, nr2, coe1,l,irtmp1,irtmp2;
  float dt,fmin,fmax,lamda;
  float ampb;

  ifstream swq;
  swq.open("primary_estimation_MCG_forward.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn9>>fn10>>fn5>>ns>>nr>>lt>>dt>>fmin>>fmax>>coe1>>l;
  swq.close();

  cout<<"Fna of Weights Real is===="<<fn1<<endl; 
  cout<<"Fna of Weights Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of CSG Real is===="<<fn9<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn10<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<" , "<<dt<<endl;
  cout<<"minimum and maximum frequency to be calculated are===="<<fmin<<" , "<<fmax<<endl;
  cout<<"MCG Aperture is===="<<l<<endl; 

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1,irtmp; 
  ltt=2*lt;
  
  nr2=nr*nr;

  float *omega;
  omega=alloc1float(ltt);
  for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
  for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));
 
  ifmin=int(fmin*dt*ltt/1000);
  ifmax=int(fmax*dt*ltt/1000)+1;

  cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;
 
  complex<float> **usf;
  usf=alloc2complex(ltt,nr);
  complex<float> **urf;
  urf=alloc2complex(ltt,nr);
  complex<float> **uswf;
  uswf=alloc2complex(ltt,nr);
  complex<float> **wf;
  wf=alloc2complex(nr,nr);
 
  complex<float> **pp0;
  pp0=alloc2complex(nr,nr);

  complex<float> *mf1d;
  mf1d=alloc1complex(nr);
  complex<float> **mf2d;
  mf2d=alloc2complex(ltt,nr);
  float **mt;
  mt=alloc2float(lt,nr);

  complex<float> *b;
  b=alloc1complex(nr);

  complex<float> **a;
  a=alloc2complex(nr,nr);

   float *amp1;
   amp1=alloc1float(nr*nr);
   float *amp2;
   amp2=alloc1float(nr);

  fftwf_complex *in1,*out1;
  fftwf_plan p1;
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);
 
  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl; 

  ifstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;

  ifstream swq3;
  swq3.open(fn3,ios::binary);
  if(!swq3)
       cout<<"cannot open "<<fn3<<endl;

  ifstream swq4;
  swq4.open(fn4,ios::binary);
  if(!swq4)
       cout<<"cannot open "<<fn4<<endl;

  ifstream swq9;
  swq9.open(fn9,ios::binary);
  if(!swq9)
       cout<<"cannot open "<<fn9<<endl;

  ifstream swq10;
  swq10.open(fn10,ios::binary);
  if(!swq10)
       cout<<"cannot open "<<fn10<<endl;

  ofstream swq5;
  swq5.open(fn5,ios::binary);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;
 
//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
           mf2d[ir][it]=(0.0,0.0);

      swq9.seekg(0,ios::beg);
      swq10.seekg(0,ios::beg);

      for(ir=0;ir<is;ir++)
        {
           swq9.seekg(nr*ltt*4,ios::cur);
           swq10.seekg(nr*ltt*4,ios::cur);
        }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           swq9.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real())); 
           swq10.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
         }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
            usf[ir][it].real()*=coe1;
            usf[ir][it].imag()*=coe1;
         } 
       
      for(it=ifmin;it<ifmax;it++)
        {//0000
          for(ir=0;ir<nr;ir++)
            b[ir]=usf[ir][it];

          for(ir=0;ir<nr;ir++)
           for(ir1=0;ir1<nr;ir1++)
             pp0[ir][ir1]=(0.0,0.0);

          ampb=0.0;
          for(ir=0;ir<nr;ir++)
           ampb+=b[ir].real()*b[ir].real()+b[ir].imag()*b[ir].imag();

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude===="<<ampb<<endl;

          if(ampb!=0)
          {//8888
           swq3.seekg(0,ios::beg);
           swq4.seekg(0,ios::beg);

           for(ir=0;ir<nr;ir++)
              {
                for(ir1=0;ir1<nr;ir1++)
                 for(it1=0;it1<ltt;it1++)
                  {
                    swq3.read((char*)&urf[ir1][it1].real(),sizeof(urf[ir1][it1].real()));
                    swq4.read((char*)&urf[ir1][it1].imag(),sizeof(urf[ir1][it1].imag()));
                  }

                for(ir1=0;ir1<nr;ir1++)
                 for(it1=0;it1<ltt;it1++)
                  {
                     urf[ir1][it1].real()*=coe1;
                     urf[ir1][it1].imag()*=coe1;
                  }

                //form matrix CRG
                for(ir1=0;ir1<nr;ir1++)
                  pp0[ir][ir1]=urf[ir1][it]; 
              }
          
          for(ir=0;ir<nr;ir++)
            for(ir1=0;ir1<nr;ir1++)
              wf[ir][ir1]=(0.0,0.0);

          swq1.seekg(0,ios::beg);
          swq2.seekg(0,ios::beg);

          for(ir=0;ir<nr;ir++)
             {
                for(ir1=0;ir1<nr;ir1++)
                 for(it1=0;it1<ltt;it1++)
                  {
                    swq1.read((char*)&uswf[ir1][it1].real(),sizeof(uswf[ir1][it1].real()));
                    swq2.read((char*)&uswf[ir1][it1].imag(),sizeof(uswf[ir1][it1].imag()));
                  }

                for(ir1=0;ir1<nr;ir1++)
           	 for(it1=0;it1<ltt;it1++)
             	  {
                    uswf[ir1][it1].real()*=coe1;
                    uswf[ir1][it1].imag()*=coe1;
                  } 
        
                for(ir1=0;ir1<nr;ir1++)
                  wf[ir][ir1]=uswf[ir1][it];
             } 

         for(ir=0;ir<nr;ir++)
          for(ir1=0;ir1<nr;ir1++)       
             a[ir][ir1]=pp0[ir][ir1]*wf[ir][ir1];

            if(it==51)
              {
                 cout<<"    "<<is+1<<" Shot, "<<it+1<<"  Frequency Slice A is===="<<endl;
                 for(ir=0;ir<nr;ir++)
                    {
                       cout<<"===="<<ir+1<<endl;
                       for(ir1=0;ir1<nr;ir1++)
                          cout<<a[ir][ir1]<<"     ";
                       cout<<endl;
                    }
                 cout<<"    "<<is+1<<" Shot, "<<it+1<<"  Frequency Slice Weight Matrix is===="<<endl;
                 for(ir=0;ir<nr;ir++)
                    {
                       cout<<"===="<<ir+1<<endl;
                       for(ir1=0;ir1<nr;ir1++)
                          cout<<wf[ir][ir1]<<"     ";
                       cout<<endl;
                    }
              }

         complex_matrix_multiply_vector(mf1d, a , b , nr, nr);

          for(ir=0;ir<nr;ir++)
            mf2d[ir][it]=mf1d[ir];

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;
        }//8888

       else
        {
           for(ir=0;ir<nr;ir++)
             mf2d[ir][it]=(0.0,0.0);
           cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;
        }

    }//0000

      for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mf2d[ir][it].real()=mf2d[ir][ltt-it].real();
          mf2d[ir][it].imag()=-mf2d[ir][ltt-it].imag();
        }

    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mf2d[ir][it].real();
         in1[it][1]=mf2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mt[ir][it]=out1[it][0]/ltt;
     }

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
         mt[ir][it]*=-1.0;

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq5.write((char*)&mt[ir][it],sizeof(mt[ir][it]));

      cout<<is+1<<" Shot Forward Prediction Done ..."<<endl;


    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











