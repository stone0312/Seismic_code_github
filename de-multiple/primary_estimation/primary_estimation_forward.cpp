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
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256];
  int ns,nr,lt,lim,flag;
  float dt,fmin,fmax,min,max;

  ifstream swq;
  swq.open("primary_estimation_forward.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn5>>ns>>nr>>lt>>dt>>fmin>>fmax>>lim>>flag;
  swq.close();

  cout<<"Fna of CSG Real is===="<<fn1<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<" , "<<dt<<endl;
  cout<<"minimum and maximum frequency to be calculated are===="<<fmin<<" , "<<fmax<<endl;

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1; 
  ltt=2*lt;
 
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
 
  complex<float> **mf;
  mf=alloc2complex(ltt,nr);
  complex<float> *af;
  af=alloc1complex(nr*nr);
  complex<float> **pp0;
  pp0=alloc2complex(nr*nr,nr);

  complex<float> *p0;
  p0=alloc1complex(nr);
  complex<float> **p;
  p=alloc2complex(nr,nr);

  complex<float> *mf1d;
  mf1d=alloc1complex(nr);
  float **mt;
  mt=alloc2float(lt,nr);

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

  ofstream swq5;
  swq5.open(fn5,ios::binary);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;

  if(flag==1)
    cout<<"With MCG Constraints ..."<<endl; 
  else
    cout<<"Without MCG Constraints ..."<<endl; 
 
//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
           mf[ir][it]=(0.0,0.0);

      swq1.seekg(0,ios::beg);
      swq2.seekg(0,ios::beg);

      for(ir=0;ir<is;ir++)
        {
           swq1.seekg(nr*ltt*4,ios::cur);
           swq2.seekg(nr*ltt*4,ios::cur);
        }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real())); 
           swq2.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
         }
       
      for(it=ifmin;it<ifmax;it++)
        {
          for(ir=0;ir<nr*nr;ir++)
             {
               af[ir].real()=1.0;
               af[ir].imag()=0.0;
             }

          for(ir=0;ir<nr;ir++)
             p0[ir]=usf[ir][it];

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

                //form matrix CRG
                for(ir1=0;ir1<nr;ir1++)
                  p[ir][ir1]=urf[ir1][it]; 

              }

//form matrix PP0
          for(ir=0;ir<nr;ir++)       
            {
                for(ir1=ir*nr;ir1<(ir+1)*nr;ir1++)
                    pp0[ir][ir1]=p0[ir1-ir*nr]*p[ir][ir1-ir*nr];
            }

//With MCG aperature constraints
       if(flag==1)
        { 
/* 
          for(ir=0;ir<nr;ir++)
            {
               if(ir<is)
                 {
                    for(ir1=ir*nr;ir1<ir*nr+ir;ir1++)
                       pp0[ir][ir1]=(0.0,0.0);
                    for(ir1=ir*nr+ir;ir1<ir*nr+is+1;ir1++)
                      {
                         pp0[ir][ir1].real()/=(is-ir+1);
                         pp0[ir][ir1].imag()/=(is-ir+1);
                      }
                    for(ir1=ir*nr+is+1;ir1<ir*nr+nr;ir1++)
                       pp0[ir][ir1]=(0.0,0.0);
                 } 
               else
                 {
                    for(ir1=ir*nr;ir1<ir*nr+is;ir1++)
                       pp0[ir][ir1]=(0.0,0.0);
                    for(ir1=ir*nr+is;ir1<ir*nr+ir+1;ir1++)
                      {
                         pp0[ir][ir1].real()/=(ir-is+1);
                         pp0[ir][ir1].imag()/=(ir-is+1);
                      }
                    for(ir1=ir*nr+ir+1;ir1<ir*nr+nr;ir1++)
                       pp0[ir][ir1]=(0.0,0.0);
                 } 
            }
*/

          for(ir=0;ir<nr;ir++)
            {
               if(abs(ir-is)<lim)
                 {//1111
                    if(ir<is)
                      {
                          for(ir1=ir*nr;ir1<ir*nr+ir;ir1++)
	                    pp0[ir][ir1]=(0.0,0.0);
                    	  for(ir1=ir*nr+ir;ir1<ir*nr+is+1;ir1++)
                            {
                               pp0[ir][ir1].real()/=(is-ir+1);
                               pp0[ir][ir1].imag()/=(is-ir+1);
                            }
                          for(ir1=ir*nr+is+1;ir1<ir*nr+nr;ir1++)
                            pp0[ir][ir1]=(0.0,0.0);
                      }
                    else
                     {
                         for(ir1=ir*nr;ir1<ir*nr+is;ir1++)
                           pp0[ir][ir1]=(0.0,0.0);
                         for(ir1=ir*nr+is;ir1<ir*nr+ir+1;ir1++)
                           {
                             pp0[ir][ir1].real()/=(ir-is+1);
                             pp0[ir][ir1].imag()/=(ir-is+1);
                           }
                         for(ir1=ir*nr+ir+1;ir1<ir*nr+nr;ir1++)
                           pp0[ir][ir1]=(0.0,0.0);
                     }
                 }//1111 
               else
                 {//2222
                    for(ir1=ir*nr;ir1<ir*nr+(int)((ir+is-lim)/2);ir1++)
                      pp0[ir][ir1]=(0.0,0.0);
                    for(ir1=ir*nr+(int)((ir+is-lim)/2);ir1<ir*nr+(int)((ir+is+lim)/2);ir1++)
                      {
                        pp0[ir][ir1].real()/=lim;
                        pp0[ir][ir1].imag()/=lim;
                      }
                    for(ir1=ir*nr+(int)((ir+is+lim)/2);ir1<ir*nr+nr;ir1++)
                      pp0[ir][ir1]=(0.0,0.0);                
                 }//2222 
            }


          complex_matrix_multiply_vector(mf1d, pp0, af, nr, nr*nr);

          for(ir=0;ir<nr;ir++)
            mf[ir][it]=mf1d[ir];

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Done..."<<endl;

        } 

     else
      {
         cout<<"Without MCG Constraints ..."<<endl; 

         complex_matrix_multiply_vector(mf1d, pp0, af, nr, nr*nr);

          for(ir=0;ir<nr;ir++)
            mf[ir][it]=mf1d[ir];

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Done..."<<endl;
      }
    }


      for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mf[ir][it].real()=mf[ir][ltt-it].real();
          mf[ir][it].imag()=-mf[ir][ltt-it].imag();
        }
 
    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mf[ir][it].real();
         in1[it][1]=mf[ir][it].imag();
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











