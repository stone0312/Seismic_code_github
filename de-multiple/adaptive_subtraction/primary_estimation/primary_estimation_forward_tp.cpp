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
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn55[256];
  int ns,nr,np,lt;
  float dt,fmin,fmax,min,max;

  ifstream swq;
  swq.open("primary_estimation_forward_tp.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn55>>fn5>>ns>>nr>>np>>lt>>dt>>fmin>>fmax>>min>>max;
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

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1,ip; 
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
  usf=alloc2complex(ltt,nr*np);
  complex<float> **urf;
  urf=alloc2complex(ltt,ns*np);
  complex<float> *urfsp;
  urfsp=alloc1complex(ns*np);

  complex<float> *b;
  b=alloc1complex(np);
  complex<float> **mf;
  mf=alloc2complex(ltt,nr);

  for(ir=0;ir<nr;ir++)
    for(it=0;it<ltt;it++)
      mf[ir][it]=(0.0,0.0);

//ax=b
  complex<float> **a;
  a=alloc2complex(nr*np,nr);

  complex<float> *x;
  x=alloc1complex(nr*np);

  complex<float> **mpf;
  mpf=alloc2complex(ltt,np);

  float **mpt;
  mpt=alloc2float(lt,np);
 
  float *mt;
  mt=alloc1float(lt);

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

  ofstream swq55;
  swq55.open(fn55,ios::binary);
  if(!swq55)
       cout<<"cannot open "<<fn55<<endl;

  ofstream swq5;
  swq5.open(fn5,ios::binary);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;

  for(is=0;is<ns;is++)
//  for(is=0;is<1;is++)
    {
      for(ir=0;ir<np;ir++)
        for(it=0;it<ltt;it++)
           mpf[ir][it]=(0.0,0.0);

      for(ir=0;ir<nr*np;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real())); 
           swq2.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
         }
       
           swq3.seekg(0,ios::beg);
           swq4.seekg(0,ios::beg);

           for(ir=0;ir<nr;ir++)
              {//888
                for(ir1=0;ir1<nr*np;ir1++)
                 for(it1=0;it1<ltt;it1++)
                  {
                    swq3.read((char*)&urf[ir1][it1].real(),sizeof(urf[ir1][it1].real()));
                    swq4.read((char*)&urf[ir1][it1].imag(),sizeof(urf[ir1][it1].imag()));
                  }
                //sparse constraints
                for(it=ifmin;it<ifmax;it++)
        	{//369
          	 for(ir1=0;ir1<nr*np;ir1++)
            	  x[ir1]=0.0;

          	for(ir1=0;ir1<nr;ir1++)
            	 {
                   for(ip=0;ip<np;ip++)
                     x[ip*nr+ir1]=usf[ir1*np+ip][it];
            	 }
           
          	for(ir1=0;ir1<np;ir1++)
            	 for(ir2=0;ir2<nr*np;ir2++)
                    a[ir1][ir2]=(0.0,0.0); 

                for(ir1=0;ir1<nr;ir1++)
                  {
                     for(ip=0;ip<np;ip++)
                       urfsp[ir1*np+ip]=urf[ir1*np+np-1-ip][it];
                  }

                //form matrix a
                for(ip=0;ip<np;ip++)  
                  {
                     for(ir1=ip*nr;ir1<(ip+1)*nr;ir1++)
                       a[ip][ir1]=urfsp[(ir1-ip*nr)*np+ip];                         
                  }

                complex_matrix_multiply_vector(b,a,x,np,nr*np);

                for(ip=0;ip<np;ip++)
                  mpf[ip][it]=b[ip];

                }//369 

             for(ip=0;ip<np;ip++)
       	       for(it=ltt/2+1;it<ltt;it++)
        	{
          	 mpf[ip][it].real()=mpf[ip][ltt-it].real();
          	 mpf[ip][it].imag()=-mpf[ip][ltt-it].imag();
        	}
 
    	     for(ip=0;ip<np;ip++)
     	       {
      		for(it=0;it<ltt;it++)
       		{
         	 in1[it][0]=mpf[ip][it].real();
         	 in1[it][1]=mpf[ip][it].imag();
       		}

      		fftwf_execute(p1);

      		for(it=0;it<lt;it++)
       		 mpt[ip][it]=out1[it][0]/ltt;
     	       }

              for(ip=0;ip<np;ip++)
       	       for(it=0;it<lt;it++)
        	swq55.write((char*)&mpt[ip][it],sizeof(mpt[ip][it]));

//inverse tau-p transform
              for(it=0;it<lt;it++)
                mt[it]=0.0;

              for(it=0;it<lt;it++)
               {
                 for(ip=0;ip<np;ip++)
                   mt[it]+=mpt[ip][it];                  
               }

       	      for(it=0;it<lt;it++)
         	 mt[it]*=-1.0;

       	       for(it=0;it<lt;it++)
        	swq5.write((char*)&mt[it],sizeof(mt[it]));

//      	     cout<<ir+1<<" Trace Forward Prediction Done ..."<<endl;
 
       }//888

      cout<<is+1<<" Shot Forward Prediction Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











