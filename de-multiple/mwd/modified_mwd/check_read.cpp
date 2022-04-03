#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

int main()
{
   char fn1[256],fn2[256],fn3[256];
   int nx,ns,nr,lt,nsr;
   
   int is,ir,it,iz,is1,ir1,is2,ir2;

    ifstream swqq;
    swqq.open("check_read.par");
    if(!swqq)
      {
          cout<<"cannot open check_read.par"<<endl;
          abort();
      }
    swqq>>fn1>>fn2>>fn3>>nx>>ns>>nr>>lt;
    swqq.close();
 
    nsr=ns*nr;

    cout<<"fna of real part of fft of original data is==== "<<fn1<<endl;
    cout<<"fna of imaginary part of fft of original data is==== "<<fn2<<endl;
    cout<<"fna of checked data is==== "<<fn3<<endl;
    cout<<"No. of vel model width is==== "<<nx<<endl;
    cout<<"No. of shots is==== "<<ns<<endl;
    cout<<"No. of traces per shot is==== "<<nr<<endl;
    cout<<"No. of samples is==== "<<lt<<endl;

    complex <float> **d;
    d=alloc2complex(nr,ns);

    complex <float> **u;
    u=alloc2complex(lt,nsr);
  
    float **uu;
    uu=alloc2float(lt,nsr);

//define the plans and arrays when fftw is complemented
    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
 
    p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

    float mem;
    mem=0.0;
    mem+=(nsr*lt*3+ns*nr*2)*4.0/1024.0/1024.0;
    cout<<"Memory needed to be allocated is==== "<<mem<<"MB"<<endl;

//zeroing all arrays;
    for(is=0;is<nsr;is++)
      for(it=0;it<lt;it++)
        {  
          u[is][it]=(0.0,0.0);
          uu[is][it]=0.0;
        }

   for(is=0;is<ns;is++)
    for(iz=0;iz<nr;iz++)
       d[is][iz]=(0.0,0.0);
 
//read the original data
    ifstream swq1;
    swq1.open(fn1,ios::binary);
    if(!swq1)
       {
          cout<<"cannot open "<<fn1<<endl;
          abort();
       }
    ifstream swq11;
    swq11.open(fn2,ios::binary);
    if(!swq11)
       {
          cout<<"cannot open "<<fn2<<endl;
          abort();
       }

    ofstream swq;
    swq.open(fn3,ios::binary);
    if(!swq)
       {
          cout<<"cannot open "<<fn3<<endl;
          abort();
       } 

    for(it=0;it<lt;it++)    
        {
              if((it+1)%5==0)
                 cout<<it<<"frequency slices have been finished!"<<endl;

//forming the data matrix for a certain frequcncy
              for(is=0;is<ns;is++)
                {
                   for(ir=0;ir<nr;ir++)
                    {
                      swq1.seekg(((is*nr+ir)*lt+it)*4,ios::beg);  
                      swq11.seekg(((is*nr+ir)*lt+it)*4,ios::beg);  
                      swq1.read((char*)&(d[is][ir].real()),sizeof(d[is][ir].real()));
                      swq11.read((char*)&(d[is][ir].imag()),sizeof(d[is][ir].imag()));    
                    }                
                }

              for(is=0;is<ns;is++)
                {
                   for(ir=0;ir<nr;ir++)
                    {
                      u[is*nr+ir][it].real()=d[is][ir].real();
                      u[is*nr+ir][it].imag()=d[is][ir].imag(); 
                    }
                }
         }

      for(is=0;is<nsr;is++)
        {
            for(it=0;it<lt;it++)
               {
                  in2[it][0]=u[is][it].real();
                  in2[it][1]=u[is][it].imag();
               }

            fftwf_execute(p2);
          
            for(it=0;it<lt;it++)
               uu[is][it]=out2[it][0];
        }

       for(is=0;is<nsr;is++)
          for(it=0;it<lt;it++)
            swq.write((char*)&(uu[is][it]),sizeof(uu[is][it]));

      swq1.close();
      swq11.close();
      swq.close(); 

   
  return 0;
}

















