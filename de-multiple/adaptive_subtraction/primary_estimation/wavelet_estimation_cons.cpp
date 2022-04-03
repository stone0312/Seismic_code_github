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

//a=b*c m*n,n*p--->m*p
int complex_matrix_multiply_matrix(complex<float> **a, complex<float>**b, complex<float>**c, int m, int n, int p)
{
   int ix,iy,iz;
   for(ix=0;ix<m;ix++)
     for(iy=0;iy<p;iy++)
       a[ix][iy]=(0.0,0.0);

   for(ix=0;ix<m;ix++)
    {
      for(iy=0;iy<p;iy++)
        {
          for(iz=0;iz<n;iz++)
            a[ix][iy]+=b[ix][iz]*c[iz][iy];
        }
    }

  return 0;

}


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

int complex_CG(complex <float> **a, complex <float> *x1, complex <float> *x2,complex <float> *f,complex <float> *p1,complex <float> *p2,complex <float> *r1,complex <float> *r2,complex <float> *ap, complex <float> *apcong,int np, float err)
{
  int ip,ix,iz;

  float alfa1,alfa2,beta1,beta2;
  float r1norm1,r2norm1,r1norm,p1ap,r2norm;
  complex<float> p1ap1_tmp;

  float norm2min=0.00000000001;
  int iter_max=1000;

  for(ip=0;ip<np;ip++)
    {
      x1[ip].real()=0.0;
      x1[ip].imag()=0.0;
      x2[ip].real()=0.0;
      x2[ip].imag()=0.0;
      r1[ip].real()=0.0;
      r1[ip].imag()=0.0;
      r2[ip].real()=0.0;
      r2[ip].imag()=0.0;
      p1[ip].real()=0.0;
      p1[ip].imag()=0.0;
      p2[ip].real()=0.0;
      p2[ip].imag()=0.0;
    }

  for(ip=0;ip<np;ip++)
    {
      r1[ip]=f[ip];
      p1[ip]=r1[ip];
    }

  int iteration=1;

  while(iteration<iter_max)
    {
       r1norm1=0.0;
       for(ip=0;ip<np;ip++)
         r1norm1+=(r1[ip].real()*r1[ip].real()+r1[ip].imag()*r1[ip].imag());
       r1norm=r1norm1;

       complex_matrix_multiply_vector(ap,a,p1,np,np);


       p1ap1_tmp=(0.0,0.0);

       for(ip=0;ip<np;ip++)
       {
         apcong[ip].real()=ap[ip].real();
         apcong[ip].imag()=-ap[ip].imag();
       }

      for(ip=0;ip<np;ip++)
        p1ap1_tmp+=p1[ip]*apcong[ip];
      p1ap=p1ap1_tmp.real();


       if(fabs(p1ap)<norm2min)
         iteration=iter_max;
       else
        {
          alfa1=r1norm/p1ap;

//          cout<<iteration<<" , "<<alfa1<<endl;

          for(ip=0;ip<np;ip++)
           {
             x2[ip].real()=x1[ip].real()+alfa1*p1[ip].real();
             x2[ip].imag()=x1[ip].imag()+alfa1*p1[ip].imag();
           }

          for(ip=0;ip<np;ip++)
           {
             r2[ip].real()=r1[ip].real()-alfa1*ap[ip].real();
             r2[ip].imag()=r1[ip].imag()-alfa1*ap[ip].imag();
           }

          r2norm1=0.0;
          for(ip=0;ip<np;ip++)
            r2norm1+=(r2[ip].real()*r2[ip].real()+r2[ip].imag()*r2[ip].imag());
          r2norm=r2norm1;

          cout<<"Iteration Time===="<<iteration<<" , Error is==== "<<sqrt(r2norm)<<endl;

          beta2=r2norm/r1norm;
          for(ip=0;ip<np;ip++)
            {
              p2[ip].real()=r2[ip].real()+beta2*p1[ip].real();
              p2[ip].imag()=r2[ip].imag()+beta2*p1[ip].imag();
            }
          for(ip=0;ip<np;ip++)
            {
              x1[ip]=x2[ip];
              r1[ip]=r2[ip];
              p1[ip]=p2[ip];
            }

          if(sqrt(r2norm)<err)
             iteration=iter_max;
          else
             iteration+=1;


        }

    }

  return 0;

}


int main()
{
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn11[256],fn22[256],fn6[256],fn7[256],fn8[256],fn9[256];
  int ns,nr,lt,iter,iter_max;
  float dt,fmin,fmax,min,max,lamda1;
  int n1,n2,irmin,irmax,mcg_thre;
  float amp11;
  float err=0.000000001;

  ifstream swq;
  swq.open("wavelet_estimation.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn11>>fn22>>fn6>>fn7>>fn8>>fn9>>ns>>nr>>lt>>dt>>fmin>>fmax>>min>>max>>lamda1>>iter_max>>mcg_thre;
  swq.close();

  cout<<"Fna of CSG Real is===="<<fn1<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<" , "<<dt<<endl;
  cout<<"minimum and maximum frequency to be calculated are===="<<fmin<<" , "<<fmax<<endl;

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1; 
  ltt=2*lt;

  float namp,namp1,namp2,namp3;
  float *amp1;
  amp1=alloc1float(nr*nr);
  float *amp2;
  amp2=alloc1float(nr);
  float *amp3;
  amp3=alloc1float(nr*nr);
 
  float *omega;
  omega=alloc1float(ltt);
  for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
  for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));
 
  ifmin=int(fmin*dt*ltt/1000);
  ifmax=int(fmax*dt*ltt/1000)+1;

  cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;
 
  complex<float> **msf;
  msf=alloc2complex(ltt,nr);

  complex<float> **usf;
  usf=alloc2complex(ltt,nr);
  complex<float> **urf;
  urf=alloc2complex(ltt,ns);

  complex<float> *wi;
  wi=alloc1complex(nr);

  complex<float> *afor;
  afor=alloc1complex(nr);
  complex<float> *mfor1d;
  mfor1d=alloc1complex(nr);
  complex<float> **mfor2d;
  mfor2d=alloc2complex(ltt,nr);
  float **mfor2dt;
  mfor2dt=alloc2float(lt,nr);

  complex<float> *mcal1d;
  mcal1d=alloc1complex(nr);
  complex<float> **mcal2d;
  mcal2d=alloc2complex(ltt,nr);
  float **mcal2dt;
  mcal2dt=alloc2float(lt,nr);

  complex<float> **wi2;
  wi2=alloc2complex(ltt,nr);
  float **wit;
  wit=alloc2float(lt,nr);

  complex<float> **diff;
  diff=alloc2complex(ltt,nr);
  float **difft;
  difft=alloc2float(lt,nr);

  complex<float> **mf;
  mf=alloc2complex(ltt,nr);
  for(ir=0;ir<nr;ir++)
    for(it=0;it<ltt;it++)
      mf[ir][it]=(0.0,0.0);

//ax=b
  complex<float> **a;
  a=alloc2complex(nr,nr);
  complex<float> **araw;
  araw=alloc2complex(nr,nr);
  complex<float> *x;
  x=alloc1complex(nr);
  complex<float> *b;
  b=alloc1complex(nr);
 
  float **mt;
  mt=alloc2float(lt,nr);

//arrays for primary estimation
  complex<float> **a1;   
  a1=alloc2complex(nr,nr);
  complex<float> *x1;
  x1=alloc1complex(nr);
  complex<float> *w;
  w=alloc1complex(nr);
  complex<float> *b1;
  b1=alloc1complex(nr);

  complex<float> **at;
  at=alloc2complex(nr,nr);
  complex<float> **ata;
  ata=alloc2complex(nr,nr);
  complex<float> *atb;
  atb=alloc1complex(nr);


//arrays for inverse wavelet estimation
  complex<float> **a2;
  a2=alloc2complex(nr,nr);
  complex<float> *x2;
  x2=alloc1complex(nr);
  complex<float> *b2;
  b2=alloc1complex(nr);

   complex<float> *p11;
   p11=alloc1complex(nr);
   complex<float> *p2;
   p2=alloc1complex(nr);
   complex<float> *r1;
   r1=alloc1complex(nr);
   complex<float> *r2;
   r2=alloc1complex(nr);
   complex<float> *ap;
   ap=alloc1complex(nr);
   complex<float> *apcong;
   apcong=alloc1complex(nr);


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

  ifstream swq11;
  swq11.open(fn11,ios::binary);
  if(!swq11)
       cout<<"cannot open "<<fn11<<endl; 

  ifstream swq22;
  swq22.open(fn22,ios::binary);
  if(!swq22)
       cout<<"cannot open "<<fn22<<endl;

  ifstream swq3;
  swq3.open(fn3,ios::binary);
  if(!swq3)
       cout<<"cannot open "<<fn3<<endl;

  ifstream swq4;
  swq4.open(fn4,ios::binary);
  if(!swq4)
       cout<<"cannot open "<<fn4<<endl;

  ofstream swq6;
  swq6.open(fn6,ios::binary);
  if(!swq6)
       cout<<"cannot open "<<fn6<<endl;

  ofstream swq7;
  swq7.open(fn7,ios::binary);
  if(!swq7)
       cout<<"cannot open "<<fn7<<endl;

  ofstream swq8;
  swq8.open(fn8,ios::binary);
  if(!swq8)
       cout<<"cannot open "<<fn8<<endl;

  ofstream swq9;
  swq9.open(fn9,ios::binary);
  if(!swq9)
       cout<<"cannot open "<<fn9<<endl;

  cout<<"Primary Estimation Begins..."<<endl;

//  for(is=0;is<ns;is++)
  for(is=0;is<1;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           mcal2d[ir][it]=(0.0,0.0);
           mfor2d[ir][it]=(0.0,0.0);
         }
      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
           wi2[ir][it]=(0.0,0.0);

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real())); 
           swq2.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
           swq11.read((char*)&msf[ir][it].real(),sizeof(msf[ir][it].real())); 
           swq22.read((char*)&msf[ir][it].imag(),sizeof(msf[ir][it].imag())); 
         }
       
    for(it=ifmin;it<ifmax;it++)
     {
       for(ir=0;ir<nr;ir++)
        b1[ir]=usf[ir][it];
       for(ir=0;ir<nr;ir++)
        b[ir]=msf[ir][it];
   
       amp11=0.0;
       for(ir=0;ir<nr;ir++)
         amp11+=pow(b[ir].real(),2)+pow(b[ir].imag(),2);

       if(amp11<pow(10,-12))
        {
           for(ir=0;ir<nr;ir++)
            { 
              mcal2d[ir][it]=(0.0,0.0);
              mfor2d[ir][it]=(0.0,0.0);
            }
        }
 
       else
        {//1111
          for(ir1=0;ir1<nr;ir1++)
            for(ir2=0;ir2<nr;ir2++)
             {
               araw[ir1][ir2]=(0.0,0.0); 
               a[ir1][ir2]=(0.0,0.0); 
             }

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

//                for(ir1=0;ir1<nr;ir1++)
                for(ir1=0;ir1<ir+1;ir1++)
                  araw[ir][ir1]=urf[ir1][it]; 
              }
/*
            for(ir=0;ir<nr;ir++)
             {

                if(ir<mcg_thre)
                  {
                     irmin=0;
                     irmax=ir;
                  }
                else
                  {
                     irmin=(int)(ir*min);
                     irmax=(int)(ir*max+0.5);
                  }
               for(ir1=irmin;ir1<irmax+1;ir1++)
                 a[ir][ir1]=araw[ir][ir1];
             }
*/

          for(ir=0;ir<nr;ir++)
             {
                if(ir<mcg_thre)
                  {
                    for(ir1=0;ir1<ir+1;ir1++)
                      a[ir][ir1]=araw[ir][ir1];      
                  }
                else
                  {
                    for(ir1=0;ir1<int(ir/2)+1;ir1++)
                      {
                         a[ir][ir1].real()=araw[ir][ir1].real()*((float)(ir1)/(float)(int(ir/2)));
                         a[ir][ir1].imag()=araw[ir][ir1].imag()*((float)(ir1)/(float)(int(ir/2)));
                      }
                    for(ir1=int(ir/2)+1;ir1<ir+1;ir1++)
                      {
                         a[ir][ir1].real()=araw[ir][ir1].real()*(1.0-(float)(ir1-int(ir/2))/(float)(int(ir/2)-1));
                         a[ir][ir1].imag()=araw[ir][ir1].imag()*(1.0-(float)(ir1-int(ir/2))/(float)(int(ir/2)-1));
                      }
                  }
             }

           namp3=0.0;
           for(ir1=0;ir1<nr;ir1++)
            {
              for(ir2=0;ir2<nr;ir2++)
                amp3[ir1*nr+ir2]=sqrt(a[ir1][ir2].real()*a[ir1][ir2].real()+a[ir1][ir2].imag()*a[ir1][ir2].imag());
            }
           for(ir1=0;ir1<nr*nr;ir1++)
             {
               if(amp3[ir1]>namp3)
                  namp3=amp3[ir1];
               else
                  namp3=namp3;
             }
         
         if(namp3!=0.0)
          {//3333   
           for(ir1=0;ir1<nr;ir1++)
            for(ir2=0;ir2<nr;ir2++)
               {
                 a[ir1][ir2].real()/=namp3;
                 a[ir1][ir2].imag()/=namp3;
               }

           for(ir1=0;ir1<nr;ir1++)
             {
               for(ir2=0;ir2<nr;ir2++)
                 a1[ir2][ir1]=a[ir2][ir1]*b1[ir1];
             }


           if(it>55&&it<59)
             {
                cout<<"===="<<it+1<<endl;
                cout<<"====Vector B is===="<<endl;
                for(ir1=0;ir1<nr;ir1++)
                   cout<<b[ir1]<<"   ";
                cout<<endl; 

                cout<<"====Matrix A is===="<<endl;
                for(ir1=0;ir1<nr;ir1++)
                  {
                      cout<<ir1<<"===="<<endl;
                      for(ir2=0;ir2<nr;ir2++)
                         cout<<a1[ir1][ir2]<<"    ";
                      cout<<endl;
                  }

             }

           for(ir1=0;ir1<nr;ir1++)
             {
                afor[ir1].real()=1.0;
                afor[ir1].imag()=0.0;
             }

           complex_matrix_multiply_vector(mfor1d,a1,afor,nr,nr);

           for(ir=0;ir<nr;ir++)
              mfor2d[ir][it]=mfor1d[ir]; 

           for(ir1=0;ir1<nr;ir1++)
             for(ir2=0;ir2<nr;ir2++)
               {
                  at[ir1][ir2].real()=a1[ir2][ir1].real();
                  at[ir1][ir2].imag()=-a1[ir2][ir1].imag();
               }

           complex_matrix_multiply_matrix(ata,at,a1,nr,nr,nr); 

           for(ir1=0;ir1<nr;ir1++)
            {
              for(ir2=0;ir2<nr;ir2++)
                amp1[ir1*nr+ir2]=sqrt(ata[ir1][ir2].real()*ata[ir1][ir2].real()+ata[ir1][ir2].imag()*ata[ir1][ir2].imag());
            }
            
           namp1=0.0;
           for(ir1=0;ir1<nr*nr;ir1++)
             {
               if(amp1[ir1]>namp1)
                  namp1=amp1[ir1];
               else
                  namp1=namp1;
             }

          cout<<"Frequency Slice "<<it+1<<" , "<<namp1<<endl;

          for(ir1=0;ir1<nr;ir1++)
            for(ir2=0;ir2<nr;ir2++)
              {
                 ata[ir1][ir2].real()/=namp1;
                 ata[ir1][ir2].imag()/=namp1;
              }

           for(ir1=0;ir1<nr;ir1++)
             ata[ir1][ir1].real()+=lamda1;

           complex_matrix_multiply_vector(atb,at,b,nr,nr);
           
           for(ir1=0;ir1<nr;ir1++)
            amp2[ir1]=sqrt(atb[ir1].real()*atb[ir1].real()+atb[ir1].imag()*atb[ir1].imag());

           namp2=0.0;
           for(ir1=0;ir1<nr;ir1++)
             {
               if(amp2[ir1]>namp2)
                  namp2=amp2[ir1];
               else
                  namp2=namp2;
             }

          cout<<"Frequency Slice "<<it+1<<" , "<<namp2<<endl;
         
            for(ir1=0;ir1<nr;ir1++)
              {
                 atb[ir1].real()/=namp2;
                 atb[ir1].imag()/=namp2;
              }

       if(namp1!=0.0&&namp2!=0.0)
         {//2222

           namp=namp2/namp1; 
           cout<<"Frequency Slice "<<it+1<<" , "<<namp<<endl;

           if(it>55&&it<59)
             {
                cout<<"===="<<it+1<<endl;
                cout<<"====Vector ATB is===="<<endl;
                for(ir1=0;ir1<nr;ir1++)
                   cout<<atb[ir1]<<"   ";
                cout<<endl; 

                cout<<"====Matrix ATA is===="<<endl;
                for(ir1=0;ir1<nr;ir1++)
                  {
                      cout<<ir1<<"===="<<endl;
                      for(ir2=0;ir2<nr;ir2++)
                         cout<<ata[ir1][ir2]<<"    ";
                      cout<<endl;
                  }
             }

           complex_CG(ata, x1, wi, atb, p11, p2, r1, r2, ap, apcong, nr, err);

           cout<<"Old Inverse Wavelet is===="<<endl;
           for(ir1=0;ir1<nr;ir1++)
              cout<<wi[ir1]<<"   ";
           cout<<endl;


           for(ir1=0;ir1<nr;ir1++)
              {
                wi[ir1].real()*=namp;
                wi[ir1].imag()*=namp;
              }

           cout<<"New Inverse Wavelet is===="<<endl;
           for(ir1=0;ir1<nr;ir1++)
              cout<<wi[ir1]<<"   ";
           cout<<endl;

           for(ir=0;ir<nr;ir++)
              wi2[ir][it]=wi[ir]; 

           complex_matrix_multiply_vector(mcal1d,a1,wi,nr,nr);
           
           for(ir=0;ir<nr;ir++)
              mcal2d[ir][it]=mcal1d[ir]; 

         }//2222
        else
         {
            for(ir=0;ir<nr;ir++)
              wi2[ir][it]=(0.0,0.0);
            for(ir=0;ir<nr;ir++)
              mcal2d[ir][it]=(0.0,0.0);
         }

        }//3333
        else
         {
            for(ir=0;ir<nr;ir++)
              wi2[ir][it]=(0.0,0.0);
            for(ir=0;ir<nr;ir++)
              mcal2d[ir][it]=(0.0,0.0);
         }

       }//1111  
      
          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Done..."<<endl;

      } 

      for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mcal2d[ir][it].real()=mcal2d[ir][ltt-it].real();
          mcal2d[ir][it].imag()=-mcal2d[ir][ltt-it].imag();
          mfor2d[ir][it].real()=mfor2d[ir][ltt-it].real();
          mfor2d[ir][it].imag()=-mfor2d[ir][ltt-it].imag();
          wi2[ir][it].real()=wi2[ir][ltt-it].real();
          wi2[ir][it].imag()=-wi2[ir][ltt-it].imag();
        }
 

   for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mfor2d[ir][it].real();
         in1[it][1]=mfor2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mfor2dt[ir][it]=out1[it][0];
     }

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq9.write((char*)&mfor2dt[ir][it],sizeof(mfor2dt[ir][it]));

   for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=wi2[ir][it].real();
         in1[it][1]=wi2[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       wit[ir][it]=out1[it][0];
     }


      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq6.write((char*)&wit[ir][it],sizeof(wit[ir][it]));

   for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mcal2d[ir][it].real();
         in1[it][1]=mcal2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mcal2dt[ir][it]=out1[it][0];
     }
      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq7.write((char*)&mcal2dt[ir][it],sizeof(mcal2dt[ir][it]));

   for(ir=0;ir<nr;ir++) 
     for(it=0;it<ltt;it++)
        diff[ir][it]=mcal2d[ir][it]-msf[ir][it];

   for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=diff[ir][it].real();
         in1[it][1]=diff[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       difft[ir][it]=out1[it][0];
     }
      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq8.write((char*)&difft[ir][it],sizeof(difft[ir][it]));

      cout<<is+1<<" Shot Wavelet Estimation Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











