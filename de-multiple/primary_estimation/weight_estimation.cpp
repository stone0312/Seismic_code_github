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

int complex_CG(complex <float> **a, complex <float> *x1, complex <float> *x2,complex <float> *f,complex <float> *p1,complex <float> *p2,complex <float> *r1,complex <float> *r2,complex <float> *ap, complex <float> *apcong,int np, float err)
{
  int ip,ix,iz;

  float alfa1,alfa2,beta1,beta2;
  float r1norm1,r2norm1,r1norm,p1ap,r2norm;
  complex<float> p1ap1_tmp;

  float norm2min=0.00000001;
  int iter_max=30;

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
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256];
  int ns, nr,lt, nr2, coe1;
  float dt,fmin,fmax,lamda;
  float err=0.00000000000001;
  float ampb;
  float amp_ata,amp_atb;

  ifstream swq;
  swq.open("weight_estimation.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn9>>fn10>>fn5>>fn6>>fn7>>fn8>>ns>>nr>>lt>>dt>>fmin>>fmax>>lamda>>coe1;
  swq.close();

  cout<<"Fna of CSG Real is===="<<fn1<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of True Multiple Real is===="<<fn9<<endl; 
  cout<<"Fna of True Multiple Imaginary is===="<<fn10<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
  cout<<"Fna of Estimated Weights is===="<<fn6<<endl; 
  cout<<"Fna of Estimated multiples is===="<<fn7<<endl; 
  cout<<"Fna of Difference between estimated multiples and true multiples is===="<<fn8<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<" , "<<dt<<endl;
  cout<<"minimum and maximum frequency to be calculated are===="<<fmin<<" , "<<fmax<<endl;

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1; 
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
 
  complex<float> **mf;
  mf=alloc2complex(ltt,nr);
  complex<float> *af;
  af=alloc1complex(nr2);
  complex<float> **pp0;
  pp0=alloc2complex(nr2,nr);

  complex<float> *p0;
  p0=alloc1complex(nr);
  complex<float> **p;
  p=alloc2complex(nr,nr);

  complex<float> *mf1d;
  mf1d=alloc1complex(nr);
  float **mt;
  mt=alloc2float(lt,nr);

  complex<float> **a;
  a=alloc2complex(nr2,nr);
  complex<float> **at;
  at=alloc2complex(nr,nr2);
  complex<float> **ata;
  ata=alloc2complex(nr2,nr2);

  complex<float> *b;
  b=alloc1complex(nr);
  complex<float> *atb;
  atb=alloc1complex(nr2);
  
  complex<float> *wf1d;
  wf1d=alloc1complex(nr2);
  complex<float> **wf2d;
  wf2d=alloc2complex(ltt,nr2);
  float **wt2d;
  wt2d=alloc2float(lt,nr2);

  complex<float> **mftrue;
  mftrue=alloc2complex(ltt,nr);

  complex<float> *mcalf1d;
  mcalf1d=alloc1complex(nr);
  complex<float> **mcalf2d;
  mcalf2d=alloc2complex(ltt,nr);
  float **mcalt2d;
  mcalt2d=alloc2float(lt,nr);
  complex<float> **diff2d;
  diff2d=alloc2complex(ltt,nr);
  float **dift2d;
  dift2d=alloc2float(lt,nr);

   complex<float> *x1;
   x1=alloc1complex(nr2);
   complex<float> *p11;
   p11=alloc1complex(nr2);
   complex<float> *p2;
   p2=alloc1complex(nr2);
   complex<float> *r1;
   r1=alloc1complex(nr2);
   complex<float> *r2;
   r2=alloc1complex(nr2);
   complex<float> *ap;
   ap=alloc1complex(nr2);
   complex<float> *apcong;
   apcong=alloc1complex(nr2);

   float *amp1;
   amp1=alloc1float(nr2*nr2);
   float *amp2;
   amp2=alloc1float(nr2);

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
 
//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
           mf[ir][it]=(0.0,0.0);

      swq1.seekg(0,ios::beg);
      swq2.seekg(0,ios::beg);
      swq9.seekg(0,ios::beg);
      swq10.seekg(0,ios::beg);

      for(ir=0;ir<is;ir++)
        {
           swq1.seekg(nr*ltt*4,ios::cur);
           swq2.seekg(nr*ltt*4,ios::cur);
           swq9.seekg(nr*ltt*4,ios::cur);
           swq10.seekg(nr*ltt*4,ios::cur);
        }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real())); 
           swq2.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
           swq9.read((char*)&mftrue[ir][it].real(),sizeof(mftrue[ir][it].real())); 
           swq10.read((char*)&mftrue[ir][it].imag(),sizeof(mftrue[ir][it].imag())); 
         }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
            usf[ir][it].real()*=coe1;
            usf[ir][it].imag()*=coe1;
            mftrue[ir][it].real()*=coe1;
            mftrue[ir][it].imag()*=coe1;
         } 
       
      for(it=ifmin;it<ifmax;it++)
        {//0000
          for(ir=0;ir<nr2;ir++)
             {
               af[ir].real()=1.0;
               af[ir].imag()=0.0;
             }

          for(ir=0;ir<nr;ir++)
            { 
              p0[ir]=usf[ir][it];
              b[ir]=mftrue[ir][it];
            }

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
                  p[ir][ir1]=urf[ir1][it]; 
              }

//form matrix PP0
          for(ir=0;ir<nr;ir++)       
            {
                for(ir1=ir*nr;ir1<(ir+1)*nr;ir1++)
                    pp0[ir][ir1]=p0[ir1-ir*nr]*p[ir][ir1-ir*nr];
            }

         complex_matrix_multiply_vector(mf1d, pp0, af, nr, nr2);

          for(ir=0;ir<nr;ir++)
            mf[ir][it]=mf1d[ir];

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

          for(ir=0;ir<nr;ir++)       
            for(ir1=0;ir1<nr2;ir1++)
               a[ir][ir1]=pp0[ir][ir1];

         for(ir=0;ir<nr2;ir++)
          for(ir1=0;ir1<nr;ir1++)
            {       
              at[ir][ir1].real()=a[ir1][ir].real();
              at[ir][ir1].imag()=-a[ir1][ir].imag();
            }
         complex_matrix_multiply_vector(atb, at, b, nr2, nr);

         amp_atb=0.0;

         for(ir=0;ir<nr2;ir++)
           amp2[ir]=sqrt(atb[ir].real()*atb[ir].real()+atb[ir].imag()*atb[ir].imag());

         for(ir=0;ir<nr2;ir++)
             {
               if(amp2[ir]>amp_atb)
                  amp_atb=amp2[ir];
               else
                  amp_atb=amp_atb;
             }

         cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATB===="<<amp_atb<<endl;

         complex_matrix_multiply_matrix(ata, at, a, nr2, nr, nr2);

         amp_ata=0.0;   

         for(ir=0;ir<nr2;ir++)
          {
            for(ir1=0;ir1<nr2;ir1++)
              amp1[ir*nr2+ir1]=sqrt(ata[ir][ir1].real()*ata[ir][ir1].real()+ata[ir][ir1].imag()*ata[ir][ir1].imag());
          }

         for(ir=0;ir<nr2*nr2;ir++)
          {
            if(amp1[ir]>amp_ata)
               amp_ata=amp1[ir];
            else
               amp_ata=amp_ata;
          }

         cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATA===="<<amp_ata<<endl;

         if(amp_atb!=0.0&&amp_ata!=0.0)
          {//6666
            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Begins..."<<endl;         

             for(ir=0;ir<nr2;ir++)  
               for(ir1=0;ir1<nr2;ir1++)
                 {
                    ata[ir][ir1].real()/=amp_ata;
                    ata[ir][ir1].imag()/=amp_ata;
                 }  

            for(ir=0;ir<nr2;ir++) 
              ata[ir][ir].real()+=lamda; 

            for(ir=0;ir<nr2;ir++)  
               {
                 atb[ir].real()/=amp_atb;
                 atb[ir].imag()/=amp_atb;
               }

            complex_CG(ata, x1, wf1d, atb, p11, p2, r1, r2, ap, apcong, nr2, err); 

            for(ir=0;ir<nr2;ir++) 
              {
                wf1d[ir].real()*=amp_ata/amp_atb;
                wf1d[ir].imag()*=amp_ata/amp_atb;
              }
 
            for(ir=0;ir<nr2;ir++)
             wf2d[ir][it]=wf1d[ir];       

            complex_matrix_multiply_vector(mcalf1d, a, wf1d, nr, nr2);

            for(ir=0;ir<nr;ir++)
              mcalf2d[ir][it]=mcalf1d[ir];

            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;
          }//6666

          else
           {//7777
             for(ir=0;ir<nr;ir++)
               mf[ir][it]=(0.0,0.0);
             cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

             for(ir=0;ir<nr2;ir++)
              wf2d[ir][it]=(0.0,0.0);

             for(ir=0;ir<nr;ir++)
               mcalf2d[ir][it]=(0.0,0.0);

             cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;         
           }//7777

        }//8888

       else
        {
          for(ir=0;ir<nr;ir++)
            mf[ir][it]=(0.0,0.0);
          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

           for(ir=0;ir<nr2;ir++)
            wf2d[ir][it]=(0.0,0.0);

          for(ir=0;ir<nr;ir++)
            mcalf2d[ir][it]=(0.0,0.0);

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;         
        }

    }//0000

      for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mf[ir][it].real()=mf[ir][ltt-it].real();
          mf[ir][it].imag()=-mf[ir][ltt-it].imag();
          mcalf2d[ir][it].real()=mcalf2d[ir][ltt-it].real();
          mcalf2d[ir][it].imag()=-mcalf2d[ir][ltt-it].imag();
        }
      for(ir=0;ir<nr2;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          wf2d[ir][it].real()=wf2d[ir][ltt-it].real();
          wf2d[ir][it].imag()=-wf2d[ir][ltt-it].imag();
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

    for(ir=0;ir<nr2;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=wf2d[ir][it].real();
         in1[it][1]=wf2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       wt2d[ir][it]=out1[it][0]/ltt;
     }

   for(ir=0;ir<nr2;ir++)
       for(it=0;it<lt;it++)
        swq6.write((char*)&wt2d[ir][it],sizeof(wt2d[ir][it]));

    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mcalf2d[ir][it].real();
         in1[it][1]=mcalf2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mcalt2d[ir][it]=out1[it][0]/ltt;
     }

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq7.write((char*)&mcalt2d[ir][it],sizeof(mcalt2d[ir][it]));

    for(ir=0;ir<nr;ir++)
     for(it=0;it<ltt;it++)
       diff2d[ir][it]=mcalf2d[ir][it]-mf[ir][it];

    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=diff2d[ir][it].real();
         in1[it][1]=diff2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       dift2d[ir][it]=out1[it][0]/ltt;
     }

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq8.write((char*)&dift2d[ir][it],sizeof(dift2d[ir][it]));

      cout<<is+1<<" Shot Inverse Estimation Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}










