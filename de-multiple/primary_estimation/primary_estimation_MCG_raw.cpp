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

//  float norm2min=0.00000000001;
  float norm2min=pow(10.0,-20);
  int iter_max=1;

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


//      cout<<"P1AP===="<<fabs(p1ap)<<endl;


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

//      cout<<"ITE Time ===="<<iteration<<endl;


    }

  return 0;

}

int main()
{
  char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256];
  int ns, nr,lt, nr2, coe1,l,irtmp1,irtmp2;
  float dt,fmin,fmax,lamda;
  float err=0.00000000000001;
  float ampb;
  float amp_ata,amp_atb;

  ifstream swq;
  swq.open("primary_estimation_MCG.par");
  swq>>fn1>>fn2>>fn3>>fn4>>fn9>>fn10>>fn5>>fn6>>ns>>nr>>lt>>dt>>fmin>>fmax>>lamda>>coe1>>l;
  swq.close();

  cout<<"Fna of Weights Real is===="<<fn1<<endl; 
  cout<<"Fna of Weights Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of CSG Real is===="<<fn9<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn10<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
  cout<<"Fna of Estimated Primary is===="<<fn6<<endl; 
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

  complex<float> *p01df;
  p01df=alloc1complex(nr);
  complex<float> **p02df;
  p02df=alloc2complex(ltt,nr);
  float **p02dt;
  p02dt=alloc2float(lt,nr);


  complex<float> **a;
  a=alloc2complex(nr,nr);
  complex<float> **at;
  at=alloc2complex(nr,nr);
  complex<float> **ata;
  ata=alloc2complex(nr,nr);
  complex<float> *atb;
  atb=alloc1complex(nr);

   complex<float> *x1;
   x1=alloc1complex(nr);
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

   float *amp1;
   amp1=alloc1float(nr*nr);
   float *amp2;
   amp2=alloc1float(nr);

  int *irbeg;
  irbeg=alloc1int(nr);
  int *irend;
  irend=alloc1int(nr);

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
 
//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         { 
           p02df[ir][it]=(0.0,0.0);
           mf2d[ir][it]=(0.0,0.0);
         }

      for(ir=0;ir<nr;ir++)
        for(ir1=0;ir1<nr;ir1++)
          pp0[ir][ir1]=(0.0,0.0);
       if(is<(l-1)/2)
         {
           for(ir=0;ir<(l-1)/2;ir++)
             {
                irbeg[ir]=0;
                irend[ir]=l-1;
             }
           for(ir=(l-1)/2;ir<nr;ir++)
            {
               irtmp=int(float(is+ir)/2.0);
               irbeg[ir]=irtmp-(l-1)/2;
               irend[ir]=irtmp+(l-1)/2;
            }
         }
       else if(is>nr-(l-1)/2-1&&ir<nr)
         {
           for(ir=0;ir<nr-(l-1)/2-1;ir++)
             {
               irtmp=int(float(is+ir)/2.0);
               irbeg[ir]=irtmp-(l-1)/2;
               irend[ir]=irtmp+(l-1)/2;
             }
           for(ir=nr-(l-1)/2-1;ir<nr;ir++)
            {
                irbeg[ir]=nr-l;
                irend[ir]=nr-1;
            }
         }
        else
         {
            for(ir=0;ir<nr;ir++)
              {
                 irtmp=int(float(is+ir)/2.0);
                 irbeg[ir]=irtmp-(l-1)/2;
                 irend[ir]=irtmp+(l-1)/2;
              }

         }

      swq1.seekg(0,ios::beg);
      swq2.seekg(0,ios::beg);
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
                  wf[ir][ir1]=uswf[ir1][it];
             } 

          for(ir=0;ir<nr;ir++)
           for(ir1=0;ir1<nr;ir1++)
             {
                wf[ir][ir1].real()*=coe1;
                wf[ir][ir1].imag()*=coe1;
             } 


         for(ir=0;ir<nr;ir++)
          for(ir1=0;ir1<nr;ir1++)       
             a[ir][ir1]=(0.0,0.0);
  
        for(ir=0;ir<nr;ir++)       
            {
               irtmp1=irbeg[ir];
               irtmp2=irend[ir];
              for(ir1=irtmp1;ir1<irtmp2;ir1++)
                 a[ir][ir1]=pp0[ir][ir1];
            }

         for(ir=0;ir<nr;ir++)
          for(ir1=0;ir1<nr;ir1++)       
             a[ir][ir1]*=wf[ir][ir1];

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
/*
         for(ir=0;ir<nr;ir++)
            a[ir][ir].real()+=1.0;

         for(ir=0;ir<nr;ir++)
          for(ir1=0;ir1<nr;ir1++)       
           {
             at[ir][ir1].real()=a[ir1][ir].real();
             at[ir][ir1].imag()=-a[ir1][ir].imag();
           }

         complex_matrix_multiply_vector(atb, at, b, nr, nr);

         amp_atb=0.0;

         for(ir=0;ir<nr;ir++)
           amp2[ir]=sqrt(atb[ir].real()*atb[ir].real()+atb[ir].imag()*atb[ir].imag());

         for(ir=0;ir<nr;ir++)
             {
               if(amp2[ir]>amp_atb)
                  amp_atb=amp2[ir];
               else
                  amp_atb=amp_atb;
             }

         cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATB===="<<amp_atb<<endl;

         complex_matrix_multiply_matrix(ata, at, a, nr, nr, nr);

         amp_ata=0.0;   

         for(ir=0;ir<nr;ir++)
          {
            for(ir1=0;ir1<nr;ir1++)
              amp1[ir*nr+ir1]=sqrt(ata[ir][ir1].real()*ata[ir][ir1].real()+ata[ir][ir1].imag()*ata[ir][ir1].imag());
          }

         for(ir=0;ir<nr*nr;ir++)
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

             for(ir=0;ir<nr;ir++)  
               for(ir1=0;ir1<nr;ir1++)
                 {
                    ata[ir][ir1].real()/=amp_ata;
                    ata[ir][ir1].imag()/=amp_ata;
                 }  

            for(ir=0;ir<nr;ir++) 
              ata[ir][ir].real()+=lamda; 

            for(ir=0;ir<nr;ir++)  
               {
                 atb[ir].real()/=amp_atb;
                 atb[ir].imag()/=amp_atb;
               }

            complex_CG(ata, x1, p01df, atb, p11, p2, r1, r2, ap, apcong, nr, err); 

            for(ir=0;ir<nr;ir++) 
              {
                p01df[ir].real()*=amp_atb/amp_ata/coe1;
                p01df[ir].imag()*=amp_atb/amp_ata/coe1;
              }

            for(ir=0;ir<nr;ir++)
               p02df[ir][it]=p01df[ir]; 

            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Done..."<<endl;
          }//6666
          else
           {//7777
             for(ir=0;ir<nr;ir++)
               mf2d[ir][it]=(0.0,0.0);
             cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

             for(ir=0;ir<nr;ir++)
               p02df[ir][it]=(0.0,0.0);

             cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Done..."<<endl;         
           }//7777

*/
        }//8888

       else
        {
           for(ir=0;ir<nr;ir++)
             mf2d[ir][it]=(0.0,0.0);
           cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

           for(ir=0;ir<nr;ir++)
               p02df[ir][it]=(0.0,0.0);

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Done..."<<endl;         
        }

    }//0000

      for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mf2d[ir][it].real()=mf2d[ir][ltt-it].real();
          mf2d[ir][it].imag()=-mf2d[ir][ltt-it].imag();
          p02df[ir][it].real()=p02df[ir][ltt-it].real();
          p02df[ir][it].imag()=-p02df[ir][ltt-it].imag();
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

    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=p02df[ir][it].real();
         in1[it][1]=p02df[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       p02dt[ir][it]=out1[it][0]/ltt;
     }

   for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq6.write((char*)&p02dt[ir][it],sizeof(p02dt[ir][it]));


      cout<<is+1<<" Shot Primary Estimation Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











