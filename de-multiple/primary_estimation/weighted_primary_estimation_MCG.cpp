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
//  int iter_max=300;
  int iter_max=2000;

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

void polygonalFilter(float *f, float *amps, int npoly, int nfft, float dt, float *filter)

#define PIBY2   1.57079632679490
{
        int *intfr;             /* .... integerizations of f            */
        int icount,ifs;         /* loop counting variables              */
        int taper=0;            /* flag counter                         */
        int nf;                 /* number of frequencies (incl Nyq)     */
        int nfm1;               /* nf-1                                 */
        float onfft;            /* reciprocal of nfft                   */
        float df;               /* frequency spacing (from dt)          */


        intfr=alloc1int(npoly);

        nf = nfft/2 + 1;
        nfm1 = nf - 1;
        onfft = 1.0 / nfft;

        /* Compute array of integerized frequencies that define the filter*/
        df = onfft / dt;
        for(ifs=0; ifs < npoly ; ++ifs) {
//                intfr[ifs] = NINT(f[ifs]/df);
                intfr[ifs] = (int)(f[ifs]/df);
                if (intfr[ifs] > nfm1) intfr[ifs] = nfm1;
        }

        /* Build filter, with scale, and taper specified by amps[] values*/
        /* Do low frequency end first*/
        for(icount=0; icount < intfr[0] ; ++icount)
                filter[icount] = amps[0] * onfft;

        /* now do the middle frequencies */
        for(ifs=0 ; ifs<npoly-1 ; ++ifs){
           if(amps[ifs] < amps[ifs+1]) {
                ++taper;
                for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
                    float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
                    float s = sin(c*(icount - intfr[ifs] + 1));
                    float adiff = amps[ifs+1] - amps[ifs];
                    filter[icount] = (amps[ifs] + adiff*s*s) * onfft;
                }
           } else if (amps[ifs] > amps[ifs+1]) {
                ++taper;
                for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
                           float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
                           float s = sin(c*(intfr[ifs+1] - icount + 1));
                           float adiff = amps[ifs] - amps[ifs+1];
                           filter[icount] = (amps[ifs+1] + adiff*s*s) * onfft;
                  }
           } else
                if(!(taper)){
                for(icount=intfr[ifs]; icount <= intfr[ifs+1]; ++icount)
                           filter[icount] = amps[ifs] * onfft;
                } else {
                for(icount=intfr[ifs]+1; icount <= intfr[ifs+1]; ++icount)
                           filter[icount] = amps[ifs] * onfft;
                }
        }

        /* finally do the high frequency end */
        for(icount=intfr[npoly-1]+1; icount<nf; ++icount){
                filter[icount] = amps[npoly-1] * onfft;
        }

}




int main()
{
  char fn1[256],fn2[256],fn111[256],fn222[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fn12[256];
  int ns, nr,lt, nr2, coe1,l,nrl;
  float dt,fmin,fmax,lamda,lfbeg,lfend,hfbeg,hfend,dt1;
  float err=0.00000000000001;
  float ampb;
  float amp_ata,amp_atb;

  int npoly=4;

  ifstream swq;
  swq.open("weighted_primary_estimation_MCG.par");
  swq>>fn1>>fn2>>fn111>>fn222>>fn3>>fn4>>fn9>>fn10>>fn5>>fn6>>fn7>>fn8>>fn11>>fn12>>ns>>nr>>lt>>dt>>fmin>>fmax>>lamda>>coe1>>l>>lfbeg>>lfend>>hfbeg>>hfend;
  swq.close();

  cout<<"Fna of CSG Real is===="<<fn1<<endl; 
  cout<<"Fna of CSG Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of Initial Primary Real is===="<<fn1<<endl; 
  cout<<"Fna of Initial Primary Imaginary is===="<<fn2<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of Initial Multiple Real is===="<<fn9<<endl; 
  cout<<"Fna of Initial Multiple Imaginary is===="<<fn10<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
  cout<<"Fna of Estimated Weights is===="<<fn6<<endl; 
  cout<<"Fna of Estimated multiples is===="<<fn7<<endl; 
  cout<<"Fna of Difference between estimated multiples and true multiples is===="<<fn8<<endl; 
  cout<<"Fna of New Estimated multiples for the forward problem of primary estimation is===="<<fn11<<endl; 
  cout<<"Fna of New Estimated primary is===="<<fn12<<endl; 
  cout<<"No. of Shots is===="<<ns<<endl;
  cout<<"No. of Traces per Shot is===="<<nr<<endl;
  cout<<"lt and dt are===="<<lt<<" , "<<dt<<"ms"<<endl;
  cout<<"minimum and maximum frequency to be calculated are===="<<fmin<<" , "<<fmax<<endl;
  cout<<"MCG Aperture is===="<<l<<endl; 
  cout<<"Frequency Slope is===="<<lfbeg<<"Hz , "<<lfend<<"Hz , "<<hfbeg<<"Hz , "<<hfend<<"Hz"<<endl;


  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1,irtmp,irtmp1,len1,len2; 
  ltt=2*lt;
  dt1=dt/1000;
 
  nrl=nr*l;
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

   len1=10;
   len2=20;

   float *f;
   f=alloc1float(npoly);
   float *amps;
   amps=alloc1float(npoly);
   float *filter;
   filter=alloc1float(lt);
   
   f[0]=lfbeg;
   f[1]=lfend;
   f[2]=hfbeg;
   f[3]=hfend;

   amps[0]=0.0;
   amps[1]=1.0;
   amps[2]=1.0;
   amps[3]=0.0;

   polygonalFilter(f, amps,  npoly, ltt,  dt1, filter);
   
  complex<float> **usf;
  usf=alloc2complex(ltt,nr);
  complex<float> **usf1;
  usf1=alloc2complex(ltt,nr);
  complex<float> **urf;
  urf=alloc2complex(ltt,nr);
 
  complex<float> **mf;
  mf=alloc2complex(ltt,nr);
  complex<float> *af;
  af=alloc1complex(nrl);
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

  int *irbeg;
  irbeg=alloc1int(nr);
  int *irend;
  irend=alloc1int(nr);

  complex<float> *wf1d;
  wf1d=alloc1complex(nrl);
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

  complex<float> **a;
  a=alloc2complex(nrl,nr);
  complex<float> **at;
  at=alloc2complex(nr,nrl);
  complex<float> **ata;
  ata=alloc2complex(nrl,nrl);
  complex<float> *b;
  b=alloc1complex(nr);
  complex<float> *atb;
  atb=alloc1complex(nrl);
  
   complex<float> *x1;
   x1=alloc1complex(nrl);
   complex<float> *p11;
   p11=alloc1complex(nrl);
   complex<float> *p2;
   p2=alloc1complex(nrl);
   complex<float> *r1;
   r1=alloc1complex(nrl);
   complex<float> *r2;
   r2=alloc1complex(nrl);
   complex<float> *ap;
   ap=alloc1complex(nrl);
   complex<float> *apcong;
   apcong=alloc1complex(nrl);

   float *amp1;
   amp1=alloc1float(nrl*nrl);
   float *amp2;
   amp2=alloc1float(nrl);

  complex<float> **anew;
  anew=alloc2complex(nr,nr);
  complex<float> *mf1dnew;
  mf1dnew=alloc1complex(nr);
  complex<float> **mf2dnew;
  mf2dnew=alloc2complex(ltt,nr);
  float **mt2dnew;
  mt2dnew=alloc2float(lt,nr);

  complex<float> **a1;
  a1=alloc2complex(nr,nr);
  complex<float> **at1;
  at1=alloc2complex(nr,nr);
  complex<float> **ata1;
  ata1=alloc2complex(nr,nr);
  complex<float> *b1;
  b1=alloc1complex(nr);
  complex<float> *atb1;
  atb1=alloc1complex(nr);
  
   complex<float> *x11;
   x11=alloc1complex(nr);
   complex<float> *p111;
   p111=alloc1complex(nr);
   complex<float> *p21;
   p21=alloc1complex(nr);
   complex<float> *r11;
   r11=alloc1complex(nr);
   complex<float> *r21;
   r21=alloc1complex(nr);
   complex<float> *ap1;
   ap1=alloc1complex(nr);
   complex<float> *apcong1;
   apcong1=alloc1complex(nr);

   float amp3max,amp4max,amp5max,amp6max;

   float *amp3;
   amp3=alloc1float(nr*nr);
   float *amp4;
   amp4=alloc1float(nr);
   float *amp5;
   amp5=alloc1float(nr*nr);
   float *amp6;
   amp6=alloc1float(nr);

  complex<float> **pf2d;
  pf2d=alloc2complex(ltt,nr);
  complex<float> *pf1d;
  pf1d=alloc1complex(nr);

  float **pt2d;
  pt2d=alloc2float(lt,nr);

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

  ifstream swq111;
  swq111.open(fn111,ios::binary);
  if(!swq111)
       cout<<"cannot open "<<fn111<<endl; 

  ifstream swq222;
  swq222.open(fn222,ios::binary);
  if(!swq222)
       cout<<"cannot open "<<fn222<<endl;

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
 
  ofstream swq11;
  swq11.open(fn11,ios::binary);
  if(!swq11)
       cout<<"cannot open "<<fn11<<endl;
 
  ofstream swq12;
  swq12.open(fn12,ios::binary);
  if(!swq12)
       cout<<"cannot open "<<fn11<<endl;
 
//  for(is=0;is<ns;is++)
  for(is=99;is<100;is++)
    {
       for(ir=0;ir<nr2;ir++)
        for(it=0;it<ltt;it++)
          wf2d[ir][it]=(0.0,0.0);

       for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           pf2d[ir][it]=(0.0,0.0);    
           mcalf2d[ir][it]=(0.0,0.0);    
           mf2dnew[ir][it]=(0.0,0.0);    
         }

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

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
           mf[ir][it]=(0.0,0.0);

      swq1.seekg(0,ios::beg);
      swq2.seekg(0,ios::beg);
      swq111.seekg(0,ios::beg);
      swq222.seekg(0,ios::beg);
      swq9.seekg(0,ios::beg);
      swq10.seekg(0,ios::beg);

      for(ir=0;ir<is;ir++)
        {
           swq1.seekg(nr*ltt*4,ios::cur);
           swq2.seekg(nr*ltt*4,ios::cur);
           swq111.seekg(nr*ltt*4,ios::cur);
           swq222.seekg(nr*ltt*4,ios::cur);
           swq9.seekg(nr*ltt*4,ios::cur);
           swq10.seekg(nr*ltt*4,ios::cur);
        }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real()));  //recorded csg 
           swq2.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag())); 
           swq111.read((char*)&usf1[ir][it].real(),sizeof(usf1[ir][it].real()));  //primary 
           swq222.read((char*)&usf1[ir][it].imag(),sizeof(usf1[ir][it].imag())); 
           swq9.read((char*)&mftrue[ir][it].real(),sizeof(mftrue[ir][it].real())); 
           swq10.read((char*)&mftrue[ir][it].imag(),sizeof(mftrue[ir][it].imag())); 
         }

      for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
         {
            usf[ir][it].real()*=coe1;
            usf[ir][it].imag()*=coe1;
            usf1[ir][it].real()*=coe1;
            usf1[ir][it].imag()*=coe1;
            mftrue[ir][it].real()*=coe1;
            mftrue[ir][it].imag()*=coe1;
         } 

      for(it=0;it<ltt/2;it++)
         {
           for(ir=0;ir<nr2;ir++)
            wf2d[ir][it]=(0.0,0.0);

           for(ir=0;ir<nr;ir++)
            {  
              pf2d[ir][it]=(0.0,0.0);
              mcalf2d[ir][it]=(0.0,0.0);
              mf2dnew[ir][it]=(0.0,0.0);
            }
         }

      for(it=ifmin;it<ifmax;it++)
        {//0000
          for(ir=0;ir<nrl;ir++)
             {
               af[ir].real()=1.0;
               af[ir].imag()=0.0;
             }

          for(ir=0;ir<nr;ir++)
            { 
              p0[ir]=usf1[ir][it];
              b1[ir]=usf[ir][it];
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

          for(ir=0;ir<nr;ir++)       
            for(ir1=0;ir1<nrl;ir1++)
               a[ir][ir1]=(0.0,0.0);

          for(ir=0;ir<nr;ir++)       
            {
               irtmp1=irbeg[ir];
               for(ir1=ir*l;ir1<(ir+1)*l;ir1++)
                  a[ir][ir1]=pp0[ir][ir*nr+irtmp1+ir1-ir*l];    
            } 

         complex_matrix_multiply_vector(mf1d, a, af, nr, nrl);

          for(ir=0;ir<nr;ir++)
            mf[ir][it]=mf1d[ir];

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

         for(ir=0;ir<nrl;ir++)
          for(ir1=0;ir1<nr;ir1++)       
           {
             at[ir][ir1].real()=a[ir1][ir].real();
             at[ir][ir1].imag()=-a[ir1][ir].imag();
           }

         complex_matrix_multiply_vector(atb, at, b, nrl, nr);

         amp_atb=0.0;

         for(ir=0;ir<nrl;ir++)
           amp2[ir]=sqrt(atb[ir].real()*atb[ir].real()+atb[ir].imag()*atb[ir].imag());

         for(ir=0;ir<nrl;ir++)
             {
               if(amp2[ir]>amp_atb)
                  amp_atb=amp2[ir];
               else
                  amp_atb=amp_atb;
             }

         cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATB===="<<amp_atb<<endl;

         complex_matrix_multiply_matrix(ata, at, a, nrl, nr, nrl);

         amp_ata=0.0;   

         for(ir=0;ir<nrl;ir++)
          {
            for(ir1=0;ir1<nrl;ir1++)
              amp1[ir*nrl+ir1]=sqrt(ata[ir][ir1].real()*ata[ir][ir1].real()+ata[ir][ir1].imag()*ata[ir][ir1].imag());
          }

         for(ir=0;ir<nrl*nrl;ir++)
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

             for(ir=0;ir<nrl;ir++)  
               for(ir1=0;ir1<nrl;ir1++)
                 {
                    ata[ir][ir1].real()/=amp_ata;
                    ata[ir][ir1].imag()/=amp_ata;
                 }  

            for(ir=0;ir<nrl;ir++) 
              ata[ir][ir].real()+=lamda; 

            for(ir=0;ir<nrl;ir++)  
               {
                 atb[ir].real()/=amp_atb;
                 atb[ir].imag()/=amp_atb;
               }

            complex_CG(ata, x1, wf1d, atb, p11, p2, r1, r2, ap, apcong, nrl, err); 

            for(ir=0;ir<nrl;ir++) 
              {
                wf1d[ir].real()*=amp_atb/amp_ata/coe1;
                wf1d[ir].imag()*=amp_atb/amp_ata/coe1;
              }

            complex_matrix_multiply_vector(mcalf1d, a, wf1d, nr, nrl);

            for(ir=0;ir<nr;ir++)
              mcalf2d[ir][it]=mcalf1d[ir];

            for(ir=0;ir<nr;ir++)
               {
                  irtmp1=irbeg[ir];
                  for(ir1=ir*nr+irtmp1;ir1<ir*nr+irtmp1+l;ir1++)
                    wf2d[ir1][it]=wf1d[ir*l+ir1-ir*nr-irtmp1]; 
               }

            for(ir=0;ir<nr;ir++)
               {
                 for(ir1=0;ir1<nr;ir1++)
                    anew[ir][ir1]=p[ir][ir1]*wf2d[ir*nr+ir1][it];
               }


            complex_matrix_multiply_vector(mf1dnew, anew, p0, nr, nr);

            for(ir=0;ir<nr;ir++)
              mf2dnew[ir][it]=mf1dnew[ir];

            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;
            
            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Begins..."<<endl;

            for(ir=0;ir<nr;ir++)
              for(ir1=0;ir1<nr;ir1++)
                 a1[ir][ir1]=anew[ir][ir1];               

            for(ir=0;ir<nr;ir++)
              for(ir1=0;ir1<nr;ir1++)
                {
                    a1[ir][ir1].real()*=coe1;
                    a1[ir][ir1].imag()*=coe1;
                }

            for(ir=0;ir<nr;ir++)
                a1[ir][ir].real()+=1.0;
/*
            if(it!=0)
              {
                 cout<<"    "<<is+1<<" Shot, "<<it+1<<"  Frequency Slice ANEW is===="<<endl;
                 for(ir=0;ir<nr;ir++)
                    {
                       cout<<"===="<<ir+1<<endl;
                       for(ir1=0;ir1<nr;ir1++)
                          cout<<a1[ir][ir1]<<"     ";
                       cout<<endl;
                    }

                 cout<<"    "<<is+1<<" Shot, "<<it+1<<"  Frequency Slice B is===="<<endl;
                 for(ir=0;ir<nr;ir++)                    
                   cout<<b1[ir]<<"     "; 
                 cout<<endl;                  
              }
*/

            for(ir=0;ir<nr;ir++)
              for(ir1=0;ir1<nr;ir1++)       
               {
                 at1[ir][ir1].real()=a1[ir1][ir].real();
                 at1[ir][ir1].imag()=-a1[ir1][ir].imag();
               }

            complex_matrix_multiply_matrix(ata1, at1, a1, nr, nr, nr);

            for(ir=0;ir<nr;ir++)
             {
               for(ir1=0;ir1<nr;ir1++)
                  amp5[ir*nr+ir1]=sqrt(ata1[ir][ir1].real()*ata1[ir][ir1].real()+ata1[ir][ir1].imag()*ata1[ir][ir1].imag()); 
             } 

            amp5max=0.0;
            for(ir=0;ir<nr*nr;ir++)
             {
              if(amp5[ir]>amp5max)
                amp5max=amp5[ir];
              else
                amp5max=amp5max;
             }

            for(ir=0;ir<nr;ir++)
              for(ir1=0;ir1<nr;ir1++)
                {
                  ata1[ir][ir1].real()/=amp5max;
                  ata1[ir][ir1].imag()/=amp5max;
                }

            complex_matrix_multiply_vector(atb1, at1, b1, nr, nr);

            for(ir=0;ir<nr;ir++)
             amp6[ir]=sqrt(atb1[ir].real()*atb1[ir].real()+atb1[ir].imag()*atb1[ir].imag()); 

            amp6max=0.0;
            for(ir=0;ir<nr;ir++)
             {
              if(amp6[ir]>amp6max)
                amp6max=amp6[ir];
              else
                amp6max=amp6max;
             }
            for(ir=0;ir<nr;ir++)
                {
                  atb1[ir].real()/=amp6max;
                  atb1[ir].imag()/=amp6max;
                }
/*
            if(it==51)
              {
                 cout<<"    "<<is+1<<" Shot, "<<it+1<<"  Frequency Slice ATA is===="<<endl;
                 for(ir=0;ir<nr;ir++)
                    {
                       cout<<"===="<<ir+1<<endl;
                       for(ir1=0;ir1<nr;ir1++)
                          cout<<ata1[ir][ir1]<<"     ";
                       cout<<endl;
                    }

                 cout<<"    "<<is+1<<" Shot, "<<it+1<<"  Frequency Slice ATB is===="<<endl;
                 for(ir=0;ir<nr;ir++)                    
                   cout<<atb1[ir]<<"     "; 
                 cout<<endl;                  
              }
*/

            complex_CG(ata1, x11, pf1d , atb1, p111, p21, r11, r21, ap1, apcong1, nr, err); 

            for(ir=0;ir<nr;ir++)
               pf2d[ir][it]=pf1d[ir];

            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Done..."<<endl;

          }//6666

          else
           {//7777
             for(ir=0;ir<nr;ir++)
               mf[ir][it]=(0.0,0.0);
             cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

             for(ir=0;ir<nr2;ir++)
              wf2d[ir][it]=(0.0,0.0);

             for(ir=0;ir<nr;ir++)
              {  
                pf2d[ir][it]=(0.0,0.0);
                mcalf2d[ir][it]=(0.0,0.0);
                mf2dnew[ir][it]=(0.0,0.0);
              }

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
              {  
                pf2d[ir][it]=(0.0,0.0);
                mcalf2d[ir][it]=(0.0,0.0);
                mf2dnew[ir][it]=(0.0,0.0);
              }

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;         
        }

    }//0000



      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
         {
            pf2d[ir][it].real()*=filter[it];
            pf2d[ir][it].imag()*=filter[it];
         }
    


      for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mf[ir][it].real()=mf[ir][ltt-it].real();
          mf[ir][it].imag()=-mf[ir][ltt-it].imag();
          mcalf2d[ir][it].real()=mcalf2d[ir][ltt-it].real();
          mcalf2d[ir][it].imag()=-mcalf2d[ir][ltt-it].imag();
          mf2dnew[ir][it].real()=mf2dnew[ir][ltt-it].real();
          mf2dnew[ir][it].imag()=-mf2dnew[ir][ltt-it].imag();
          pf2d[ir][it].real()=pf2d[ir][ltt-it].real();
          pf2d[ir][it].imag()=-pf2d[ir][ltt-it].imag();
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
         in1[it][0]=mcalf2d[ir][it].real()*coe1;
         in1[it][1]=mcalf2d[ir][it].imag()*coe1;
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
        {
           mftrue[ir][it].real()/=coe1;
           mftrue[ir][it].imag()/=coe1;
        }

    for(ir=0;ir<nr;ir++)
     for(it=0;it<ltt;it++)
       diff2d[ir][it]=mcalf2d[ir][it]-mftrue[ir][it];

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

    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mf2dnew[ir][it].real()*coe1;
         in1[it][1]=mf2dnew[ir][it].imag()*coe1;
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mt2dnew[ir][it]=out1[it][0]/ltt;
     }

    for(ir=0;ir<nr;ir++)
     for(it=0;it<lt;it++)
       swq11.write((char*)&mt2dnew[ir][it],sizeof(mt2dnew[ir][it]));

    for(ir=0;ir<nr;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=pf2d[ir][it].real();
         in1[it][1]=pf2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       pt2d[ir][it]=out1[it][0]/ltt;
     }

    for(ir=0;ir<nr;ir++)
     for(it=0;it<lt;it++)
       swq12.write((char*)&pt2d[ir][it],sizeof(pt2d[ir][it]));

      cout<<is+1<<" Shot Inverse Estimation Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











