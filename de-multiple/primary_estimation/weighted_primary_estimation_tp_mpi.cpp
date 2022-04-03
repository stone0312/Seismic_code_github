using namespace std;
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"
#include "mpi.h"

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

int main(int argc,char **argv)
{
  char fn1[256],fn2[256],fn13[256],fn14[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fn12[256],fn15[256],fn16[256],fn17[256],fn18[256],fn19[256],fn20[256],fn21[256],fn22[256],fn23[256],fn24[256];

  int ns, nr,lt, nr2, coe1,l,nrl,np;
  float dt,fmin,fmax,lamda,lfbeg,lfend,hfbeg,hfend,dt1;
  float err=0.00000000000001;
  float ampb;
  float amp_ata,amp_atb;

  int npp,myid;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&npp);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  int npoly=4;

  ifstream swq;
  swq.open("weighted_primary_estimation_tp_mpi.par");
  swq>>fn1>>fn2>>fn13>>fn14>>fn3>>fn4>>fn9>>fn10>>fn5>>fn6>>fn7>>fn8>>fn11>>fn12>>fn15>>fn16>>fn17>>fn18>>fn19>>fn20>>fn21>>fn22>>fn23>>fn24>>ns>>nr>>lt>>np>>dt>>fmin>>fmax>>lamda>>coe1>>l>>lfbeg>>lfend>>hfbeg>>hfend;
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

  int ltt,ifmin,ifmax,is,ir,ir1,ir2,it,it1,irtmp,irtmp1,len1,len2,ip,nrnp,lnrnp,nnrnp; 
  ltt=2*lt;
  dt1=dt/1000;
 
  nrl=nr*l;
  nr2=nr*nr;
  nrnp=nr*np;
  lnrnp=l*nr*np;
  nnrnp=nr*nr*np;

  float *omega;
  omega=alloc1float(ltt);
  for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
  for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));
 
  ifmin=int(fmin*dt*ltt/1000);
  ifmax=int(fmax*dt*ltt/1000)+1;

  cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

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
   
  complex<float> **usfp;
  usfp=alloc2complex(ltt,nrnp);

  complex<float> **usfp1;
  usfp1=alloc2complex(ltt,nrnp);

  complex<float> **urfp;
  urfp=alloc2complex(ltt,nrnp);

  complex<float> **urfp1;
  urfp1=alloc2complex(np,nr);

  complex<float> **urfp2;
  urfp2=alloc2complex(np,nr);
 
  complex<float> **mfp;
  mfp=alloc2complex(ltt,nrnp);
  complex<float> *af;
  af=alloc1complex(lnrnp);
  complex<float> **pp0;
  pp0=alloc2complex(nrnp,nrnp);

  complex<float> *p0;
  p0=alloc1complex(nrnp);
  complex<float> **p;
  p=alloc2complex(nrnp,nrnp);

  complex<float> *mfp1d;
  mfp1d=alloc1complex(nrnp);
  float **mtp;
  mtp=alloc2float(lt,nrnp);

  int *irbeg;
  irbeg=alloc1int(nr);
  int *irend;
  irend=alloc1int(nr);

//arrays for weight estimation
  complex<float> *wpf1d;
  wpf1d=alloc1complex(lnrnp);
  complex<float> **wpf2d;
  wpf2d=alloc2complex(ltt,nnrnp);
  complex<float> **wpf2dtmp;
  wpf2dtmp=alloc2complex(ltt,nnrnp);

  complex<float> **mfptrue;
  mfptrue=alloc2complex(ltt,nrnp);

  complex<float> *mcalfp1d;
  mcalfp1d=alloc1complex(nrnp);

  complex<float> **a;
  a=alloc2complex(lnrnp,nrnp);
  complex<float> **at;
  at=alloc2complex(nrnp,lnrnp);
  complex<float> **ata;
  ata=alloc2complex(lnrnp,lnrnp);
  complex<float> *b;
  b=alloc1complex(nrnp);
  complex<float> *atb;
  atb=alloc1complex(lnrnp);
  
   complex<float> *x1;
   x1=alloc1complex(lnrnp);
   complex<float> *p11;
   p11=alloc1complex(lnrnp);
   complex<float> *p2;
   p2=alloc1complex(lnrnp);
   complex<float> *r1;
   r1=alloc1complex(lnrnp);
   complex<float> *r2;
   r2=alloc1complex(lnrnp);
   complex<float> *ap;
   ap=alloc1complex(lnrnp);
   complex<float> *apcong;
   apcong=alloc1complex(lnrnp);

   float *amp1;
   amp1=alloc1float(lnrnp*lnrnp);
   float *amp2;
   amp2=alloc1float(lnrnp);

  complex<float> **anew;
  anew=alloc2complex(nrnp,nrnp);
  complex<float> *mfp1dnew;
  mfp1dnew=alloc1complex(nrnp);

  complex<float> **a1;
  a1=alloc2complex(nrnp,nrnp);
  complex<float> **at1;
  at1=alloc2complex(nrnp,nrnp);
  complex<float> **ata1;
  ata1=alloc2complex(nrnp,nrnp);
  complex<float> *b1;
  b1=alloc1complex(nrnp);
  complex<float> *atb1;
  atb1=alloc1complex(nrnp);
  
   complex<float> *x11;
   x11=alloc1complex(nrnp);
   complex<float> *p111;
   p111=alloc1complex(nrnp);
   complex<float> *p21;
   p21=alloc1complex(nrnp);
   complex<float> *r11;
   r11=alloc1complex(nrnp);
   complex<float> *r21;
   r21=alloc1complex(nrnp);
   complex<float> *ap1;
   ap1=alloc1complex(nrnp);
   complex<float> *apcong1;
   apcong1=alloc1complex(nrnp);

   float amp3max,amp4max,amp5max,amp6max;

   float *amp3;
   amp3=alloc1float(nrnp*nrnp);
   float *amp4;
   amp4=alloc1float(nrnp);
   float *amp5;
   amp5=alloc1float(nrnp*nrnp);
   float *amp6;
   amp6=alloc1float(nrnp);

  complex<float> *ppf1d;
  ppf1d=alloc1complex(nrnp);

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

  ifstream swq13;
  swq13.open(fn13,ios::binary);
  if(!swq13)
       cout<<"cannot open "<<fn13<<endl; 

  ifstream swq14;
  swq14.open(fn14,ios::binary);
  if(!swq14)
       cout<<"cannot open "<<fn14<<endl;

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

  fstream swq5;
  swq5.open(fn5,ios::binary|ios::out);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;

  fstream swq6;
  swq6.open(fn6,ios::binary|ios::out);
  if(!swq6)
       cout<<"cannot open "<<fn6<<endl;

  fstream swq7;
  swq7.open(fn7,ios::binary|ios::out);
  if(!swq7)
       cout<<"cannot open "<<fn7<<endl;

  fstream swq8;
  swq8.open(fn8,ios::binary|ios::out);
  if(!swq8)
       cout<<"cannot open "<<fn8<<endl;
 
  fstream swq11;
  swq11.open(fn11,ios::binary|ios::out);
  if(!swq11)
       cout<<"cannot open "<<fn11<<endl;
 
  fstream swq12;
  swq12.open(fn12,ios::binary|ios::out);
  if(!swq12)
       cout<<"cannot open "<<fn11<<endl;
 
  fstream swq155; //write calculated mfp1d real
  swq155.open(fn15,ios::binary|ios::out);
  if(!swq155)
       cout<<"cannot open "<<fn15<<endl;
 
  fstream swq166;  // write calculated mfp1d imag
  swq166.open(fn16,ios::binary|ios::out);
  if(!swq166)
       cout<<"cannot open "<<fn16<<endl;
 
  fstream swq177;   //write calculated wpf1d real
  swq177.open(fn17,ios::binary|ios::out);
  if(!swq177)
       cout<<"cannot open "<<fn17<<endl;
 
  fstream swq188;   //write calculated wpf1d imag
  swq188.open(fn18,ios::binary|ios::out);
  if(!swq188)
       cout<<"cannot open "<<fn18<<endl;
 
  fstream swq199; //write calculated mcalfp1d real
  swq199.open(fn19,ios::binary|ios::out);
  if(!swq199)
       cout<<"cannot open "<<fn19<<endl;
 
  fstream swq200;  //write calculated mcalfp1d imag
  swq200.open(fn20,ios::binary|ios::out);
  if(!swq200)
       cout<<"cannot open "<<fn20<<endl;
 
  fstream swq211;   //write calculated ppf1d real
  swq211.open(fn21,ios::binary|ios::out);
  if(!swq211)
       cout<<"cannot open "<<fn21<<endl;
 
  fstream swq222;   //write calculated ppf1d imag
  swq222.open(fn22,ios::binary|ios::out);
  if(!swq222)
       cout<<"cannot open "<<fn22<<endl;
 
  fstream swq233;   //write calculated mcalfp1dnew real
  swq233.open(fn23,ios::binary|ios::out);
  if(!swq233)
       cout<<"cannot open "<<fn23<<endl;
 
  fstream swq244;   //write calculated mcalfp1dnew imag
  swq244.open(fn24,ios::binary|ios::out);
  if(!swq244)
       cout<<"cannot open "<<fn24<<endl;
 
//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
       for(ir=0;ir<nnrnp;ir++)
        for(it=0;it<ltt;it++)
          wpf2dtmp[ir][it]=(0.0,0.0);

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
      swq13.seekg(0,ios::beg);
      swq14.seekg(0,ios::beg);
      swq9.seekg(0,ios::beg);
      swq10.seekg(0,ios::beg);

      for(ir=0;ir<is;ir++)
        {
           swq1.seekg(nrnp*ltt*4,ios::cur);
           swq2.seekg(nrnp*ltt*4,ios::cur);
           swq13.seekg(nrnp*ltt*4,ios::cur);
           swq14.seekg(nrnp*ltt*4,ios::cur);
           swq9.seekg(nrnp*ltt*4,ios::cur);
           swq10.seekg(nrnp*ltt*4,ios::cur);
        }

      for(ir=0;ir<nrnp;ir++)
        for(it=0;it<ltt;it++)
         {
           swq1.read((char*)&usfp[ir][it].real(),sizeof(usfp[ir][it].real()));  //recorded csg 
           swq2.read((char*)&usfp[ir][it].imag(),sizeof(usfp[ir][it].imag())); 
           swq13.read((char*)&usfp1[ir][it].real(),sizeof(usfp1[ir][it].real()));  //primary 
           swq14.read((char*)&usfp1[ir][it].imag(),sizeof(usfp1[ir][it].imag())); 
           swq9.read((char*)&mfptrue[ir][it].real(),sizeof(mfptrue[ir][it].real())); 
           swq10.read((char*)&mfptrue[ir][it].imag(),sizeof(mfptrue[ir][it].imag())); 
         }

      for(ir=0;ir<nrnp;ir++)
        for(it=0;it<ltt;it++)
         {
            usfp[ir][it].real()*=coe1;
            usfp[ir][it].imag()*=coe1;
            usfp1[ir][it].real()*=coe1;
            usfp1[ir][it].imag()*=coe1;
            mfptrue[ir][it].real()*=coe1;
            mfptrue[ir][it].imag()*=coe1;
         } 

      for(it=ifmin+myid;it<ifmax;it+=npp)
        {//0000

          for(ir=0;ir<lnrnp;ir++)
             {
               af[ir].real()=1.0;
               af[ir].imag()=0.0;
             }

          for(ir=0;ir<nrnp;ir++)
            { 
              p0[ir]=usfp1[ir][it];
              b1[ir]=usfp[ir][it];
              b[ir]=mfptrue[ir][it];
            }

          for(ir1=0;ir1<nrnp;ir1++)
            for(ir2=0;ir2<nrnp;ir2++)
              {
                p[ir1][ir2]=(0.0,0.0);
                pp0[ir1][ir2]=(0.0,0.0);
              }

          for(ir1=0;ir1<nrnp;ir1++)
            for(ir2=0;ir2<lnrnp;ir2++)
                a[ir1][ir2]=(0.0,0.0);
 
         ampb=0.0;
          for(ir=0;ir<nrnp;ir++)
           ampb+=b[ir].real()*b[ir].real()+b[ir].imag()*b[ir].imag();

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude===="<<ampb<<endl;

          if(ampb!=0)
          {//8888
           swq3.seekg(0,ios::beg);
           swq4.seekg(0,ios::beg);

           for(ir=0;ir<nr;ir++)
              {
                for(ir1=0;ir1<nrnp;ir1++)
                 for(it1=0;it1<ltt;it1++)
                  {
                    swq3.read((char*)&urfp[ir1][it1].real(),sizeof(urfp[ir1][it1].real()));
                    swq4.read((char*)&urfp[ir1][it1].imag(),sizeof(urfp[ir1][it1].imag()));
                  }

                for(ir1=0;ir1<nr;ir1++)
                 {
                   for(ip=0;ip<np;ip++)
                     urfp1[ir1][ip]=urfp[ir1*np+ip][it];
                 }
                for(ir1=0;ir1<nr;ir1++)
                 {
                   for(ip=0;ip<np;ip++)
                     urfp2[ir1][ip]=urfp1[ir1][np-ip];
                 }

                for(ir1=ir*np;ir1<(ir+1)*np;ir1++)
                 {
                    for(ip=0;ip<np;ip++)
                      {
                         for(ir2=ip*nr;ir2<(ip+1)*nr;ir2++)
                           p[ir1][ir2]=urfp2[ir2-ip*nr][ir1-ir*np];
                      }
                 }
              }

          for(ir=0;ir<nr;ir++)
            {
              for(ir1=ir*np;ir1<(ir+1)*np;ir1++)
                {
                  for(ir2=(ir1-ir*np)*nr;ir2<(ir1-ir*np+1)*nr;ir2++)
                    pp0[ir1][ir2]=p[ir1][ir2]*p0[ir2];
                }
            }

//form matrix A nrnp*lnrnp
          for(ir=0;ir<nrnp;ir++)       
            for(ir1=0;ir1<lnrnp;ir1++)
               a[ir][ir1]=(0.0,0.0);

         for(ir=0;ir<nr;ir++)
            {
               irtmp1=irbeg[ir];
               for(ir1=ir*np;ir1<(ir+1)*np;ir1++)
                 {
                    for(ir2=ir1*l;ir2<(ir1+1)*l;ir2++)
                      a[ir1][ir2]=pp0[ir1][(ir1-ir*np)*nr+irtmp1+ir2-ir1*l];
                 }
            }

         complex_matrix_multiply_vector(mfp1d, a, af, nrnp, lnrnp);

         swq155.seekg(0,ios::beg);
         swq166.seekg(0,ios::beg);

         for(it1=ifmin;it1<it;it1++)
            {
               swq155.seekg(nrnp*4,ios::cur);
               swq166.seekg(nrnp*4,ios::cur);
            }

        for(ir=0;ir<nrnp;ir++)
           {
             swq155.write((char*)&mfp1d[ir].real(),sizeof(mfp1d[ir].real()));
             swq166.write((char*)&mfp1d[ir].imag(),sizeof(mfp1d[ir].imag()));
           } 

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

         for(ir=0;ir<lnrnp;ir++)
          for(ir1=0;ir1<nrnp;ir1++)       
           {
             at[ir][ir1].real()=a[ir1][ir].real();
             at[ir][ir1].imag()=-a[ir1][ir].imag();
           }

         complex_matrix_multiply_vector(atb, at, b, lnrnp, nrnp);

         amp_atb=0.0;

         for(ir=0;ir<lnrnp;ir++)
           amp2[ir]=sqrt(atb[ir].real()*atb[ir].real()+atb[ir].imag()*atb[ir].imag());

         for(ir=0;ir<lnrnp;ir++)
             {
               if(amp2[ir]>amp_atb)
                  amp_atb=amp2[ir];
               else
                  amp_atb=amp_atb;
             }

         cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATB===="<<amp_atb<<endl;

         complex_matrix_multiply_matrix(ata, at, a, lnrnp, nrnp, lnrnp);

         amp_ata=0.0;   

         for(ir=0;ir<lnrnp;ir++)
          {
            for(ir1=0;ir1<lnrnp;ir1++)
              amp1[ir*nrl+ir1]=sqrt(ata[ir][ir1].real()*ata[ir][ir1].real()+ata[ir][ir1].imag()*ata[ir][ir1].imag());
          }

         for(ir=0;ir<lnrnp*lnrnp;ir++)
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

             for(ir=0;ir<lnrnp;ir++)  
               for(ir1=0;ir1<lnrnp;ir1++)
                 {
                    ata[ir][ir1].real()/=amp_ata;
                    ata[ir][ir1].imag()/=amp_ata;
                 }  

            for(ir=0;ir<lnrnp;ir++) 
              ata[ir][ir].real()+=lamda; 

            for(ir=0;ir<lnrnp;ir++)  
               {
                 atb[ir].real()/=amp_atb;
                 atb[ir].imag()/=amp_atb;
               }

            complex_CG(ata, x1, wpf1d, atb, p11, p2, r1, r2, ap, apcong, lnrnp, err); 

            for(ir=0;ir<lnrnp;ir++) 
              {
                wpf1d[ir].real()*=amp_atb/amp_ata/coe1;
                wpf1d[ir].imag()*=amp_atb/amp_ata/coe1;
              }

         swq177.seekg(0,ios::beg);
         swq188.seekg(0,ios::beg);

         for(it1=ifmin;it1<it;it1++)
            {
               swq177.seekg(lnrnp*4,ios::cur);
               swq188.seekg(lnrnp*4,ios::cur);
            }

        for(ir=0;ir<lnrnp;ir++)
           {
             swq177.write((char*)&wpf1d[ir].real(),sizeof(wpf1d[ir].real()));
             swq188.write((char*)&wpf1d[ir].imag(),sizeof(wpf1d[ir].imag()));
           } 

         complex_matrix_multiply_vector(mcalfp1d, a, wpf1d, nrnp, lnrnp);

         swq199.seekg(0,ios::beg);
         swq200.seekg(0,ios::beg);

         for(it1=ifmin;it1<it;it1++)
            {
               swq199.seekg(nrnp*4,ios::cur);
               swq200.seekg(nrnp*4,ios::cur);
            }

        for(ir=0;ir<nrnp;ir++)
           {
             swq199.write((char*)&mcalfp1d[ir].real(),sizeof(mcalfp1d[ir].real()));
             swq200.write((char*)&mcalfp1d[ir].imag(),sizeof(mcalfp1d[ir].imag()));
           } 

            for(ir=0;ir<nr;ir++)
              {
                irtmp1=irbeg[ir];
                for(ir1=0;ir1<np;ir1++)
                 {
                   for(ir2=ir*nr*np+ir1*nr+irtmp1;ir2<ir*nr*np+ir1*nr+irtmp1+l;ir2++)
                      wpf2dtmp[ir2][it]=wpf1d[ir*l*np+ir1*l+ir2-ir*nr*np-ir1*nr-irtmp1];
                 }
              }

            for(ir=0;ir<nr;ir++)
               {
                 for(ir1=ir*np;ir1<(ir+1)*np;ir1++)
                   {
                    for(ir2=(ir1-ir*np)*nr;ir2<(ir1-ir*np+1)*nr;ir2++) 
                      anew[ir1][ir2]=p[ir1][ir2]*wpf2dtmp[ir1*nr+ir2-(ir1-ir*np)*nr][it];
                   }
               }

            complex_matrix_multiply_vector(mfp1dnew, anew, p0, nrnp, nrnp);

         swq233.seekg(0,ios::beg);
         swq244.seekg(0,ios::beg);

         for(it1=ifmin;it1<it;it1++)
            {
               swq233.seekg(nrnp*4,ios::cur);
               swq244.seekg(nrnp*4,ios::cur);
            }

        for(ir=0;ir<nrnp;ir++)
           {
             swq233.write((char*)&mfp1dnew[ir].real(),sizeof(mfp1dnew[ir].real()));
             swq244.write((char*)&mfp1dnew[ir].imag(),sizeof(mfp1dnew[ir].imag()));
           } 

            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;
            
            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Begins..."<<endl;

            for(ir=0;ir<nrnp;ir++)
              for(ir1=0;ir1<nrnp;ir1++)
                 a1[ir][ir1]=anew[ir][ir1];               

            for(ir=0;ir<nrnp;ir++)
              for(ir1=0;ir1<nrnp;ir1++)
                {
                    a1[ir][ir1].real()*=coe1;
                    a1[ir][ir1].imag()*=coe1;
                }

            for(ir=0;ir<nrnp;ir++)
                a1[ir][ir].real()+=1.0;

            for(ir=0;ir<nrnp;ir++)
              for(ir1=0;ir1<nrnp;ir1++)       
               {
                 at1[ir][ir1].real()=a1[ir1][ir].real();
                 at1[ir][ir1].imag()=-a1[ir1][ir].imag();
               }

            complex_matrix_multiply_matrix(ata1, at1, a1, nrnp, nrnp, nrnp);

            for(ir=0;ir<nrnp;ir++)
             {
               for(ir1=0;ir1<nrnp;ir1++)
                  amp5[ir*nr+ir1]=sqrt(ata1[ir][ir1].real()*ata1[ir][ir1].real()+ata1[ir][ir1].imag()*ata1[ir][ir1].imag()); 
             } 

            amp5max=0.0;
            for(ir=0;ir<nrnp*nrnp;ir++)
             {
              if(amp5[ir]>amp5max)
                amp5max=amp5[ir];
              else
                amp5max=amp5max;
             }

            for(ir=0;ir<nrnp;ir++)
              for(ir1=0;ir1<nrnp;ir1++)
                {
                  ata1[ir][ir1].real()/=amp5max;
                  ata1[ir][ir1].imag()/=amp5max;
                }

            complex_matrix_multiply_vector(atb1, at1, b1, nrnp, nrnp);

            for(ir=0;ir<nrnp;ir++)
             amp6[ir]=sqrt(atb1[ir].real()*atb1[ir].real()+atb1[ir].imag()*atb1[ir].imag()); 

            amp6max=0.0;
            for(ir=0;ir<nrnp;ir++)
             {
              if(amp6[ir]>amp6max)
                amp6max=amp6[ir];
              else
                amp6max=amp6max;
             }
            for(ir=0;ir<nrnp;ir++)
                {
                  atb1[ir].real()/=amp6max;
                  atb1[ir].imag()/=amp6max;
                }

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATA===="<<amp5max<<endl;
          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice, Amplitude ATb===="<<amp6max<<endl;

            complex_CG(ata1, x11, ppf1d , atb1, p111, p21, r11, r21, ap1, apcong1, nrnp, err); 

         swq211.seekg(0,ios::beg);
         swq222.seekg(0,ios::beg);

         for(it1=ifmin;it1<it;it1++)
            {
               swq211.seekg(nrnp*4,ios::cur);
               swq222.seekg(nrnp*4,ios::cur);
            }

        for(ir=0;ir<nrnp;ir++)
           {
             swq211.write((char*)&ppf1d[ir].real(),sizeof(ppf1d[ir].real()));
             swq222.write((char*)&ppf1d[ir].imag(),sizeof(ppf1d[ir].imag()));
           } 

            cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Primary Estimation Done..."<<endl;
          }//6666

          else
           {//7777
 
             for(ir=0;ir<nrnp;ir++)
               mfp1d[ir]=(0.0,0.0);

             swq155.seekg(0,ios::beg);
             swq166.seekg(0,ios::beg);

             for(it1=ifmin;it1<it;it1++)
              {
               swq155.seekg(nrnp*4,ios::cur);
               swq166.seekg(nrnp*4,ios::cur);
              }

            for(ir=0;ir<nrnp;ir++)
             {
               swq155.write((char*)&mfp1d[ir].real(),sizeof(mfp1d[ir].real()));
               swq166.write((char*)&mfp1d[ir].imag(),sizeof(mfp1d[ir].imag()));
             } 

             for(ir=0;ir<lnrnp;ir++)
               wpf1d[ir]=(0.0,0.0);

             swq177.seekg(0,ios::beg);
             swq188.seekg(0,ios::beg);

             for(it1=ifmin;it1<it;it1++)
              {
               swq177.seekg(lnrnp*4,ios::cur);
               swq188.seekg(lnrnp*4,ios::cur);
              }

            for(ir=0;ir<lnrnp;ir++)
             {
               swq177.write((char*)&wpf1d[ir].real(),sizeof(wpf1d[ir].real()));
               swq188.write((char*)&wpf1d[ir].imag(),sizeof(wpf1d[ir].imag()));
             } 

             for(ir=0;ir<nrnp;ir++)
               mfp1dnew[ir]=(0.0,0.0);

             for(it1=ifmin;it1<it;it1++)
              {
               swq233.seekg(nrnp*4,ios::cur);
               swq244.seekg(nrnp*4,ios::cur);
              }

             for(ir=0;ir<nrnp;ir++)
              {
                swq233.write((char*)&mfp1dnew[ir].real(),sizeof(mfp1dnew[ir].real()));
                swq244.write((char*)&mfp1dnew[ir].imag(),sizeof(mfp1dnew[ir].imag()));
              } 

             for(ir=0;ir<nrnp;ir++)
               mcalfp1d[ir]=(0.0,0.0);

              swq199.seekg(0,ios::beg);
              swq200.seekg(0,ios::beg);

             for(it1=ifmin;it1<it;it1++)
              {
               swq199.seekg(nrnp*4,ios::cur);
               swq200.seekg(nrnp*4,ios::cur);
              }

            for(ir=0;ir<nrnp;ir++)
             {
               swq199.write((char*)&mcalfp1d[ir].real(),sizeof(mcalfp1d[ir].real()));
               swq200.write((char*)&mcalfp1d[ir].imag(),sizeof(mcalfp1d[ir].imag()));
             } 

             for(ir=0;ir<nrnp;ir++)
               ppf1d[ir]=(0.0,0.0);

            swq211.seekg(0,ios::beg);
            swq222.seekg(0,ios::beg);

            for(it1=ifmin;it1<it;it1++)
             {
               swq211.seekg(nrnp*4,ios::cur);
               swq222.seekg(nrnp*4,ios::cur);
             }

           for(ir=0;ir<nrnp;ir++)
            {
             swq211.write((char*)&ppf1d[ir].real(),sizeof(ppf1d[ir].real()));
             swq222.write((char*)&ppf1d[ir].imag(),sizeof(ppf1d[ir].imag()));
            } 

             cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;         
           }//7777


        }//8888

       else
        {
             for(ir=0;ir<nrnp;ir++)
               mfp1d[ir]=(0.0,0.0);

             swq155.seekg(0,ios::beg);
             swq166.seekg(0,ios::beg);

             for(it1=ifmin;it1<it;it1++)
              {
               swq155.seekg(nrnp*4,ios::cur);
               swq166.seekg(nrnp*4,ios::cur);
              }

            for(ir=0;ir<nrnp;ir++)
             {
               swq155.write((char*)&mfp1d[ir].real(),sizeof(mfp1d[ir].real()));
               swq166.write((char*)&mfp1d[ir].imag(),sizeof(mfp1d[ir].imag()));
             } 

             for(ir=0;ir<nrnp;ir++)
               mfp1dnew[ir]=(0.0,0.0);

             for(it1=ifmin;it1<it;it1++)
              {
               swq233.seekg(nrnp*4,ios::cur);
               swq244.seekg(nrnp*4,ios::cur);
              }

             for(ir=0;ir<nrnp;ir++)
              {
                swq233.write((char*)&mfp1dnew[ir].real(),sizeof(mfp1dnew[ir].real()));
                swq244.write((char*)&mfp1dnew[ir].imag(),sizeof(mfp1dnew[ir].imag()));
              } 

             for(ir=0;ir<lnrnp;ir++)
               wpf1d[ir]=(0.0,0.0);

             swq177.seekg(0,ios::beg);
             swq188.seekg(0,ios::beg);

             for(it1=ifmin;it1<it;it1++)
              {
               swq177.seekg(lnrnp*4,ios::cur);
               swq188.seekg(lnrnp*4,ios::cur);
              }

            for(ir=0;ir<lnrnp;ir++)
             {
               swq177.write((char*)&wpf1d[ir].real(),sizeof(wpf1d[ir].real()));
               swq188.write((char*)&wpf1d[ir].imag(),sizeof(wpf1d[ir].imag()));
             } 

             for(ir=0;ir<nrnp;ir++)
               mcalfp1d[ir]=(0.0,0.0);

              swq199.seekg(0,ios::beg);
              swq200.seekg(0,ios::beg);

             for(it1=ifmin;it1<it;it1++)
              {
               swq199.seekg(nrnp*4,ios::cur);
               swq200.seekg(nrnp*4,ios::cur);
              }

            for(ir=0;ir<nrnp;ir++)
             {
               swq199.write((char*)&mcalfp1d[ir].real(),sizeof(mcalfp1d[ir].real()));
               swq200.write((char*)&mcalfp1d[ir].imag(),sizeof(mcalfp1d[ir].imag()));
             } 

             for(ir=0;ir<nrnp;ir++)
               ppf1d[ir]=(0.0,0.0);

            swq211.seekg(0,ios::beg);
            swq222.seekg(0,ios::beg);

            for(it1=ifmin;it1<it;it1++)
             {
               swq211.seekg(nrnp*4,ios::cur);
               swq222.seekg(nrnp*4,ios::cur);
             }

           for(ir=0;ir<nrnp;ir++)
            {
             swq211.write((char*)&ppf1d[ir].real(),sizeof(ppf1d[ir].real()));
             swq222.write((char*)&ppf1d[ir].imag(),sizeof(ppf1d[ir].imag()));
            } 

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Weights Estimation Done..."<<endl;         
        }

    }//0000

    MPI_Barrier(MPI_COMM_WORLD);

    swq155.close();
    swq166.close();
    swq177.close();
    swq188.close();
    swq199.close();
    swq200.close();
    swq211.close();
    swq222.close();
    swq233.close();
    swq244.close();

      cout<<is+1<<" Shot Inverse Estimation Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











