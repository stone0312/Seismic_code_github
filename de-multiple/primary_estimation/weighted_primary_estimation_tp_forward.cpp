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
  char fn1[256],fn2[256],fn111[256],fn222[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fn12[256];
  int ns, nr,lt, nr2, coe1,l,nrl,np;
  float dt,fmin,fmax,lamda,lfbeg,lfend,hfbeg,hfend,dt1;

  int npoly=4;

  ifstream swq;
  swq.open("weighted_primary_estimation_tp_forward.par");
  swq>>fn111>>fn222>>fn3>>fn4>>fn5>>ns>>nr>>lt>>np>>dt>>fmin>>fmax>>lamda>>coe1>>l>>lfbeg>>lfend>>hfbeg>>hfend;
  swq.close();

  cout<<"Fna of Initial Primary Real is===="<<fn111<<endl; 
  cout<<"Fna of Initial Primary Imaginary is===="<<fn222<<endl; 
  cout<<"Fna of CRG Real is===="<<fn3<<endl; 
  cout<<"Fna of CRG Imaginary is===="<<fn4<<endl; 
  cout<<"Fna of Forward Prediction is===="<<fn5<<endl; 
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
  complex<float> *af1;
  af1=alloc1complex(nrnp);
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

  float **mt;
  mt=alloc2float(lt,nr);

  complex<float> **a;
  a=alloc2complex(lnrnp,nrnp);

  fftwf_complex *in1,*out1;
  fftwf_plan p1;
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);
 
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

  ofstream swq5;
  swq5.open(fn5,ios::binary);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;
 
          for(ir=0;ir<lnrnp;ir++)
             {
               af[ir].real()=1.0;
               af[ir].imag()=0.0;
             }

          for(ir=0;ir<nrnp;ir++)
             {
               af1[ir].real()=1.0;
               af1[ir].imag()=0.0;
             }

//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
       for(ir=0;ir<nrnp;ir++)
        for(it=0;it<ltt;it++)
           mfp[ir][it]=(0.0,0.0);

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

      swq111.seekg(0,ios::beg);
      swq222.seekg(0,ios::beg);

      for(ir=0;ir<is;ir++)
        {
           swq111.seekg(nrnp*ltt*4,ios::cur);
           swq222.seekg(nrnp*ltt*4,ios::cur);
        }

      for(ir=0;ir<nrnp;ir++)
        for(it=0;it<ltt;it++)
         {
           swq111.read((char*)&usfp1[ir][it].real(),sizeof(usfp1[ir][it].real()));  //primary 
           swq222.read((char*)&usfp1[ir][it].imag(),sizeof(usfp1[ir][it].imag())); 
         }

      for(ir=0;ir<nrnp;ir++)
        for(it=0;it<ltt;it++)
         {
            usfp1[ir][it].real()*=coe1;
            usfp1[ir][it].imag()*=coe1;
         } 

      for(it=ifmin;it<ifmax;it++)
        {//0000
          for(ip=0;ip<np;ip++)
           {
            for(ir=0;ir<nr;ir++)
               p0[ip*nr+ir]=usfp1[ir*np+ip][it];
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
//         complex_matrix_multiply_vector(mfp1d, pp0, af1, nrnp, nrnp);

          for(ir=0;ir<nrnp;ir++)
            mfp[ir][it]=mfp1d[ir];

          cout<<"    "<<is+1<<" Shot, "<<it+1<<" Frequency Slice Multiple Prediction Done..."<<endl;

    }//0000
    
      for(ir=0;ir<nrnp;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mfp[ir][it].real()=mfp[ir][ltt-it].real();
          mfp[ir][it].imag()=-mfp[ir][ltt-it].imag();
        }

    for(ir=0;ir<nrnp;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mfp[ir][it].real();
         in1[it][1]=mfp[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mtp[ir][it]=out1[it][0]/ltt;
     }

      for(ir=0;ir<nrnp;ir++)
       for(it=0;it<lt;it++)
         mtp[ir][it]*=-1.0;

//inverse tau-p
    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mt[ir][it]=0.0; 

    for(it=0;it<lt;it++)
      {
        for(ir=0;ir<nr;ir++)
           {
               for(ip=0;ip<np;ip++)
                  mt[ir][it]+=mtp[ir*np+ip][it];
           }
      }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mt[ir][it]/=np; 

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq5.write((char*)&mt[ir][it],sizeof(mt[ir][it]));

      cout<<is+1<<" Shot Forward Prediction Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











