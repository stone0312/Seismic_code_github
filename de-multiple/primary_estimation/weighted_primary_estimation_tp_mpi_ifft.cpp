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
  swq.open("weighted_primary_estimation_tp_mpi_ifft.par");
  swq>>fn5>>fn6>>fn7>>fn8>>fn11>>fn12>>fn15>>fn16>>fn17>>fn18>>fn19>>fn20>>fn21>>fn22>>fn23>>fn24>>ns>>nr>>lt>>np>>dt>>fmin>>fmax>>lamda>>coe1>>l>>lfbeg>>lfend>>hfbeg>>hfend;
  swq.close();

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
 
  complex<float> **mfp;
  mfp=alloc2complex(ltt,nrnp);
  complex<float> *af;
  af=alloc1complex(lnrnp);

  complex<float> **p;
  p=alloc2complex(nrnp,nrnp);

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
  float **wt2d;
  wt2d=alloc2float(lt,nnrnp);

  float **mt;
  mt=alloc2float(lt,nr);

  complex<float> **mfptrue;
  mfptrue=alloc2complex(ltt,nrnp);

  complex<float> **mcalfp2d;
  mcalfp2d=alloc2complex(ltt,nrnp);
  float **mcalpt2d;
  mcalpt2d=alloc2float(lt,nrnp);
  float **mcalt2d;
  mcalt2d=alloc2float(lt,nr);
  complex<float> **diffp2d;
  diffp2d=alloc2complex(ltt,nrnp);
  float **diftp2d;
  diftp2d=alloc2float(lt,nrnp);
  float **dift2d;
  dift2d=alloc2float(lt,nr);

  complex<float> **mfp2dnew;
  mfp2dnew=alloc2complex(ltt,nrnp);
  float **mtp2dnew;
  mtp2dnew=alloc2float(lt,nrnp);
  float **mt2dnew;
  mt2dnew=alloc2float(lt,nr);

  complex<float> **ppf2d;
  ppf2d=alloc2complex(ltt,nrnp);

  float **ppt2d;
  ppt2d=alloc2float(lt,nrnp);
  float **pt2d;
  pt2d=alloc2float(lt,nr);

  fftwf_complex *in1,*out1;
  fftwf_plan p1;
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
  p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);

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
 
//  for(is=0;is<ns;is++)
  for(is=49;is<50;is++)
    {
       for(ir=0;ir<nnrnp;ir++)
        for(it=0;it<ltt;it++)
          wpf2dtmp[ir][it]=(0.0,0.0);

       for(ir=0;ir<nrnp;ir++)
        for(it=0;it<ltt;it++)
         {
           ppf2d[ir][it]=(0.0,0.0);    
           mcalfp2d[ir][it]=(0.0,0.0);    
           mfp2dnew[ir][it]=(0.0,0.0);    
           mfp[ir][it]=(0.0,0.0);
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

  ifstream swq15; //read calculated mfp1d
  swq15.open(fn15,ios::binary);
  if(!swq15)
       cout<<"cannot open "<<fn15<<endl;
 
  ifstream swq16;  //read calculated wpf1d
  swq16.open(fn16,ios::binary);
  if(!swq16)
       cout<<"cannot open "<<fn16<<endl;
 
  ifstream swq17;   //read calculated mfp1dnew
  swq17.open(fn17,ios::binary);
  if(!swq17)
       cout<<"cannot open "<<fn17<<endl;
 
  ifstream swq18;   //read calculated ppf1d
  swq18.open(fn18,ios::binary);
  if(!swq18)
       cout<<"cannot open "<<fn18<<endl;
 
  ifstream swq19; //read calculated mfp1d
  swq19.open(fn19,ios::binary);
  if(!swq19)
       cout<<"cannot open "<<fn19<<endl;
 
  ifstream swq20;  //read calculated wpf1d
  swq20.open(fn20,ios::binary);
  if(!swq20)
       cout<<"cannot open "<<fn20<<endl;
 
  ifstream swq21;   //read calculated mfp1dnew
  swq21.open(fn21,ios::binary);
  if(!swq21)
       cout<<"cannot open "<<fn21<<endl;
 
  ifstream swq22;   //read calculated ppf1d
  swq22.open(fn22,ios::binary);
  if(!swq22)
       cout<<"cannot open "<<fn22<<endl;

  ifstream swq23;   //read calculated mfp1dnew
  swq23.open(fn23,ios::binary);
  if(!swq23)
       cout<<"cannot open "<<fn23<<endl;
 
  ifstream swq24;   //read calculated ppf1d
  swq24.open(fn24,ios::binary);
  if(!swq24)
       cout<<"cannot open "<<fn24<<endl;

  for(ir=0;ir<nrnp;ir++)
    for(it=0;it<ltt;it++)
      {
         mfp[ir][it]=(0.0,0.0);
         mcalfp2d[ir][it]=(0.0,0.0);
         mfp2dnew[ir][it]=(0.0,0.0);
         ppf2d[ir][it]=(0.0,0.0);
      }

  for(ir=0;ir<nnrnp;ir++)
    for(it=0;it<ltt;it++)
       wpf2d[ir][it]=(0.0,0.0);

  for(it=ifmin;it<ifmax;it++)
    {
       for(ir=0;ir<nrnp;ir++)
          {
             swq15.read((char*)&mfp[ir][it].real(),sizeof(mfp[ir][it].real()));
             swq16.read((char*)&mfp[ir][it].imag(),sizeof(mfp[ir][it].imag()));
             swq19.read((char*)&mcalfp2d[ir][it].real(),sizeof(mcalfp2d[ir][it].real()));
             swq20.read((char*)&mcalfp2d[ir][it].imag(),sizeof(mcalfp2d[ir][it].imag()));
             swq21.read((char*)&ppf2d[ir][it].real(),sizeof(ppf2d[ir][it].real()));
             swq22.read((char*)&ppf2d[ir][it].imag(),sizeof(ppf2d[ir][it].imag()));
             swq23.read((char*)&mfp2dnew[ir][it].real(),sizeof(mfp2dnew[ir][it].real()));
             swq24.read((char*)&mfp2dnew[ir][it].imag(),sizeof(mfp2dnew[ir][it].imag()));
          }
       for(ir=0;ir<lnrnp;ir++)
          {
             swq17.read((char*)&wpf1d[ir].real(),sizeof(wpf1d[ir].real()));
             swq18.read((char*)&wpf1d[ir].imag(),sizeof(wpf1d[ir].imag()));
          }

       for(ir=0;ir<nr;ir++)
          {
            irtmp1=irbeg[ir];
            for(ir1=0;ir1<np;ir1++)
              {
                 for(ir2=ir*nr*np+ir1*nr+irtmp1;ir2<ir*nr*np+ir1*nr+irtmp1+l;ir2++)
                   wpf2d[ir2][it]=wpf1d[ir*l*np+ir1*l+ir2-ir*nr*np-ir1*nr-irtmp1];
              }
          }
    }

      for(ir=0;ir<nrnp;ir++)
       for(it=0;it<lt;it++)
         {
            ppf2d[ir][it].real()*=filter[it];
            ppf2d[ir][it].imag()*=filter[it];
         }
    
      for(ir=0;ir<nrnp;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          mfp[ir][it].real()=mfp[ir][ltt-it].real();
          mfp[ir][it].imag()=-mfp[ir][ltt-it].imag();
          mcalfp2d[ir][it].real()=mcalfp2d[ir][ltt-it].real();
          mcalfp2d[ir][it].imag()=-mcalfp2d[ir][ltt-it].imag();
          mfp2dnew[ir][it].real()=mfp2dnew[ir][ltt-it].real();
          mfp2dnew[ir][it].imag()=-mfp2dnew[ir][ltt-it].imag();
          ppf2d[ir][it].real()=ppf2d[ir][ltt-it].real();
          ppf2d[ir][it].imag()=-ppf2d[ir][ltt-it].imag();
        }

      for(ir=0;ir<nnrnp;ir++)
       for(it=ltt/2+1;it<ltt;it++)
        {
          wpf2d[ir][it].real()=wpf2d[ir][ltt-it].real();
          wpf2d[ir][it].imag()=-wpf2d[ir][ltt-it].imag();
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
               for(ip=ir*np;ip<(ir+1)*np;ip++)
                  mt[ir][it]+=mtp[ip][it];
           }
      }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mt[ir][it]/=np; 

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq5.write((char*)&mt[ir][it],sizeof(mt[ir][it]));

      cout<<is+1<<" Shot Forward Prediction Done ..."<<endl;

    for(ir=0;ir<nnrnp;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=wpf2d[ir][it].real();
         in1[it][1]=wpf2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       wt2d[ir][it]=out1[it][0]/ltt;
     }

   for(ir=0;ir<nnrnp;ir++)
       for(it=0;it<lt;it++)
        swq6.write((char*)&wt2d[ir][it],sizeof(wt2d[ir][it]));

    for(ir=0;ir<nrnp;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mcalfp2d[ir][it].real()*coe1;
         in1[it][1]=mcalfp2d[ir][it].imag()*coe1;
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mtp[ir][it]=out1[it][0]/ltt;
     }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mcalt2d[ir][it]=0.0; 

    for(it=0;it<lt;it++)
      {
        for(ir=0;ir<nr;ir++)
           {
               for(ip=ir*np;ip<(ir+1)*np;ip++)
                  mcalt2d[ir][it]+=mtp[ip][it];
           }
      }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mcalt2d[ir][it]/=np; 

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq7.write((char*)&mcalt2d[ir][it],sizeof(mcalt2d[ir][it]));

    for(ir=0;ir<nrnp;ir++)
      for(it=0;it<ltt;it++)
        {
           mfptrue[ir][it].real()/=coe1;
           mfptrue[ir][it].imag()/=coe1;
        }

    for(ir=0;ir<nrnp;ir++)
     for(it=0;it<ltt;it++)
       diffp2d[ir][it]=mcalfp2d[ir][it]-mfptrue[ir][it];

    for(ir=0;ir<nrnp;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=diffp2d[ir][it].real();
         in1[it][1]=diffp2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       diftp2d[ir][it]=out1[it][0]/ltt;
     }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         dift2d[ir][it]=0.0; 

    for(it=0;it<lt;it++)
      {
        for(ir=0;ir<nr;ir++)
           {
               for(ip=ir*np;ip<(ir+1)*np;ip++)
                  dift2d[ir][it]+=diftp2d[ip][it];
           }
      }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         dift2d[ir][it]/=np; 

      for(ir=0;ir<nr;ir++)
       for(it=0;it<lt;it++)
        swq8.write((char*)&dift2d[ir][it],sizeof(dift2d[ir][it]));


    for(ir=0;ir<nrnp;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=mfp2dnew[ir][it].real()*coe1;
         in1[it][1]=mfp2dnew[ir][it].imag()*coe1;
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       mtp2dnew[ir][it]=out1[it][0]/ltt;
     }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mt2dnew[ir][it]=0.0; 

    for(it=0;it<lt;it++)
      {
        for(ir=0;ir<nr;ir++)
           {
               for(ip=ir*np;ip<(ir+1)*np;ip++)
                  mt2dnew[ir][it]+=mtp2dnew[ip][it];
           }
      }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         mt2dnew[ir][it]/=np; 

    for(ir=0;ir<nr;ir++)
     for(it=0;it<lt;it++)
       swq11.write((char*)&mt2dnew[ir][it],sizeof(mt2dnew[ir][it]));

    for(ir=0;ir<nrnp;ir++)
     {
      for(it=0;it<ltt;it++)
       {
         in1[it][0]=ppf2d[ir][it].real();
         in1[it][1]=ppf2d[ir][it].imag();
       }

      fftwf_execute(p1);

      for(it=0;it<lt;it++)
       ppt2d[ir][it]=out1[it][0]/ltt;
     }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         pt2d[ir][it]=0.0; 

    for(it=0;it<lt;it++)
      {
        for(ir=0;ir<nr;ir++)
           {
               for(ip=ir*np;ip<(ir+1)*np;ip++)
                  pt2d[ir][it]+=ppt2d[ip][it];
           }
      }

    for(ir=0;ir<nr;ir++)
      for(it=0;it<lt;it++)
         pt2d[ir][it]/=np; 

    for(ir=0;ir<nr;ir++)
     for(it=0;it<lt;it++)
       swq12.write((char*)&pt2d[ir][it],sizeof(pt2d[ir][it]));

      cout<<is+1<<" Shot Inverse Estimation Done ..."<<endl;

    }
  
  cout<<"All Done!"<<endl;

  return 0;

}











