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


int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn11[256],fn22[256],fn33[256],fn44[256],fn5[256],fn6[256],fn7[256];
   int nshot,nr,lt;
   float dt,fmax,fmin,st_beg,st_end;
   int ltt;
   int ifmax,ifmin,nsr,nsr1;
   int is,ir,it,ix,itmp;
  
   ifstream swq;
   swq.open("usrmg_qualitive_mcg_multis.par");
   if(!swq)
      {
         cout<<"cannot open usrmg_qualitive_mcg_multis.par"<<endl;
         abort();
      }

   swq>>fn1>>fn11>>fn2>>fn22>>fn3>>fn33>>fn4>>fn44>>fn5>>fn6>>fn7>>nshot>>nr>>lt>>dt>>fmax>>fmin>>st_beg>>st_end;
   swq.close();

   ltt=2*lt;
   nsr=nshot*nr;
   nsr1=nshot*(2*nr-1);

   float *omega;
   omega=alloc1float(ltt);
  
//calculate the omega and kx
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   ifmin=int(fmin*dt*ltt/1000);
   ifmax=int(fmax*dt*ltt/1000)+1;

   cout<<"Frequency slices to be calculated is==== "<<ifmax-ifmin+1<<endl;

   complex<float> *su;
   su=alloc1complex(nsr);
   complex<float> *su1;
   su1=alloc1complex(nr);
//   cout<<"=========="<<endl;

   complex<float> *sd;
   sd=alloc1complex(nsr);
    complex<float> *sd1;
   sd1=alloc1complex(nr);

   complex<float> *ru;
   ru=alloc1complex(nsr1);
   complex<float> *ru1;
   ru1=alloc1complex(nshot);

   complex<float> *rd;
   rd=alloc1complex(nsr1);
   complex<float> *rd1;
   rd1=alloc1complex(nshot);

   complex<float> *srmf;
   srmf=alloc1complex(nr);
 
   complex<float> *srmf1;
   srmf1=alloc1complex(ltt);
   
   float *srm;
   srm=alloc1float(ltt);

   float *srm1;
   srm1=alloc1float(lt);
 
   complex<float> mcg1;

   complex<float> *mcg2;
   mcg2=alloc1complex(ltt);

   float *mcg3;
   mcg3=alloc1float(ltt);
//   cout<<"=========="<<endl;

   fftwf_complex *in1,*out1;
   fftwf_plan p1;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_ESTIMATE);

//   cout<<"=========="<<endl;

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
      {
         cout<<"cannot open "<<fn1<<endl;
         abort();
      } 
   ifstream swq11;
   swq11.open(fn11,ios::binary);
   if(!swq11)
      {
         cout<<"cannot open "<<fn11<<endl;
         abort();
      }
   
   ifstream swq2;
   swq2.open(fn2,ios::binary);
   if(!swq2)
      {
         cout<<"cannot open "<<fn2<<endl;
         abort();
      } 
   ifstream swq22;
   swq22.open(fn22,ios::binary);
   if(!swq22)
      {
         cout<<"cannot open "<<fn22<<endl;
         abort();
      } 

   ifstream swq3;
   swq3.open(fn3,ios::binary);
   if(!swq3)
      {
         cout<<"cannot open "<<fn3<<endl;
         abort();
      } 
   ifstream swq33;
   swq33.open(fn33,ios::binary);
   if(!swq33)
      {
         cout<<"cannot open "<<fn33<<endl;
         abort();
      } 

   ifstream swq4;
   swq4.open(fn4,ios::binary);
   if(!swq4)
      {
         cout<<"cannot open "<<fn4<<endl;
         abort();
      } 
   ifstream swq44;
   swq44.open(fn44,ios::binary);
   if(!swq44)
      {
         cout<<"cannot open "<<fn44<<endl;
         abort();
      }
//   cout<<"=========="<<endl;
   
   ofstream swq6;
   swq6.open(fn6,ios::binary);
   if(!swq6)
      {
         cout<<"cannot open "<<fn6<<endl;
         abort();
      }

   ofstream swq7;
   swq7.open(fn7,ios::binary);
   if(!swq7)
      {
         cout<<"cannot open "<<fn7<<endl;
         abort();
      }

   ofstream swq5;
   swq5.open(fn5,ios::binary);
   if(!swq5)
      {
         cout<<"cannot open "<<fn5<<endl;
         abort();
      }

   cout<<"SRM Prediction of Positive Frequency Starts===="<<endl;

   for(it=0;it<lt+1;it++)
      {
          if((it+1)%10==0)
              cout<<it+1<<"th frequency slices prediction done!"<<endl;            

          if(it<ifmin)
            {
                for(is=0;is<nsr;is++)
                   srmf[is]=(0.0,0.0);
                for(is=0;is<nsr;is++)
                   {
                      swq6.write((char*)&srmf[is].real(),sizeof(srmf[is].real()));
                      swq7.write((char*)&srmf[is].imag(),sizeof(srmf[is].imag()));
                   }
            }
          else if(it>ifmax&&it<lt+1)
            {
                for(is=0;is<nsr;is++)
                   srmf[is]=(0.0,0.0);
                for(is=0;is<nsr;is++)
                   {
                      swq6.write((char*)&srmf[is].real(),sizeof(srmf[is].real()));
                      swq7.write((char*)&srmf[is].imag(),sizeof(srmf[is].imag()));
                   }
            }  
           else
            {
               swq1.seekg(it*4,ios::beg);
               swq1.read((char*)&(su[0].real()),sizeof(su[0].real()));
               swq11.seekg(it*4,ios::beg);
               swq11.read((char*)&(su[0].imag()),sizeof(su[0].imag()));
               swq2.seekg(it*4,ios::beg);
               swq2.read((char*)&(sd[0].real()),sizeof(sd[0].real()));
               swq22.seekg(it*4,ios::beg);
               swq22.read((char*)&(sd[0].imag()),sizeof(sd[0].imag()));
               swq3.seekg(it*4,ios::beg);
               swq3.read((char*)&(ru[0].real()),sizeof(ru[0].real()));
               swq33.seekg(it*4,ios::beg);
               swq33.read((char*)&(ru[0].imag()),sizeof(ru[0].imag()));
               swq4.seekg(it*4,ios::beg);
               swq4.read((char*)&(rd[0].real()),sizeof(rd[0].real()));
               swq44.seekg(it*4,ios::beg);
               swq44.read((char*)&(rd[0].imag()),sizeof(rd[0].imag()));

               for(is=1;is<nsr;is++)
                 {
                    swq1.seekg((ltt-1)*4,ios::cur);
                    swq1.read((char*)&(su[is].real()),sizeof(su[is].real()));
                    swq11.seekg((ltt-1)*4,ios::cur);
                    swq11.read((char*)&(su[is].imag()),sizeof(su[is].imag()));
                    swq2.seekg((ltt-1)*4,ios::cur);
                    swq2.read((char*)&(sd[is].real()),sizeof(sd[is].real()));
                    swq22.seekg((ltt-1)*4,ios::cur);
                    swq22.read((char*)&(sd[is].imag()),sizeof(sd[is].imag()));
                 }

              for(is=1;is<nsr1;is++)
                 {
                    swq3.seekg((ltt-1)*4,ios::cur);
                    swq3.read((char*)&(ru[is].real()),sizeof(ru[is].real()));
                    swq33.seekg((ltt-1)*4,ios::cur);
                    swq33.read((char*)&(ru[is].imag()),sizeof(ru[is].imag()));
                    swq4.seekg((ltt-1)*4,ios::cur);
                    swq4.read((char*)&(rd[is].real()),sizeof(rd[is].real()));
                    swq44.seekg((ltt-1)*4,ios::cur);
                    swq44.read((char*)&(rd[is].imag()),sizeof(rd[is].imag()));
                 }
 
                for(is=0;is<nshot;is++)
                  {
                     for(ir=0;ir<nr;ir++)
                       {
                          su1[ir].real()=su[is*nr+ir].real();
                          su1[ir].imag()=su[is*nr+ir].imag();
                          sd1[ir].real()=sd[is*nr+ir].real();
                          sd1[ir].imag()=sd[is*nr+ir].imag();
                       }

                     for(ir=0;ir<nr;ir++)
                       {
                          itmp=is*nr+ir;
                          srmf[itmp]=(0.0,0.0);

                          for(ix=0;ix<nshot;ix++)
                             {
                                ru1[ix].real()=ru[is*nshot+ir*nshot+ix].real();
                                ru1[ix].imag()=ru[is*nshot+ir*nshot+ix].imag();
                                rd1[ix].real()=rd[is*nshot+ir*nshot+ix].real();
                                rd1[ix].imag()=rd[is*nshot+ir*nshot+ix].imag();
                             }
                  
                          if(ir>=is)
                            {
                                 for(ix=is;ix<is+int((ir-is)*st_beg);ix++)
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[itmp]+=mcg1;
                                   }
                                 for(ix=is+int((ir-is)*st_beg);ix<is+int((ir-is)*st_end)+1;ix++)
                                   {
//                                     mcg1=su1[ix]*rd1[ix]+sd1[ix]*ru1[ix];
                                     mcg1=su1[ix]*rd1[ix];
                                     srmf[itmp]+=mcg1;
                                   }

                                 for(ix=is+int((ir-is)*st_end)+1;ix<ir+1;ix++)
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[itmp]+=mcg1;
                                   }

                               srmf[itmp].real()/=(int((ir-is)*st_end)-int((ir-is)*st_beg)+1);
                               srmf[itmp].imag()/=(int((ir-is)*st_end)-int((ir-is)*st_beg)+1);
//                                 srmf[itmp].imag()/=(ir-is+1);
//                                 srmf[itmp].real()/=(ir-is+1);

                            }
                           else
                            {
                                 for(ix=ir;ix<ir+int((is-ir)*st_beg);ix++) 
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[itmp]+=mcg1;                             
                                   }
                                 for(ix=ir+int((is-ir)*st_beg);ix<ir+int((is-ir)*st_end)+1;ix++)
                                   {
//                                     mcg1=su1[ix]*rd1[ix]+sd1[ix]*ru1[ix];
                                     mcg1=su1[ix]*rd1[ix];
                                     srmf[itmp]+=mcg1;
                                   }
                                 for(ix=ir+int((is-ir)*st_end)+1;ix<is+1;ix++)
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[itmp]+=mcg1;
                                   }

                               srmf[itmp].real()/=(int((is-ir)*st_end)-int((is-ir)*st_beg)+1);
                               srmf[itmp].imag()/=(int((is-ir)*st_end)-int((is-ir)*st_beg)+1);
//                                 srmf[itmp].real()/=(is-ir+1);
//                                 srmf[itmp].imag()/=(is-ir+1);

                            }

                       }

                }  

              for(is=0;is<nsr;is++)
                 {
                    swq6.write((char*)&(srmf[is].real()),sizeof(srmf[is].real()));
                    swq7.write((char*)&(srmf[is].imag()),sizeof(srmf[is].imag()));
                 }

            }//end of else

          }//end of it for positive frequency


           swq1.close();
           swq11.close();
           swq2.close();
           swq22.close();
           swq3.close();
           swq33.close();
           swq4.close();
           swq44.close();
           swq6.close();
           swq7.close();

          cout<<"SRM Prediction of Positive Frequency Done!"<<endl;
          cout<<"SRM Prediction of Negative Frequency Starts===="<<endl;

//negtive frequency starts; 
          ifstream swq66;
          swq66.open(fn6,ios::binary);
          if(!swq66)
             {
                cout<<"cannot open "<<fn6<<endl;
                abort();
             }

          ifstream swq77;
          swq77.open(fn7,ios::binary);
          if(!swq77)
            {
                cout<<"cannot open "<<fn7<<endl;
                abort();
            }
          for(is=0;is<nsr;is++)
            {
                for(it=0;it<ltt;it++)
                    srmf1[it]=(0.0,0.0);
            
                swq66.seekg(is*4,ios::beg);
                swq66.read((char*)&(srmf1[0].real()),sizeof(srmf1[0].real()));
                swq77.seekg(is*4,ios::beg);
                swq77.read((char*)&(srmf1[0].imag()),sizeof(srmf1[0].imag()));

                for(it=1;it<lt+1;it++)
                   {
                      swq66.seekg((nsr-1)*4,ios::cur);
                      swq66.read((char*)&(srmf1[it].real()),sizeof(srmf1[it].real()));
                      swq77.seekg((nsr-1)*4,ios::cur);
                      swq77.read((char*)&(srmf1[it].imag()),sizeof(srmf1[it].imag()));
                   }
                for(it=lt+1;it<ltt;it++)
                   {
                      srmf1[it].real()=srmf1[ltt-it].real();
                      srmf1[it].imag()=-srmf1[ltt-it].imag(); 
                   }
 
                if((in1==NULL)||(out1==NULL))
                   cout<<"memory insufficient"<<endl;
                else
                  {
                     for(it=0;it<ltt;it++)
                       {
                         in1[it][0]=srmf1[it].real();
                         in1[it][1]=srmf1[it].imag();
                       }
                  }

                fftwf_execute(p1);

                for(it=0;it<ltt;it++)
                   srm[it]=out1[it][0];

                for(it=0;it<lt;it++) 
                   swq5.write((char*)&(srm[it]),sizeof(srm[it]));

            }


          swq66.close();
          swq77.close();
          swq5.close();

          cout<<"ALL DONE!"<<endl; 
 
      return 0;
}
























































































































