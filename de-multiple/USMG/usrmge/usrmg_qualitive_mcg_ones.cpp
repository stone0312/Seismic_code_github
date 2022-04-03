#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

using namespace std;

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn11[256],fn22[256],fn33[256],fn44[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn100[256],fn1000[256];
   int nshot,nr,lt,ns1,nt;
   float dt,fmax,fmin,st_beg,st_end;
   int ltt;
   int ifmax,ifmin,nsr;
   int is,ir,it,ix,itmp,icount,icount1,abs1;
  
   ifstream swq;
   swq.open("usrmg_qualitive_mcg_ones.par");
   if(!swq)
      {
         cout<<"cannot open usrmg_qualitive_mcg_ones.par"<<endl;
         abort();
      }

   swq>>fn1>>fn11>>fn2>>fn22>>fn3>>fn33>>fn4>>fn44>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>fn100>>fn1000>>nshot>>nr>>lt>>dt>>fmax>>fmin>>ns1>>nt>>st_beg>>st_end;
   swq.close();

   ltt=2*lt;
   nsr=nshot*nr;
   
   icount=(1+(nr-ns1))*(nr-ns1)/2+((1+(ns1+1))*(ns1+1)/2-1);
   cout<<"No. of total traces for MCG of shot "<<ns1<<" is===="<<icount<<endl;

   if(nt<ns1)
     icount1=((ns1+1)+(ns1-(nt-1)+1))*nt/2;
   else if(nt==ns1)
     icount1=((ns1+1)+2)*ns1/2;
   else
     icount1=((ns1+1)+1)*(ns1+1)/2+(2+(nt-1+ns1-1))*(nt-ns1)/2;

   abs1=abs(nt-ns1)+1;
   cout<<abs1<<endl;


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

   complex<float> **su;
   su=alloc2complex(ltt,nr);
   complex<float> *su1;
   su1=alloc1complex(nr);

   complex<float> **sd;
   sd=alloc2complex(ltt,nr);
    complex<float> *sd1;
   sd1=alloc1complex(nr);

   complex<float> *ru;
   ru=alloc1complex(nsr);
   complex<float> *ru1;
   ru1=alloc1complex(nshot);

   complex<float> *rd;
   rd=alloc1complex(nsr);
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

   float **mcg_cer;
   mcg_cer=alloc2float(lt,abs1);
 
   complex<float> mcg1;

   complex<float> *mcg;
   mcg=alloc1complex(icount);

   complex<float> *mcg2;
   mcg2=alloc1complex(ltt);

   float *mcg3;
   mcg3=alloc1float(ltt);

   fftwf_complex *in1,*out1;
   fftwf_plan p1;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);


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

   ofstream swq8;
   swq8.open(fn8,ios::binary);
   if(!swq8)
      {
         cout<<"cannot open "<<fn8<<endl;
         abort();
      }

   ofstream swq9;
   swq9.open(fn9,ios::binary);
   if(!swq9)
      {
         cout<<"cannot open "<<fn9<<endl;
         abort();
      }

   ofstream swq10;
   swq10.open(fn10,ios::binary);
   if(!swq10)
      {
         cout<<"cannot open "<<fn10<<endl;
         abort();
      }


   cout<<"SRM Prediction of Positive Frequency Starts===="<<endl;
   cout<<"start to read the "<<ns1<<"th shot data...."<<endl;

   swq1.seekg(0,ios::beg);
   swq11.seekg(0,ios::beg);
   swq2.seekg(0,ios::beg);
   swq22.seekg(0,ios::beg);

   for(is=0;is<ns1;is++)
     {
         swq1.seekg(nr*ltt*4,ios::cur);
         swq11.seekg(nr*ltt*4,ios::cur);
         swq2.seekg(nr*ltt*4,ios::cur);
         swq22.seekg(nr*ltt*4,ios::cur); 
     }

   for(ir=0;ir<nr;ir++)
     for(it=0;it<ltt;it++)
      {
         swq1.read((char*)&su[ir][it].real(),sizeof(su[ir][it].real()));
         swq11.read((char*)&su[ir][it].imag(),sizeof(su[ir][it].imag()));
         swq2.read((char*)&sd[ir][it].real(),sizeof(sd[ir][it].real()));
         swq22.read((char*)&sd[ir][it].imag(),sizeof(sd[ir][it].imag()));
      } 

/*
   ofstream swqaa;
   swqaa.open("/data2/swq/usrmge/depth_200m/test1.dat",ios::binary);
   for(ir=0;ir<nr;ir++)
     for(it=0;it<ltt;it++)
        swqaa.write((char*)&su[ir][it].real(),sizeof(su[ir][it].real())); 
   swqaa.close();

   ofstream swqbb;
   swqbb.open("/data2/swq/usrmge/depth_200m/test2.dat",ios::binary);
   for(ir=0;ir<nr;ir++)
     for(it=0;it<ltt;it++)
        swqbb.write((char*)&su[ir][it].imag(),sizeof(su[ir][it].imag()));
   swqbb.close();

   return 0;
*/

   cout<<ns1<<"th shot data reading done!"<<endl;

   for(it=0;it<lt+1;it++)
      {
          if((it+1)%10==0)
              cout<<it+1<<"th frequency slices prediction done!"<<endl;            

          if(it<ifmin)
            {
                for(is=0;is<nr;is++)
                   srmf[is]=(0.0,0.0);
                for(is=0;is<nr;is++)
                   {
                      swq6.write((char*)&srmf[is].real(),sizeof(srmf[is].real()));
                      swq7.write((char*)&srmf[is].imag(),sizeof(srmf[is].imag()));
                   }
                for(is=0;is<icount;is++)
                   mcg[is]=(0.0,0.0);
                for(is=0;is<icount;is++)
                   {
                      swq8.write((char*)&mcg[is].real(),sizeof(mcg[is].real()));
                      swq9.write((char*)&mcg[is].imag(),sizeof(mcg[is].imag()));
                   }
            }
          else if(it>ifmax&&it<lt+1)
            {
                for(is=0;is<nr;is++)
                   srmf[is]=(0.0,0.0);
                for(is=0;is<nr;is++)
                   {
                      swq6.write((char*)&srmf[is].real(),sizeof(srmf[is].real()));
                      swq7.write((char*)&srmf[is].imag(),sizeof(srmf[is].imag()));
                   }
                for(is=0;is<icount;is++)
                   mcg[is]=(0.0,0.0);
                for(is=0;is<icount;is++)
                   {
                      swq8.write((char*)&mcg[is].real(),sizeof(mcg[is].real()));
                      swq9.write((char*)&mcg[is].imag(),sizeof(mcg[is].imag()));
                   }
            }  
           else
            {
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
                    swq3.seekg((ltt-1)*4,ios::cur);
               	    swq3.read((char*)&(ru[is].real()),sizeof(ru[is].real()));
                    swq33.seekg((ltt-1)*4,ios::cur);
                    swq33.read((char*)&(ru[is].imag()),sizeof(ru[is].imag()));
                    swq4.seekg((ltt-1)*4,ios::cur);
                    swq4.read((char*)&(rd[is].real()),sizeof(rd[is].real()));
                    swq44.seekg((ltt-1)*4,ios::cur);
                    swq44.read((char*)&(rd[is].imag()),sizeof(rd[is].imag()));
                 }
              
                     for(ir=0;ir<nr;ir++)
                       {
                          su1[ir].real()=su[ir][it].real();
                          su1[ir].imag()=su[ir][it].imag();
                          sd1[ir].real()=sd[ir][it].real();
                          sd1[ir].imag()=sd[ir][it].imag();
                       } 

                     for(ir=0;ir<nr;ir++)
                       {
                          srmf[ir]=(0.0,0.0);

                          for(ix=0;ix<nshot;ix++)
                             {
                                ru1[ix].real()=ru[ir*nshot+ix].real();
                                ru1[ix].imag()=ru[ir*nshot+ix].imag();
                                rd1[ix].real()=rd[ir*nshot+ix].real();
                                rd1[ix].imag()=rd[ir*nshot+ix].imag();
                             }
                  
                          if(ir>=ns1)
                            {
                                 for(ix=ns1;ix<ns1+int((ir-ns1)*st_beg+0.5);ix++)
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[ir]+=mcg1;
                                     swq8.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                                     swq9.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));
                                   }
                                 for(ix=ns1+int((ir-ns1)*st_beg+0.5);ix<ns1+int((ir-ns1)*st_end+0.5)+1;ix++)
                                   {
//source and receiver-side wlrm 
                                     mcg1=su1[ix]*rd1[ix]+sd1[ix]*ru1[ix];
//source-side wlrm
//                                     mcg1=sd1[ix]*ru1[ix];
//receiver-side wlrm
//                                     mcg1=su1[ix]*rd1[ix];

                                     srmf[ir]+=mcg1;

                                     swq8.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                                     swq9.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));
                                   }
                                 for(ix=ns1+int((ir-ns1)*st_end+0.5)+1;ix<ir+1;ix++)
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[ir]+=mcg1;
                                     swq8.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                                     swq9.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));
                                   }
                                 srmf[ir].real()/=(int((ir-ns1)*st_end)-int((ir-ns1)*st_beg)+1);
                                 srmf[ir].imag()/=(int((ir-ns1)*st_end)-int((ir-ns1)*st_beg)+1);

                            }
                           else
                            {
                                 for(ix=ir;ix<ir+int((ns1-ir)*st_beg+0.5);ix++) 
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[ir]+=mcg1;                             
                                     swq8.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                                     swq9.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));
                                   }
                                 for(ix=ir+int((ns1-ir)*st_beg+0.5);ix<ir+int((ns1-ir)*st_end+0.5)+1;ix++)
                                   {
//source and receiver-side wlrm 
                                     mcg1=su1[ix]*rd1[ix]+sd1[ix]*ru1[ix];
//source-side wlrm
//                                     mcg1=sd1[ix]*ru1[ix];
//receiver-side wlrm
//                                     mcg1=su1[ix]*rd1[ix];

                                     srmf[ir]+=mcg1;
                                     
                                     swq8.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                                     swq9.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));
                                   }
                                 for(ix=ir+int((ns1-ir)*st_end+0.5)+1;ix<ns1+1;ix++)
                                   {
                                     mcg1=(0.0,0.0);
                                     srmf[ir]+=mcg1;
                                     swq8.write((char*)&(mcg1.real()),sizeof(mcg1.real()));
                                     swq9.write((char*)&(mcg1.imag()),sizeof(mcg1.imag()));
                                   }
                                 srmf[ir].real()/=(int((ns1-ir)*st_end)-int((ns1-ir)*st_beg)+1);
                                 srmf[ir].imag()/=(int((ns1-ir)*st_end)-int((ns1-ir)*st_beg)+1);

                            }

                       }

              for(is=0;is<nr;is++)
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
           swq8.close();
           swq9.close();

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
      
          ifstream swq88;
          swq88.open(fn8,ios::binary);
          if(!swq88)
             {
                cout<<"cannot open "<<fn8<<endl;
                abort();
             }

          ifstream swq99;
          swq99.open(fn9,ios::binary);
          if(!swq99)
            {
                cout<<"cannot open "<<fn9<<endl;
                abort();
            }
 
          for(is=0;is<nr;is++)
            {
                for(it=0;it<ltt;it++)
                    srmf1[it]=(0.0,0.0);
            
                swq66.seekg(is*4,ios::beg);
                swq66.read((char*)&(srmf1[0].real()),sizeof(srmf1[0].real()));
                swq77.seekg(is*4,ios::beg);
                swq77.read((char*)&(srmf1[0].imag()),sizeof(srmf1[0].imag()));

                for(it=1;it<lt+1;it++)
                   {
                      swq66.seekg((nr-1)*4,ios::cur);
                      swq66.read((char*)&(srmf1[it].real()),sizeof(srmf1[it].real()));
                      swq77.seekg((nr-1)*4,ios::cur);
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

         for(is=0;is<icount;is++)
            {
                for(it=0;it<ltt;it++)
                    mcg2[it]=(0.0,0.0);

                swq88.seekg(is*4,ios::beg);
                swq88.read((char*)&(mcg2[0].real()),sizeof(mcg2[0].real()));
                swq99.seekg(is*4,ios::beg);
                swq99.read((char*)&(mcg2[0].imag()),sizeof(mcg2[0].imag()));

                for(it=1;it<lt+1;it++)
                   {
                      swq88.seekg((icount-1)*4,ios::cur);
                      swq88.read((char*)&(mcg2[it].real()),sizeof(mcg2[it].real()));
                      swq99.seekg((icount-1)*4,ios::cur);
                      swq99.read((char*)&(mcg2[it].imag()),sizeof(mcg2[it].imag()));
                   }
                for(it=lt+1;it<ltt;it++)
                   {
                      mcg2[it].real()=mcg2[ltt-it].real();

                      if(mcg2[ltt-it].imag()==0.0) 
                        mcg2[it].imag()=mcg2[ltt-it].imag();
                      else
                        mcg2[it].imag()=-mcg2[ltt-it].imag();
                   }

                if((in1==NULL)||(out1==NULL))
                   cout<<"memory insufficient"<<endl;
                else
                  {
                     for(it=0;it<ltt;it++)
                       {
                         in1[it][0]=mcg2[it].real();
                         in1[it][1]=mcg2[it].imag();
                       }
                  }

                fftwf_execute(p1);

                for(it=0;it<ltt;it++)
                   mcg3[it]=out1[it][0];

                for(it=0;it<lt;it++)
                   swq10.write((char*)&(mcg3[it]),sizeof(mcg3[it]));

            }


          swq66.close();
          swq77.close();
          swq88.close();
          swq99.close();
          swq10.close();
          swq5.close();

          ifstream swq200;
          swq200.open(fn10,ios::binary);
          if(!swq200)
            {
              cout<<"cannot open "<<fn10<<endl;
              abort();
            }

          ifstream swq2000;
          swq2000.open(fn5,ios::binary);
          if(!swq2000)
            {
              cout<<"cannot open "<<fn5<<endl;
              abort();
           }
         
           ofstream swq100;
           swq100.open(fn100,ios::binary);
           if(!swq100)
             {
                cout<<"cannot open "<<fn100<<endl;
                abort();
             }

           ofstream swq1000;
           swq1000.open(fn1000,ios::binary);
           if(!swq1000)
             {
                cout<<"cannot open "<<fn1000<<endl;
                abort();
             }

           swq200.seekg(0,ios::beg);
           swq2000.seekg(0,ios::beg);
           for(ir=0;ir<nt;ir++)
               swq2000.seekg(lt*4,ios::cur);
           for(ir=0;ir<icount1;ir++)
               swq200.seekg(lt*4,ios::cur);

           for(it=0;it<lt;it++) 
              swq2000.read((char*)&srm1[it],sizeof(srm1[it]));
           for(it=0;it<lt;it++)
              swq1000.write((char*)&srm1[it],sizeof(srm1[it]));

            
           for(ir=0;ir<abs1;ir++)
              for(it=0;it<lt;it++)
                 swq200.read((char*)&mcg_cer[ir][it],sizeof(mcg_cer[ir][it]));

           for(ir=0;ir<abs1;ir++)
              for(it=0;it<lt;it++)
                 swq100.write((char*)&mcg_cer[ir][it],sizeof(mcg_cer[ir][it]));


          swq100.close();
          swq200.close();
          swq1000.close();
          swq2000.close();


          cout<<"ALL DONE!"<<endl; 
 
      return 0;
}
























































































































