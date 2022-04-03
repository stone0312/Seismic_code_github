#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#define pai 3.14159265

using namespace std;

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256];
   int ns,nr,lt,itmax;
   float dt,fmax,fmin,zs,zr,dsx,drx,st_beg,st_end,dep,err,v,dtau_max,dtau_min,amp1,amp2;
   int ltt,nsr,it_count;
   int ifmax,ifmin;
   int is,ir,it,ix,itmp;

   ifstream swq;
   swq.open("usrmge.par");
   if(!swq)
      {
         cout<<"cannot open usrmge.par"<<endl;
         abort();
      }

   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>ns>>nr>>zs>>zr>>dsx>>drx>>lt>>dt>>fmin>>fmax>>st_beg>>st_end>>dep>>v>>itmax>>err;
   swq.close();
 
   nsr=ns*nr;
   ltt=2*lt;

   dtau_max=1/(4*fmax);
   dtau_min=-dtau_max;


   cout<<"fna of real parts of csg is===="<<fn1<<endl;
   cout<<"fna of imaginary parts of csg is===="<<fn2<<endl;
   cout<<"fna of real parts of crg is===="<<fn3<<endl;
   cout<<"fna of imaginary parts of crg is===="<<fn4<<endl;
   cout<<"fna of real parts of srmg for positive frequency is===="<<fn5<<endl;
   cout<<"fna of imaginary parts of srmg for positive frequency is===="<<fn6<<endl;
   cout<<"fna of real parts of desrmg for positive frequency is===="<<fn7<<endl;
   cout<<"fna of imaginary parts of desrmg for positive frequency is===="<<fn8<<endl;
   cout<<"fna of final results of srmg is===="<<fn9<<endl;
   cout<<"fna of final results of desrmg is===="<<fn10<<endl;  
   cout<<"nshot is===="<<ns<<endl;
   cout<<"ntrace is===="<<nr<<endl;
   cout<<"depth of shots is===="<<zs<<"m"<<endl;
   cout<<"depth of receivers is===="<<zr<<"m"<<endl;
   cout<<"interval of shots is===="<<dsx<<"m"<<endl;
   cout<<"interval of receivers is===="<<drx<<"m"<<endl;
   cout<<"No. of time samples is===="<<lt<<endl;
   cout<<"interval of time sampling is===="<<dt<<"ms"<<endl;
   cout<<"minimum and maximum frequency to be calculated is===="<<fmin<<"Hz , "<<fmax<<"Hz"<<endl;
   cout<<"starting and ending point for starting the MCG is===="<<st_beg<<" , "<<st_end<<endl;
   cout<<"depth of water is===="<<dep<<"m"<<endl;
   cout<<"velocity of water is===="<<v<<"m/s"<<endl;
   cout<<"maximum times allowing for iteration is===="<<itmax<<endl;
   cout<<"error for iteration is===="<<err<<endl;
   cout<<"maximum and minimum perturbation is==="<<dtau_max<<"s , "<<dtau_min<<"s"<<endl;

   complex<float> lamda;
   lamda.real()=0.001;
   lamda.imag()=0.0;

   float *omega;
   omega=alloc1float(ltt);

   float **tau;
   tau=alloc2float(nr,ns);
   complex<float> gh;

   float **tau1;
   tau1=allocfloat(nr);

   complex<float> *gr;
   gr=alloc1complex(nr);
   complex<float> *gs;
   gs=alloc1complex(ns);

   complex<float> *gr1;
   gr1=alloc1complex(nr);
   complex<float> *gs1;
   gs1=alloc1complex(ns);

   float **dtau;
   dtau=alloc2float(nr,ns);

   complex<float> *us;
   us=alloc1complex(nsr); 
   complex<float> *ur;
   ur=alloc1complex(nsr);

   complex<float> *us1;
   us1=alloc1complex(nr);
   complex<float> *ur1;
   ur1=alloc1complex(ns);

   complex<float> **pm;
   pm=alloc2complex(ns,ns);
   
   complex<float> *srmg;
   srmg=alloc1complex(nsr);
 
   float l2norm;
   l2norm=0.0;

   complex<float> *desrmg;
   desrmg=alloc1complex(nsr);

   complex<float> *srmgf;
   srmgf=alloc1complex(ltt);
   complex<float> *desrmgf;
   desrmgf=alloc1complex(ltt);
 
   float *srmgt;
   srmgt=alloc1float(ltt);
   float *desrmgt;
   desrmgt=alloc1float(ltt);

   fftwf_complex *in1,*out1;
   fftwf_plan p1;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);

   fftwf_complex *in2,*out2;
   fftwf_plan p2;
   in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p2=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_BACKWARD,FFTW_MEASURE);

   float mem=0.0;
   mem+=(ltt*4+ns*nr*4+ns*nr*4+nr*4*2+ns*4*2+nsr*4*2+nsr*4*2+nr*4*2+ns*4*2+ns*ns*4*2+nsr*4*2+nsr*4*2+ltt*4*2+ltt*4*2+ltt*4+ltt*4+ltt*4*2+ltt*4*2);
   mem=mem/4/1024/1024;
   cout<<"Memory to be allocated is===="<<mem<<"MB"<<endl; 

//calculate the omega
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   ifmin=int(fmin*dt*ltt/1000);
   ifmax=int(fmax*dt*ltt/1000)+1;

   cout<<"Frequency slices to be calculated is==== "<<ifmax-ifmin+1<<endl;

   for(is=0;is<ns;is++)
     for(ir=0;ir<nr;ir++)
       tau[is][ir]=(sqrt(4*dep*dep+abs(is-ir)*dsx*abs(is-ir)*dsx)-sqrt(4*(dep-zs)*(dep-zs)+abs(is-ir)*dsx*abs(is-ir)*dsx))/v;
/*
   for(is=0;is<ns;is++)
     {
        for(ir=0;ir<nr;ir++)
           cout<<tau[is][ir]<<"  ";
        cout<<endl;
     }
   return 0;
*/

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
      {
         cout<<"cannot open "<<fn1<<endl;
         abort();
      }

   ifstream swq2;
   swq2.open(fn2,ios::binary);
   if(!swq2)
      {
         cout<<"cannot open "<<fn2<<endl;
         abort();
      }

   ifstream swq3;
   swq3.open(fn3,ios::binary);
   if(!swq3)
      {
         cout<<"cannot open "<<fn3<<endl;
         abort();
      }

   ifstream swq4;
   swq4.open(fn4,ios::binary);
   if(!swq4)
      {
         cout<<"cannot open "<<fn4<<endl;
         abort();
      }

   ofstream swq5;
   swq5.open(fn5,ios::binary);
   if(!swq5)
      {
         cout<<"cannot open "<<fn5<<endl;
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

   ofstream swq8;
   swq8.open(fn8,ios::binary);
   if(!swq8)
      {
         cout<<"cannot open "<<fn8<<endl;
         abort();
      }
 
   cout<<"SRM Prediction of Positive Frequency Starts===="<<endl;

   for(it=0;it<lt+1;it++)
//7777for 
      {

          if((it+1)%10==0)
              cout<<it+1<<"th frequency slices prediction done!"<<endl;

          if(it<ifmin)
            {
                for(is=0;is<nsr;is++)
                  {
                    srmg[is]=(0.0,0.0);
                    desrmg[is]=(0.0,0.0);
                  }
                 for(is=0;is<nsr;is++)
                   {
                      swq5.write((char*)&srmg[is].real(),sizeof(srmg[is].real()));
                      swq6.write((char*)&srmg[is].imag(),sizeof(srmg[is].imag()));
                      swq7.write((char*)&desrmg[is].real(),sizeof(desrmg[is].real()));
                      swq8.write((char*)&desrmg[is].imag(),sizeof(desrmg[is].imag()));
                   }
            }
          else if(it>ifmax&&it<lt+1)
            {
               for(is=0;is<nsr;is++)
                  {
                    srmg[is]=(0.0,0.0);
                    desrmg[is]=(0.0,0.0);
                  }
                 for(is=0;is<nsr;is++)
                   {
                      swq5.write((char*)&srmg[is].real(),sizeof(srmg[is].real()));
                      swq6.write((char*)&srmg[is].imag(),sizeof(srmg[is].imag()));
                      swq7.write((char*)&desrmg[is].real(),sizeof(desrmg[is].real()));
                      swq8.write((char*)&desrmg[is].imag(),sizeof(desrmg[is].imag()));
                   }

            }
          else
//9999else
            {
               swq1.seekg(it*4,ios::beg);
               swq1.read((char*)&(us[0].real()),sizeof(us[0].real()));
               swq2.seekg(it*4,ios::beg);
               swq2.read((char*)&(us[0].imag()),sizeof(us[0].imag()));
               swq3.seekg(it*4,ios::beg);
               swq3.read((char*)&(ur[0].real()),sizeof(ur[0].real()));
               swq4.seekg(it*4,ios::beg);
               swq4.read((char*)&(ur[0].imag()),sizeof(ur[0].imag()));

               for(is=1;is<nsr;is++)
                 {
                    swq1.seekg((ltt-1)*4,ios::cur);
                    swq1.read((char*)&(us[is].real()),sizeof(us[is].real()));
                    swq2.seekg((ltt-1)*4,ios::cur);
                    swq2.read((char*)&(us[is].imag()),sizeof(us[is].imag()));
                    swq3.seekg((ltt-1)*4,ios::cur);
                    swq3.read((char*)&(ur[is].real()),sizeof(ur[is].real()));
                    swq4.seekg((ltt-1)*4,ios::cur);
                    swq4.read((char*)&(ur[is].imag()),sizeof(ur[is].imag()));
                 }
               
               for(is=0;is<ns;is++)
                 {
                    for(ir=0;ir<nr;ir++)
                       {
                          itmp=is*nr+ir;
                          us1[ir]=us[itmp];
                       } 
     
                    for(ir=0;ir<nr;ir++)
                       {
                          for(ix=0;ix<nr;ix++)  
                            {  
                              itmp=ir*ns+ix;
                              ur1[ix]=ur[itmp];
                            }

                          for(ix=0;ix<nr;ix++)
                            {
                               gr[ix].real()=0.0;
                               gr[ix].imag()=-omega[it]*tau[is][ix];
                               gr[ix]=exp(gr[ix]);
                            }

                          for(ix=0;ix<ns;ix++)
                            {
                               gs[ix].real()=0.0;
                               gs[ix].imag()=-omega[it]*tau[ir][ix];
                               gs[ix]=exp(gs[ix]);
                            }  
                       
                           for(ix=0;ix<nr;ix++)
                            {
                               gr1[ix].real()=1.0+gr[ix].real();
                               gr1[ix].imag()=gr[ix].imag();
                               gs1[ix].real()=1.0+gs[ix].real();
                               gs1[ix].imag()=gs[ix].imag();
                            }


                          itmp=is*nr+ir;
                          srmg[itmp]=(0.0,0.0);                           
                          
                          if(ir>is)
                             {
                                if(ir-is<10)
                                   {
                                      for(ix=is;ix<ir+1;ix++)
                                        {
                                           gh=(0.0,0.0); 
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda); 
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(ir-is+1)); 
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(ir-is+1));
                                        } 
                                   } 
                                 else
                                  {
                                      for(ix=is;ix<is+int((ir-is)*st_beg);ix++)
                                           srmg[itmp]+=(0.0,0.0);
                                       for(ix=is+int((ir-is)*st_beg);ix<is+int((ir-is)*st_end);ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(int((ir-is)*st_end)-int((ir-is)*st_beg)));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(int((ir-is)*st_end)-int((ir-is)*st_beg)));
                                        }
                                       for(ix=is+int((ir-is)*st_end);ix<ir+1;ix++)
                                          srmg[ir]+=(0.0,0.0);
                                  }
                             }
                          else if(ir<is)
                             {
                                if(is-ir<10)
                                   {
                                      for(ix=ir;ix<is+1;ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(is-ir+1));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(is-ir+1));
                                        } 
                                   } 
                                 else
                                  {
                                      for(ix=ir;ix<ir+int((is-ir)*st_beg);ix++)
                                           srmg[itmp]+=(0.0,0.0);
                                       for(ix=ir+int((is-ir)*st_beg);ix<ir+int((is-ir)*st_end);ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(int((is-ir)*st_end)-int((is-ir)*st_beg)));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(int((is-ir)*st_end)-int((is-ir)*st_beg)));
                                        }
                                       for(ix=ir+int((is-ir)*st_end);ix<is+1;ix++)
                                          srmg[ir]+=(0.0,0.0);
                                  }
                             } 
                          else
                             {
                                gh=(0.0,0.0);
                                gh=(gr[is]+gs[is])/(gr1[is]*gs1[is]+lamda);
                                srmg[itmp]+=gh*us1[is]*ur1[is];
                             } 
                       }//end of ir
                 }//end of is             

                amp1=0.0;
                amp2=0.0;

                for(ix=0;ix<nsr;ix++)
                 {
                    amp1+=sqrt(us[ix].real()*us[ix].real()+us[ix].imag()*us[ix].imag());
                    amp2+=sqrt(srmg[ix].real()*srmg[ix].real()+srmg[ix].imag()*srmg[ix].imag()); 
                 }

                for(ix=0;ix<nsr;ix++)
                   {
                      srmg[ix].real()*=amp1/(amp2+0.000001);
                      srmg[ix].imag()*=amp1/(amp2+0.000001);
                   }

//calculate the initial L2-error                
                l2norm=0.0;
                for(ix=0;ix<nsr;ix++)
                  {
                    if(us[ix].real()==0.0&&us[ix].imag()==0.0)
                      l2norm+=0.0;
                    else
                      l2norm+=((pow((us[ix].real()-srmg[ix].real()),2)+pow((us[ix].imag()-srmg[ix].imag()),2))/(pow(us[ix].real(),2)+pow(us[ix].imag(),2)+0.000001));
                  }
                l2norm=l2norm/nsr;

                it_count=0;

                cout<<it<<"th frequency slice initial value done... Error is===="<<l2norm<<endl;

                if(l2norm<err)
                   {
                      for(ix=0;ix<nsr;ix++)
                         {
                            swq5.write((char*)&srmg[is].real(),sizeof(srmg[is].real()));
                            swq6.write((char*)&srmg[is].imag(),sizeof(srmg[is].imag())); 
                         }
                      for(ix=0;ix<nsr;ix++)
                         desrmg[ix]=us[ix]-srmg[ix];
                      for(ix=0;ix<nsr;ix++)
                         {
                            swq7.write((char*)&desrmg[is].real(),sizeof(desrmg[is].real()));
                            swq8.write((char*)&desrmg[is].imag(),sizeof(desrmg[is].imag()));
                         }

                      cout<<it<<"th desrmg done... "<<"Iteration Time is 1"<<endl;
                   }
//8888 else
                else
                 {
//6666 while       
                   it_count=1;
                   while(it_count<=itmax)
                     {
                       l2norm=0.0;  
//                       cout<<it_count<<"   "<<"222"<<endl; 

                      for(is=0;is<ns;is++)
                        for(ir=0;ir<nr;ir++)
                         { 
                           srand((int)time(0));
                           dtau[is][ir]=(float)(rand()/float(RAND_MAX))*(dtau_max-dtau_min)+dtau_min;
                         }
                      for(is=0;is<ns;is++)
                        for(ir=0;ir<nr;ir++)
                           tau1[is][ir]=tau[is][ir]+dtau[is][ir];
/*

                      cout<<tau1[100][100]<<endl;
                     
                       for(is=0;is<ns;is++)
                        for(ir=0;ir<nr;ir++)
                           cout<<dtau[is][ir]<<"  ";
                       cout<<endl;
*/

                      for(is=0;is<ns;is++)
                        {
                          for(ir=0;ir<nr;ir++)
                            {
                               itmp=is*nr+ir;
                               us1[ir]=us[itmp];
                            }

                          for(ir=0;ir<nr;ir++)
                            {
                              for(ix=0;ix<ns;ix++)
                               { 
                                 itmp=ir*ns+ix;
                                 ur1[ix]=ur[itmp];
                               }

                              for(ix=0;ix<nr;ix++)
                               {
                                 gr[ix].real()=0.0;
                                 gr[ix].imag()=-omega[it]*tau1[is][ix];
                                 gr[ix]=exp(gr[ix]);
                               }

                              for(ix=0;ix<ns;ix++)
                               {
                                 gs[ix].real()=0.0;
                                 gs[ix].imag()=-omega[it]*tau1[ir][ix];
                                 gs[ix]=exp(gs[ix]);
                               }

                              for(ix=0;ix<ns;ix++)
                               {
                                 gr1[ix].real()=1.0+gr[ix].real();
                                 gr1[ix].imag()=gr[ix].imag();
                                 gs1[ix].real()=1.0+gs[ix].real();
                                 gs1[ix].real()=gs[ix].imag();
                               }
/*
                             if(is==100&&ir==100)
                              {
                                cout<<gs[100]<<endl;
                                cout<<gs1[100]<<endl;
                              }
*/
                          itmp=is*nr+ir;
                          srmg[itmp]=(0.0,0.0);

                          if(ir>is)
                             {
                                if(ir-is<10)
                                   {
                                      for(ix=is;ix<ir+1;ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(ir-is+1));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(ir-is+1));
                                        }
                                   }
                                 else
                                  {
                                      for(ix=is;ix<is+int((ir-is)*st_beg);ix++)
                                           srmg[itmp]+=(0.0,0.0);
                                       for(ix=is+int((ir-is)*st_beg);ix<is+int((ir-is)*st_end);ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(int((ir-is)*st_end)-int((ir-is)*st_beg)));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(int((ir-is)*st_end)-int((ir-is)*st_beg)));
                                        }
                                       for(ix=is+int((ir-is)*st_end);ix<ir+1;ix++)
                                          srmg[ir]+=(0.0,0.0);
                                  }
                             }
                          else if(ir<is)
                             {
                                if(is-ir<10)
                                   {
                                      for(ix=ir;ix<is+1;ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(is-ir+1));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(is-ir+1));
                                        }
                                   }
                                 else
                                  {
                                      for(ix=ir;ix<ir+int((is-ir)*st_beg);ix++)
                                           srmg[itmp]+=(0.0,0.0);
                                       for(ix=ir+int((is-ir)*st_beg);ix<ir+int((is-ir)*st_end);ix++)
                                        {
                                           gh=(0.0,0.0);
                                           gh=(gr[ix]+gs[ix])/(gr1[ix]*gs1[ix]+lamda);
                                           srmg[itmp].real()+=((gh*us1[ix]*ur1[ix]).real()/(int((is-ir)*st_end)-int((is-ir)*st_beg)));
                                           srmg[itmp].imag()+=((gh*us1[ix]*ur1[ix]).imag()/(int((is-ir)*st_end)-int((is-ir)*st_beg)));
                                        }
                                       for(ix=ir+int((is-ir)*st_end);ix<is+1;ix++)
                                          srmg[ir]+=(0.0,0.0);
                                  }
                             }
                          else
                             {
                                gh=(0.0,0.0);
                                gh=(gr[is]+gs[is])/(gr1[is]*gs1[is]+lamda);
                                srmg[itmp]+=gh*us1[is]*ur1[is];
                             }
                       }//end of ir
                 }//end of is             

                amp1=0.0;
                amp2=0.0;
                for(ix=0;ix<nsr;ix++)
                   {
                       amp1+=sqrt(us[ix].real()*us[ix].real()+us[ix].imag()*us[ix].imag());
                       amp2+=sqrt(srmg[ix].real()*srmg[ix].real()+srmg[ix].imag()*srmg[ix].imag());                         
                   }

                for(ix=0;ix<nsr;ix++)
                 {
                     srmg[ix].real()*=amp1/(amp2+0.000001);
                     srmg[ix].imag()*=amp1/(amp2+0.000001);

                 }     
                cout<<us[10000]<<endl;  
                cout<<srmg[10000]<<endl; 
                cout<<(pow((us[10000].real()-srmg[10000].real()),2)+pow((us[10000].imag()-srmg[10000].imag()),2))/(pow(us[ix].real(),2)+pow(us[ix].imag(),2)+0.000001)<<endl;

//calculate the initial L2-error                
                l2norm=0.0;
                for(ix=0;ix<nsr;ix++)
                  {
                    if(us[ix].real()==0.0&&us[ix].imag()==0.0)
                      l2norm+=0.0;
                    else
                      l2norm+=((pow((us[ix].real()-srmg[ix].real()),2)+pow((us[ix].imag()-srmg[ix].imag()),2))/(pow(us[ix].real(),2)+pow(us[ix].imag(),2)+0.000001));
                  }
                l2norm=l2norm/nsr;
                cout<<it<<"th frequency ... Iteration Time is===="<<it_count<<",error is===="<<l2norm<<endl;

                if(l2norm<err)
                    {
                      for(ix=0;ix<nsr;ix++)
                        {
                            swq5.write((char*)&srmg[is].real(),sizeof(srmg[is].real()));
                            swq6.write((char*)&srmg[is].imag(),sizeof(srmg[is].imag()));
                        }
                      for(ix=0;ix<nsr;ix++)
                         desrmg[ix]=us[ix]-srmg[ix];
                      for(ix=0;ix<nsr;ix++)
                        {
                            swq7.write((char*)&desrmg[is].real(),sizeof(desrmg[is].real()));
                            swq8.write((char*)&desrmg[is].imag(),sizeof(desrmg[is].imag()));
                        }

                      cout<<it<<"th frequency ... Iteration Time is===="<<it_count<<endl;

                      it_count=itmax; 
                    }
                 else
                   { 
                    it_count+=1;
                    if(it_count==itmax)
                     {
                      for(ix=0;ix<nsr;ix++)
                        {
                            swq5.write((char*)&srmg[is].real(),sizeof(srmg[is].real()));
                            swq6.write((char*)&srmg[is].imag(),sizeof(srmg[is].imag()));
                        }   
                      for(ix=0;ix<nsr;ix++)
                         desrmg[ix]=us[ix]-srmg[ix];
                      for(ix=0;ix<nsr;ix++)
                        {
                            swq7.write((char*)&desrmg[is].real(),sizeof(desrmg[is].real()));
                            swq8.write((char*)&desrmg[is].imag(),sizeof(desrmg[is].imag()));
                        }

                      cout<<it<<"th frequency desrmg done... Iteration Time is===="<<itmax<<endl;

                     }
                  } 
             }//6666 while end 
         }//end of else 8888  
      }//end of 9999else  
   } //end of 7777for

    swq1.close();
    swq2.close();
    swq3.close();
    swq4.close();
    swq5.close();
    swq6.close();
    swq7.close();
    swq8.close();
 
    cout<<" Positive Frequency Done!"<<endl;
    cout<<" Negative Frequency Starts===="<<endl;
     
   ifstream swq55;
   swq55.open(fn5,ios::binary);
   if(!swq55)
      {
         cout<<"cannot open "<<fn5<<endl;
         abort();
      }

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

   for(is=0;is<nsr;is++)
     {
        for(it=0;it<ltt;it++)
         {
           srmgf[it]=(0.0,0.0);
           desrmgf[it]=(0.0,0.0);
         }
        swq55.seekg(is*4,ios::beg);
        swq55.read((char*)&(srmgf[0].real()),sizeof(srmgf[0].real()));
        swq66.seekg(is*4,ios::beg);
        swq66.read((char*)&(srmgf[0].imag()),sizeof(srmgf[0].imag()));
        swq77.seekg(is*4,ios::beg);
        swq77.read((char*)&(desrmgf[0].real()),sizeof(desrmgf[0].real()));
        swq88.seekg(is*4,ios::beg);
        swq88.read((char*)&(desrmgf[0].imag()),sizeof(desrmgf[0].imag()));                 

        for(it=1;it<lt+1;it++)
          {
             swq55.seekg((nsr-1)*4,ios::cur);
             swq55.read((char*)&(srmgf[it].real()),sizeof(srmgf[it].real()));
             swq66.seekg((nr-1)*4,ios::cur);
             swq66.read((char*)&(srmgf[it].imag()),sizeof(srmgf[it].imag()));
             swq77.seekg((nsr-1)*4,ios::cur);
             swq77.read((char*)&(desrmgf[it].real()),sizeof(desrmgf[it].real()));
             swq88.seekg((nr-1)*4,ios::cur);
             swq88.read((char*)&(desrmgf[it].imag()),sizeof(desrmgf[it].imag()));
          }

       for(it=lt+1;it<ltt;it++)
          {
             srmgf[it].real()=srmgf[ltt-it].real();
             srmgf[it].imag()=-srmgf[ltt-it].imag();
             desrmgf[it].real()=desrmgf[ltt-it].real();
             desrmgf[it].imag()=-desrmgf[ltt-it].imag();
          }

       if((in1==NULL)||(out1==NULL))
       cout<<"memory insufficient"<<endl;
       else
         {
            for(it=0;it<ltt;it++)
              {
                 in1[it][0]=srmgf[it].real();
                 in1[it][1]=srmgf[it].imag();
              }
         }

       fftwf_execute(p1);

       for(it=0;it<ltt;it++)
         srmgt[it]=out1[it][0];

       for(it=0;it<lt;it++)
          swq9.write((char*)&(srmgt[it]),sizeof(srmgt[it]));

       if((in2==NULL)||(out2==NULL))
       cout<<"memory insufficient"<<endl;
       else
         {
            for(it=0;it<ltt;it++)
              {
                 in2[it][0]=desrmgf[it].real();
                 in2[it][1]=desrmgf[it].imag();
              }
         }

       fftwf_execute(p2);

       for(it=0;it<ltt;it++)
         desrmgt[it]=out2[it][0];

       for(it=0;it<lt;it++)
          swq10.write((char*)&(desrmgt[it]),sizeof(desrmgt[it]));

     }

    swq55.close();
    swq66.close();
    swq77.close();
    swq88.close();
    swq9.close();
    swq10.close();

    cout<<"ALL DONE!"<<endl; 

  return 0; 

}










