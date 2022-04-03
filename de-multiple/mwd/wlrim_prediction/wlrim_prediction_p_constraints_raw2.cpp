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
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn11[256],fninmcg[256],fnoutmcg[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb, n,irmcg,np;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda,wd, min,max;
 
   ifstream swq;
   swq.open("wlrim_prediction_p_constraints.par");
   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>fn11>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>wd>>min>>max>>irmcg>>np;
   swq.close();

   cout<<"Fna of CSG Green's Function in p-f domain Real Parts is===="<<fn1<<endl; 
   cout<<"Fna of CSG Green's Function in p-f domain Imaginary Parts is===="<<fn2<<endl; 
   cout<<"Fna of CSG in p-f domain Real Parts is===="<<fn3<<endl; 
   cout<<"Fna of CSG in p-f domain Imaginary Parts is===="<<fn4<<endl; 
   cout<<"Fna of CRG in p-f domain Real Parts is===="<<fn5<<endl; 
   cout<<"Fna of CRG in p-f domain Imaginary Parts is===="<<fn6<<endl; 
   cout<<"Fna of CSG Green's Function Semblance in p-t domain is===="<<fn7<<endl; 
   cout<<"Fna of CRG Green's Function Semblance in p-t domain is===="<<fn8<<endl; 
   cout<<"Fna of Inner MCG in x-t domain is===="<<fn9<<endl; 
   cout<<"Fna of Outer MCG in x-t domain is===="<<fn10<<endl; 
   cout<<"Fna of output WLRIM in x-t domain is===="<<fn11<<endl; 
   cout<<"No. of shots is===="<<ns<<endl;
   cout<<"No. of traces per shot is===="<<nr<<endl;
   cout<<"lt and dt is===="<<lt<<" , "<<dt<<"ms"<<endl;
   cout<<"Source and Receiver Interval are===="<<ds<<"m , "<<dr<<"m"<<endl;
   cout<<"Source and Receiver Depth are===="<<sdep<<"m , "<<rdep<<"m"<<endl;
   cout<<"Water Velocity is===="<<v<<endl;
   cout<<"Maximum Theta is===="<<theta<<endl;
   cout<<"Left and Right Boundary are===="<<boul<<" , "<<bour<<endl;  
   cout<<"Maximum and Minimum Frequency are===="<<fmax<<"Hz , "<<fmin<<"Hz"<<endl;
   cout<<"Water Bottom Depth is===="<<wd<<endl;
   cout<<"Minimum and Maximum threshold parameter are===="<<min<<" , "<<max<<endl;

   ltt=2*lt;

   nrb=boul+nr+bour;

   int *ps;
   ps=alloc1int(ns);
   int *pr;
   pr=alloc1int(nr);


   int it,it1,is,ir,is1,ir1,is2,ir2,ir3,ir4,ip;

   int max_ipr,max_ips,tmp1,tmp2;
   float maxsemr,maxsems;

   ifmax=int(fmax*dt*ltt/1000)+1;
   ifmin=int(fmin*dt*ltt/1000);

   cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

   complex<float> **uspf;
   uspf=alloc2complex(ltt,np*nr);
   complex<float> **urpf;
   urpf=alloc2complex(ltt,np*ns);

   float **sems;
   sems=alloc2float(lt,np*nr);
   float **semr;
   semr=alloc2float(lt,np*nr);

   complex<float> **gpf;
   gpf=alloc2complex(ltt,np*nr);

   complex<float> **inmcgpf;
   inmcgpf=alloc2complex(ltt,np*nr);
   float **inmcgpt;
   inmcgpt=alloc2float(ltt,np*nr);
   float **inmcgt;
   inmcgt=alloc2float(ltt,nr);
   float **outmcgt;
   outmcgt=alloc2float(ltt,nr);

   complex<float> **impf;
   impf=alloc2complex(ltt,nr*np);
   float **impt;
   impt=alloc2float(ltt,nr*np);
   float **im;
   im=alloc2float(ltt,nr);

   fftwf_complex *in2,*out2;
   fftwf_plan p2;
   in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
   out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
   p2=fftwf_plan_dft_1d(nrb,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

   fftwf_complex *in5,*out5;
   fftwf_plan p5;
   in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p5=fftwf_plan_dft_1d(ltt,in5,out5,FFTW_FORWARD,FFTW_MEASURE);

   fftwf_complex *in6,*out6;
   fftwf_plan p6;
   in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
   p6=fftwf_plan_dft_1d(ltt,in6,out6,FFTW_BACKWARD,FFTW_MEASURE);

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
 
  ifstream swq5;
  swq5.open(fn5,ios::binary);
  if(!swq5)
       cout<<"cannot open "<<fn5<<endl;

  ifstream swq6;
  swq6.open(fn6,ios::binary);
  if(!swq6)
       cout<<"cannot open "<<fn6<<endl;

  ifstream swq7;
  swq7.open(fn7,ios::binary);
  if(!swq7)
       cout<<"cannot open "<<fn7<<endl;

  ifstream swq8;
  swq8.open(fn8,ios::binary);
  if(!swq8)
       cout<<"cannot open "<<fn8<<endl;
/*
  ofstream swq9;
  swq9.open(fn9,ios::binary);
  if(!swq9)
       cout<<"cannot open "<<fn9<<endl;

  ofstream swq10;
  swq10.open(fn10,ios::binary);
  if(!swq10)
       cout<<"cannot open "<<fn10<<endl;
*/
  ofstream swq11;
  swq11.open(fn11,ios::binary);
  if(!swq11)
       cout<<"cannot open "<<fn11<<endl;

  cout<<"Water-Bottom-Related Internal Multiple Prediction Begins..."<<endl;

  swq1.seekg(0,ios::beg);
  swq2.seekg(0,ios::beg);
  swq3.seekg(0,ios::beg);
  swq4.seekg(0,ios::beg);
  swq5.seekg(0,ios::beg);
  swq6.seekg(0,ios::beg);

//  for(is=0;is<ns;is++)  //A Loop
  for(is=0;is<1;is++)  //A Loop
    {
       for(ir=0;ir<nr*np;ir++)
         for(it=0;it<ltt;it++)
           impf[ir][it]=(0.0,0.0);

       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
           im[ir][it]=0.0;

       cout<<"   "<<is+1<<" shot Internal Multiple Prediction Begins..."<<endl;

       for(ir=0;ir<np*nr;ir++)
         for(it=0;it<ltt;it++)
            {
              swq3.read((char*)&uspf[ir][it].real(),sizeof(uspf[ir][it].real()));
              swq4.read((char*)&uspf[ir][it].imag(),sizeof(uspf[ir][it].imag()));
            }        

       swq5.seekg(0,ios::beg);
       swq6.seekg(0,ios::beg);
        
//       for(ir=0;ir<nr;ir++)  //E Loop
       for(ir=250;ir<260;ir++)  //E Loop
         {
           cout<<"      "<<is+1<<" shot, "<<ir+1<<" trace,  Internal Multiple Prediction Begins..."<<endl;

           for(ir1=0;ir1<ir;ir1++)
              {
                 swq5.seekg(np*ns*ltt*4,ios::cur);
                 swq6.seekg(np*ns*ltt*4,ios::cur);
              }

           for(is1=0;is1<ns*np;is1++)
            for(it=0;it<ltt;it++)
             {
               swq5.read((char*)&urpf[is1][it].real(),sizeof(urpf[is1][it].real()));
               swq6.read((char*)&urpf[is1][it].imag(),sizeof(urpf[is1][it].imag()));
             }    

             if(is<ir)
               {
                 for(ir1=is+int((ir-is)*min);ir1<is+int((ir-is)*max);ir1++) //F Loop
//                 for(ir1=is;ir1<ir+1;ir1++) //F Loop
//                 for(ir1=0;ir1<nr;ir1++) //F Loop
                  {
                     swq1.seekg(0,ios::beg);
                     swq2.seekg(0,ios::beg);
                     swq7.seekg(0,ios::beg);
                     for(ir2=0;ir2<ir1;ir2++)
                       {
                         swq1.seekg(nr*np*ltt*4,ios::cur);
                         swq2.seekg(nr*np*ltt*4,ios::cur);
                         swq7.seekg(nr*np*lt*4,ios::cur);
                       }

                     for(ir2=0;ir2<nr*np;ir2++)
                       for(it1=0;it1<ltt;it1++)
                         {
                           swq1.read((char*)&gpf[ir2][it1].real(),sizeof(gpf[ir2][it1].real()));
                           swq2.read((char*)&gpf[ir2][it1].imag(),sizeof(gpf[ir2][it1].imag()));
                         }

                      for(ir2=0;ir2<nr*np;ir2++)
                       for(it1=0;it1<ltt;it1++)
                          gpf[ir2][it1].imag()*=-1.0;                    

                     for(ir2=0;ir2<nr*np;ir2++)
                       for(it1=0;it1<lt;it1++)
                         swq7.read((char*)&semr[ir2][it1],sizeof(semr[ir2][it1]));

                     for(ir2=0;ir2<nr;ir2++)
                       {
                          max_ipr=0;
                          maxsemr=0.0;
                          for(ir3=ir2*np;ir3<(ir2+1)*np;ir3++)
                            {
                               for(it1=0;it1<lt;it1++)
                                 {
                                    if(semr[ir3][it1]>maxsemr)
                                       {
                                          maxsemr=semr[ir3][it1];
                                          max_ipr=ir3;
                                       }
                                 }
                            }

                         pr[ir2]=max_ipr-ir2*np;

                         for(ir3=ir2*np;ir3<max_ipr;ir3++)
                            for(it1=0;it1<ltt;it1++)
                              gpf[ir3][it1]=(0.0,0.0);

                         for(ir3=max_ipr+1;ir3<(ir2+1)*np;ir3++)
                            for(it1=0;it1<ltt;it1++)
                              gpf[ir3][it1]=(0.0,0.0);

                       } 

                          for(ir2=ir1+int((ir-ir1)*min);ir2<ir1+int((ir-ir1)*max);ir2++) //G Loop
//                          for(ir2=ir1;ir2<ir+1;ir2++) //G Loop
//                          for(ir2=0;ir2<nr;ir2++) //G Loop
                            {
                              swq8.seekg(0,ios::beg);
                              for(ir3=0;ir3<ir2;ir3++)
                                swq8.seekg(nr*np*lt*4,ios::cur);
                            
                              for(ir3=0;ir3<nr*np;ir3++)
                                for(it1=0;it1<lt;it1++)
                                  swq8.read((char*)&sems[ir3][it1],sizeof(sems[ir3][it1]));

                              for(ir3=0;ir3<nr;ir3++)
                                {
                                   max_ips=0;
                                   maxsems=0.0;
                                  for(ir4=ir3*np;ir4<(ir3+1)*np;ir4++)
                                   {
                                     for(it1=0;it1<lt;it1++)
                                      {
                                        if(sems[ir4][it1]>maxsems)
                                         {
                                           maxsems=sems[ir4][it1];
                                           max_ips=ir4;
                                         }
                                      }
                                   }

                                  ps[ir3]=max_ips-ir3*np;
                                } 

                              tmp1=pr[ir2];
                              tmp2=ps[ir1]; 

                              for(ip=0;ip<np;ip++) 
                               for(it=0;it<ltt;it++)
                                impf[ir*np+ip][it]=(0.0,0.0);

                              for(it=ifmin;it<ifmax;it++)
                                impf[ir*np+tmp1][it]+=uspf[ir2*np+tmp1][it]*urpf[ir1*np+tmp2][it]*gpf[ir2*np+tmp1][it];

                            }
                  }                
               }             

             else
               {
                 for(ir1=ir+int((is-ir)*min);ir1<ir+int((is-ir)*max);ir1++) //F Loop
//                 for(ir1=ir;ir1<is+1;ir1++) //F Loop
//                 for(ir1=0;ir1<nr;ir1++) //F Loop
                  {
                     swq1.seekg(0,ios::beg);
                     swq2.seekg(0,ios::beg);
                     swq7.seekg(0,ios::beg);

                     for(ir2=0;ir2<ir1;ir2++)
                       {
                         swq1.seekg(nr*np*ltt*4,ios::cur);
                         swq2.seekg(nr*np*ltt*4,ios::cur);
                         swq7.seekg(nr*np*lt*4,ios::cur);
                       }

                     for(ir2=0;ir2<nr*np;ir2++)
                       for(it1=0;it1<ltt;it1++)
                         {
                           swq1.read((char*)&gpf[ir2][it1].real(),sizeof(gpf[ir2][it1].real()));
                           swq2.read((char*)&gpf[ir2][it1].imag(),sizeof(gpf[ir2][it1].imag()));
                         }

                      for(ir2=0;ir2<nr*np;ir2++)
                       for(it1=0;it1<ltt;it1++)
                          gpf[ir2][it1].imag()*=-1.0;                    

                     for(ir2=0;ir2<nr*np;ir2++)
                       for(it1=0;it1<lt;it1++)
                         swq7.read((char*)&semr[ir2][it1],sizeof(semr[ir2][it1]));

   
                     for(ir2=0;ir2<nr;ir2++)
                       {
                          max_ipr=0;
                          maxsemr=0.0;
                          for(ir3=ir2*np;ir3<(ir2+1)*np;ir3++)
                            {
                               for(it1=0;it1<lt;it1++)
                                 {
                                    if(semr[ir3][it1]>maxsemr)
                                       {
                                          maxsemr=semr[ir3][it1];
                                          max_ipr=ir3;
                                       }
                                 }
                            }

                         pr[ir2]=max_ipr-ir2*np;

                         for(ir3=ir2*np;ir3<max_ipr;ir3++)
                            for(it1=0;it1<ltt;it1++)
                              gpf[ir3][it1]=(0.0,0.0);

                         for(ir3=max_ipr+1;ir3<(ir2+1)*np;ir3++)
                            for(it1=0;it1<ltt;it1++)
                              gpf[ir3][it1]=(0.0,0.0);

                       } 

                          for(ir2=ir1+int((is-ir1)*min);ir2<ir+int((is-ir1)*max);ir2++) //G Loop
//                          for(ir2=ir1;ir2<is+1;ir2++) //G Loop
//                          for(ir2=0;ir2<nr;ir2++) //G Loop
                            {
                              swq8.seekg(0,ios::beg);
                              
                              for(ir3=0;ir3<ir2;ir3++)
                                swq8.seekg(nr*np*lt*4,ios::cur);
                            
                              for(ir3=0;ir3<nr*np;ir3++)
                                for(it1=0;it1<lt;it1++)
                                  swq8.read((char*)&sems[ir3][it1],sizeof(sems[ir3][it1]));

                              for(ir3=0;ir3<nr;ir3++)
                                {
                                   max_ips=0;
                                   maxsems=0.0;
                                  for(ir4=ir3*np;ir4<(ir3+1)*np;ir4++)
                                   {
                                     for(it1=0;it1<lt;it1++)
                                      {
                                        if(sems[ir4][it1]>maxsems)
                                         {
                                           maxsems=sems[ir4][it1];
                                           max_ips=ir4;
                                         }
                                      }
                                   }

                                  ps[ir3]=max_ips-ir3*np;

                                } 

                              tmp1=pr[ir2];
                              tmp2=ps[ir1]; 

                              for(ip=0;ip<np;ip++) 
                               for(it=0;it<ltt;it++)
                                impf[ir*np+ip][it]=(0.0,0.0);

                              for(it=ifmin;it<ifmax;it++)
                                 impf[ir*np+tmp1][it]+=uspf[ir2*np+tmp1][it]*urpf[ir1*np+tmp2][it]*gpf[ir2*np+tmp1][it];
                            }
                  }                
               }             

           for(ip=0;ip<np;ip++)
             for(it=ltt/2+1;it<ltt;it++)
              {
                impf[ir*np+ip][it].real()=impf[ir*np+ip][ltt-it].real();               
                impf[ir*np+ip][it].imag()=-impf[ir*np+ip][ltt-it].imag();               
              }  

           for(ip=0;ip<np;ip++)
            {
             for(it=0;it<ltt;it++)
              {
                 in6[it][0]=impf[ir*np+ip][it].real();
                 in6[it][1]=impf[ir*np+ip][it].imag();
              }
             fftwf_execute(p6);

             for(it=0;it<ltt;it++)
              impt[ir*np+ip][it]=out6[it][0]/ltt;
           }

          for(it=0;it<lt;it++)
           {
             for(ip=0;ip<np;ip++) 
              im[ir][it]+=impt[ir*np+ip][it];
           }

           for(it=0;it<lt;it++)
             swq11.write((char*)&im[ir][it],sizeof(im[ir][it]));

           cout<<"      "<<is+1<<" shot, "<<ir+1<<" trace Internal Multiple Prediction Done!"<<endl;

         } //E Loop

       cout<<"   "<<is+1<<" shot Internal Multiple Prediction Done!"<<endl;

    } //A Loop

   cout<<"Internal Multiple Prediction Done!"<<endl;

   swq1.close();
   swq2.close();
   swq3.close();
   swq4.close();
   swq5.close();
   swq6.close();
/*
   swq9.close();
   swq10.close();
*/
   swq11.close();

   return 0; 

}
































