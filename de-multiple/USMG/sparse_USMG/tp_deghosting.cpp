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

int inverse_tau_p(float **u,float *ufinal,int lt,int np)
{
    int i,j;
 
    for(i=0;i<lt;i++)
      ufinal[i]=0.0;

    for(i=0;i<lt;i++)
       {
           for(j=0;j<np;j++)
               ufinal[i]+=u[j][i];
       }
    return 0;
}

int main()
{
   char fna1[256],fna2[256],fna3[256],fna33[256],fna4[256],fna44[256];
   int ns,nr,np,lt,ltt;
   float dep,dt,fmin,fmax,theta_min,theta_max,v,sem_thre,lamda;
   int ifmin,ifmax;

   ifstream swq;
   swq.open("tp_deghosting.par");
   swq>>fna1>>fna2>>fna3>>fna33>>fna4>>fna44>>ns>>nr>>np>>lt>>dt>>dep>>fmin>>fmax>>theta_min>>theta_max>>v>>sem_thre>>lamda;
   swq.close();

   cout<<"Fna of input data in the local tau-p domain is===="<<fna1<<endl;
   cout<<"Fna of semblance of input data in the local tau-p domain is===="<<fna2<<endl;
   cout<<"Fna of output deghosting upgoing data in the local tau-p domain is===="<<fna3<<endl;
   cout<<"Fna of output deghosting downgoing data in the local tau-p domain is===="<<fna33<<endl;
   cout<<"Fna of output deghosting upgoing data in the x-t domain is===="<<fna4<<endl;
   cout<<"Fna of output deghosting downgoing data in the x-t domain is===="<<fna44<<endl;

   ltt=2*lt;

   int is,ir,ip,it;

   float *omega;
   omega=alloc1float(ltt);
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   ifmin=int(fmin*dt*ltt/1000);
   ifmax=int(fmax*dt*ltt/1000)+1;

    cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;
  
    float pmax,pmin,dp;
    pmin=sin(theta_min*2*pai/360.0)/v;
    pmax=sin(theta_max*2*pai/360.0)/v;
    dp=(pmax-pmin)/(np-1);

    float *p=new float [np];
    for(ip=0;ip<np;ip++)
        p[ip]=pmin+ip*dp;

    float **upt;
    upt=alloc2float(ltt,np);

    float **upt1;
    upt1=alloc2float(ltt,np);
 
   for(ip=0;ip<np;ip++)
      for(it=0;it<ltt;it++)
        upt[ip][it]=0.0;

    float **sem;
    sem=alloc2float(ltt,np);

    for(ip=0;ip<np;ip++)
      for(it=0;it<ltt;it++)
        sem[ip][it]=0.0;

    complex<float> **upf;
    upf=alloc2complex(ltt,np);

    complex<float> **upfu;
    upfu=alloc2complex(ltt,np);

    complex<float> **upfd;
    upfd=alloc2complex(ltt,np);

    complex<float> phase1;
    complex<float> phase2;
    complex<float> phase3;
    complex<float> phase4;

    float **upu;
    upu=alloc2float(lt,np);
    float **upd;
    upd=alloc2float(lt,np);
 
    float *uut;
    uut=alloc1float(lt);
    float *udt;
    udt=alloc1float(lt);

   fftwf_complex *in1,*out1;
    fftwf_plan p1;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

    ifstream swq1;
    swq1.open(fna1,ios::binary);
    if(!swq1)
      {
        cout<<"Cannot open "<<fna1<<endl;
        return 0;
      }

    ifstream swq2;
    swq2.open(fna2,ios::binary);
    if(!swq2)
      {
        cout<<"Cannot open "<<fna2<<endl;
        return 0;
      }

    ofstream swq3;
    swq3.open(fna3,ios::binary);
    if(!swq3)
      {
        cout<<"Cannot open "<<fna3<<endl;
        return 0;
      }
    ofstream swq33;
    swq33.open(fna33,ios::binary);
    if(!swq33)
      {
        cout<<"Cannot open "<<fna33<<endl;
        return 0;
      }

    ofstream swq4;
    swq4.open(fna4,ios::binary);
    if(!swq4)
      {
        cout<<"Cannot open "<<fna4<<endl;
        return 0;
      }

    ofstream swq44;
    swq44.open(fna44,ios::binary);
    if(!swq44)
      {
        cout<<"Cannot open "<<fna44<<endl;
        return 0;
      }
/* 
    ofstream swq5;
    swq5.open("./check_read.dat",ios::binary);
*/
   for(is=0;is<ns;is++)
      {
        cout<<is+1<<" shot begin ..."<<endl;

          for(ir=0;ir<nr;ir++)
             {
                  for(ip=0;ip<np;ip++)
                   for(it=0;it<ltt;it++) 
                     {
                        upfu[ip][it].real()=0.0;
                        upfu[ip][it].imag()=0.0;
                        upfd[ip][it].real()=0.0;
                        upfd[ip][it].imag()=0.0;
                     }

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      {
                          swq1.read((char*)&upt[ip][it],sizeof(upt[ip][it]));
                          swq2.read((char*)&sem[ip][it],sizeof(sem[ip][it]));
                      }

                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                      {
                        if(sem[ip][it]<sem_thre)
                           upt[ip][it]=0.0;
                      }
/*
                 for(ip=0;ip<np;ip++)
                   for(it=0;it<lt;it++)
                     swq5.write((char*)&upt[ip][it],sizeof(upt[ip][it]));
*/
                  for(ip=0;ip<np;ip++)
                    {
                      for(it=0;it<ltt;it++)
                        {
                          in1[it][0]=upt[ip][it];
                          in1[it][1]=0.0;
                        }

                       fftwf_execute(p1);

                       for(it=0;it<ltt;it++)
                         upf[ip][it]=out1[it][0]/ltt;
                    }

                  for(ip=0;ip<np;ip++)
                    {
                      for(it=ifmin;it<ifmax;it++)
                         {
                            if(pow(omega[it]/v,2)>pow(p[ip],2))
                              {
                                  phase1.real()=0.0;
                                  phase1.imag()=-1.0*2*fabs(dep)*sqrt(pow(omega[it]/v,2)-pow(p[ip],2));

//                                  phase1.real()=0.0;
//                                  phase1.imag()=-2*omega[it]*fabs(dep)*sqrt(1.0-pow(p[ip]*v,2));
    
                                  phase2=exp(phase1);

                                //Ghost Operator
                                  phase3.real()=-1.0*phase2.real();
                                  phase3.imag()=-1.0*phase2.imag();
                               
                                  phase4.real()=1.0+phase3.real()+lamda;
                                  phase4.imag()=phase3.imag();

                                  upfu[ip][it]=upf[ip][it]/phase4;
                                  upfd[ip][it]=phase3*upf[ip][it]/phase4;
                              } 
                         }
                    }

                 for(ip=0;ip<np;ip++)
                    for(it=ltt/2+1;it<ltt;it++)
                      {
                         upfu[ip][it].real()=upfu[ip][ltt-it].real();
                         upfu[ip][it].imag()=-upfu[ip][ltt-it].imag();
                         upfd[ip][it].real()=upfd[ip][ltt-it].real();
                         upfd[ip][it].imag()=-upfd[ip][ltt-it].imag();
                      }  

                 for(ip=0;ip<np;ip++)
                    {
                      for(it=0;it<ltt;it++)
                        {
                          in2[it][0]=upfu[ip][it].real();
                          in2[it][1]=upfu[ip][it].imag();
                        }

                       fftwf_execute(p2);

                       for(it=0;it<lt;it++)
                         upu[ip][it]=out2[it][0];
                    }
 
                 for(ip=0;ip<np;ip++)
                    {
                      for(it=0;it<ltt;it++)
                        {
                          in2[it][0]=upfd[ip][it].real();
                          in2[it][1]=upfd[ip][it].imag();
                        }

                       fftwf_execute(p2);

                       for(it=0;it<lt;it++)
                         upd[ip][it]=out2[it][0];
                    }

                 for(ip=0;ip<np;ip++) 
                   for(it=0;it<lt;it++)
                     {
                        swq3.write((char*)&upu[ip][it],sizeof(upu[ip][it]));  
                        swq33.write((char*)&upd[ip][it],sizeof(upd[ip][it]));  
                     }

                 for(it=0;it<lt;it++)
                    {
                       uut[it]=0.0;
                       udt[it]=0.0;
                    }

                 inverse_tau_p(upu,uut,lt,np);
                 inverse_tau_p(upd,udt,lt,np);

                 for(it=0;it<lt;it++) 
                   {
                      swq4.write((char*)&uut[it],sizeof(uut[it]));
                      swq44.write((char*)&udt[it],sizeof(udt[it]));
                   }

               cout<<"   "<<ir+1<<" trace done ..."<<endl;

             }

          cout<<is+1<<" shot done!"<<endl;

      } 

  return 0;

}








































