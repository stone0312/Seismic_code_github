using namespace std;
#include <iostream>
#include "math.h"
#include "stdio.h"
#include <fstream>
#include "stdlib.h"
#include <complex>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

int forward_xf_local_tau_p(complex<float> **dxf, complex<float> **dpf, complex <float> *base, float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax)
{
   int ix, it, ip;

   for(ip=0;ip<np;ip++)
     for(it=0;it<lt;it++)
       dpf[ip][it]=(0.0,0.0);      

   for(it=ifmin;it<ifmax;it++)
     {
//        cout<<"ifre===="<<it+1<<endl;

        for(ip=0;ip<np;ip++)
           { 
              for(ix=0;ix<nx;ix++)
               {
                 base[ix].real()=0.0;
                 base[ix].imag()=-p[ip]*omega[it]*x[ix];
                 base[ix]=exp(base[ix]);
               }
              for(ix=0;ix<nx;ix++)
                 dpf[ip][it]+=base[ix]*dxf[ix][it]; 
           }
     }
 
   for(ip=0;ip<np;ip++)
     for(it=lt/2+1;it<lt;it++)
       {
          dpf[ip][it].real()=dpf[ip][lt-it].real();
          dpf[ip][it].imag()=-dpf[ip][lt-it].imag();
       }


   return 0;

}

int inverse_xf_local_tau_p(complex<float> **dxf, complex<float> **dpf, complex <float> *base_T, float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax)
{
   int ix, it, ip;

   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++)
       dxf[ix][it]=(0.0,0.0);      

   for(it=ifmin;it<ifmax;it++)
     {
//        cout<<"ifre===="<<it+1<<endl;

          for(ix=0;ix<nx;ix++)
           { 
             for(ip=0;ip<np;ip++)
               {
                 base_T[ip].real()=0.0;
                 base_T[ip].imag()=p[ip]*omega[it]*x[ix];
                 base_T[ip]=exp(base_T[ip]);
               }
             for(ip=0;ip<np;ip++)
               dxf[ix][it]+=base_T[ip]*dpf[ip][it]; 
           }
     }
 
   for(ix=0;ix<nx;ix++)
     for(it=lt/2+1;it<lt;it++)
       {
          dxf[ix][it].real()=dxf[ix][lt-it].real();
          dxf[ix][it].imag()=-dxf[ix][lt-it].imag();
       }

   return 0;

}

int main()
{
   char fn1[256],fn2[256],fn3[256];
   int ns, nr, lt;
   float dx, dt;
   int nx_h, np_h;
   float fmin, fmax,theta_max,theta_min,v;

   int is,ir,ix,ip,it;

   ifstream swq;
   swq.open("xf_local_slant_stack.par");
   if(!swq)
     {
        cout<<"Cannot open xf_local_slant_stack.par"<<endl;
        return 0;
     } 
   swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>dx>>dt>>nx_h>>np_h>>fmin>>fmax>>theta_max>>theta_min>>v;
   swq.close();

   int np=2*np_h+1;
   int nx=2*nx_h+1;

   float *omega;
   omega=alloc1float(lt);
   for(it=0;it<lt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*lt);
   for(it=lt/2+1;it<lt;it++)
      omega[it]=-2*pai*1000/(2*dt)+2*pai*(it-lt/2)*1000/(dt*lt);
 
   int ifmin, ifmax;
   ifmin=int(fmin*dt*lt/1000);
   ifmax=int(fmax*dt*lt/1000)+1;

   cout<<ifmin<<" ,  "<<ifmax<<endl;
   cout<<"totally "<<ifmax-ifmin<<" frequency slices need to be calculated"<<endl;
  
   float *p;
   p=alloc1float(np);
   float dp;
   dp=(sin(theta_max*pai/180.0)/v-sin(theta_min*pai/180.0)/v)/(float)(np-1);

   p[np_h]=0.0;
   for(ip=np_h+1;ip<np;ip++)
      p[ip]=(ip-np_h)*dp;
   for(ip=0;ip<np_h;ip++)
      p[ip]=-1.0*p[2*np_h-ip];
/*
   for(ip=0;ip<np;ip++)
     cout<<ip<<"===="<<p[ip]<<endl;
   return 0;
*/

   float *x;
   x=alloc1float(nx);

   x[nx_h]=0.0;
   for(ix=(nx-1)/2+1;ix<nx;ix++)
      x[ix]=(ix-((nx-1)/2))*dx;
   for(ix=0;ix<nx_h;ix++)
      x[ix]=-1.0*x[2*nx_h-ix];
/*
    for(ix=0;ix<nx;ix++)
     cout<<ix<<"===="<<x[ix]<<endl;
   return 0;
*/

   float **uraw;
   uraw=alloc2float(lt,nr);

   complex<float> **urawf;
   urawf=alloc2complex(lt,nr);

   float **upt;
   upt=alloc2float(lt,nr*np);

   complex<float> **uexf;
   uexf=alloc2complex(lt,nr+nx+1);

   complex<float> **uf;
   uf=alloc2complex(lt,nx);

   complex<float> **uipf;
   uipf=alloc2complex(lt,nx);

   float **ui;
   ui=alloc2float(lt,nx);

   complex<float> **upf;
   upf=alloc2complex(lt,np);
       
   float **upf_amp;
   upf_amp=alloc2float(lt,np);

   complex<float> *base;
   base=alloc1complex(nx); 
   complex<float> *base_T;
   base_T=alloc1complex(np); 

   fftwf_complex *in1,*out1,*in2,*out2;
   in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
   out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
   in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
   out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

   fftwf_plan p1, p2;

    p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(lt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

   ifstream swq1;
   swq1.open(fn1,ios::binary);
   if(!swq1)
     {
       cout<<"cannot open "<<fn1<<endl;
       return 0;
     }
   ofstream swq2;
   swq2.open(fn2,ios::binary);
   if(!swq2)
     {
       cout<<"cannot open "<<fn2<<endl;
       return 0;
     }
   ofstream swq3;
   swq3.open(fn3,ios::binary);
   if(!swq3)
     {
       cout<<"cannot open "<<fn3<<endl;
       return 0;
     }

   cout<<"x-f domain local tau-p begin..."<<endl;

   for(is=0;is<ns;is++) 
     {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
            swq1.read((char*)&uraw[ir][it],sizeof(uraw[ir][it]));

       for(ir=0;ir<nr+nx+1;ir++)
         for(it=0;it<lt;it++)
              uexf[ir][it]=(0.0,0.0);

       if((in1==NULL)||(out1==NULL))
             cout<<"memory insufficient"<<endl;
          else
            {
               for(ir=0;ir<nr;ir++)
                 {
                   for(it=0;it<lt;it++)
                    { 
                     in1[it][0]=uraw[ir][it];
                     in1[it][1]=0.0;
                    }
                  
                   fftwf_execute(p1);

                   for(it=0;it<lt;it++)                  
                     {
                       urawf[ir][it].real()=out1[it][0]/(float)lt;
                       urawf[ir][it].imag()=out1[it][1]/(float)lt;
                     }
                 }
            }

        for(ir=nx_h;ir<nr+nx_h;ir++)  
           for(it=0;it<lt;it++)  
              uexf[ir][it]=urawf[ir-nx_h][it];

//        for(ir=0;ir<nr;ir++)
        for(ir=0;ir<1;ir++)
          {
             for(ix=0;ix<nx;ix++)
                for(it=0;it<lt;it++)
                  uf[ix][it]=uexf[ir+ix][it];           

             forward_xf_local_tau_p(uf, upf, base, omega, p, x,  nx,  np, lt, ifmin, ifmax);

             if((in2==NULL)||(out2==NULL))
               cout<<"memory insufficient"<<endl;
             else
              { 
               for(ip=0;ip<np;ip++)
                 {
                   for(it=0;it<lt;it++)
                    {
                     in2[it][0]=upf[ip][it].real();
                     in2[it][1]=upf[ip][it].imag();
                    }

                   fftwf_execute(p2);

                   for(it=0;it<lt;it++)
                      upt[ip][it]=out2[it][0];
                 }
               }

            for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               swq2.write((char*)&upt[ip][it],sizeof(upt[ip][it]));

            cout<<"2222"<<endl;
         
            inverse_xf_local_tau_p(uipf, upf, base_T, omega, p, x,  nx,  np, lt, ifmin, ifmax);

             if((in2==NULL)||(out2==NULL))
               cout<<"memory insufficient"<<endl;
             else
              { 
               for(ix=0;ix<nx;ix++)
                 {
                   for(it=0;it<lt;it++)
                    {
                     in2[it][0]=uipf[ix][it].real();
                     in2[it][1]=uipf[ix][it].imag();
                    }

                   fftwf_execute(p2);

                   for(it=0;it<lt;it++)
                      ui[ix][it]=out2[it][0];
                 }
               }
            
            for(ix=0;ix<nx;ix++)
              for(it=0;it<lt;it++)
                swq3.write((char*)&ui[ix][it],sizeof(ui[ix][it]));

            cout<<ir+1<<" trace tp done!"<<endl;     
 
        }

       cout<<is+1<<" shot tp done!"<<endl; 
       
     }

   swq1.close();
   swq2.close();
   swq3.close();


   cout<<"All Done!"<<endl;

  return 0; 
}













































