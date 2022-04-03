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

int xf_local_tau_p(complex<float> **dxf, complex<float> **dpf, complex <float> *base, float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax)
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
                 base[ix].imag()=p[ip]*omega[it]*x[ix];
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

int cal_space_spectrum(complex<float> **dxf, complex<float> **psd, complex<float> **dxf1, complex <float> *base,complex <float> *base1, complex<float>*base_acs, complex<float>**acs, float *omega, float *kx, float *x, int nx, int lt, int ifmin, int ifmax, float *psd_amp_tmp)
{
    int ix, it, ip, ix1;

   ofstream swq88;
   swq88.open("./psd_amp_snap.dat",ios::binary);
   if(!swq88)
     {
       cout<<"cannot open "<<endl;
       return 0;
     }

   for(ip=0;ip<nx;ip++)
     for(it=0;it<lt;it++)
       psd[ip][it]=(0.0,0.0);

   for(ix=0;ix<nx;ix++)
     for(it=0;it<lt;it++)
      {      
        dxf1[ix][it].real()=dxf[ix][it].real();
        dxf1[ix][it].imag()=-dxf[ix][it].imag();
      }
       

   for(it=ifmin;it<ifmax;it++)
//   for(it=95;it<96;it++)
     {
       //calculate the acs of the it the frequency slice
       for(ix=0;ix<nx;ix++)
         for(ix1=0;ix1<nx;ix1++)
           acs[ix][ix1]=(0.0,0.0);

       for(ix=0;ix<nx;ix++)
        {
          for(ix1=0;ix1<nx;ix1++)
            acs[ix][ix1]=dxf[ix][it]*dxf1[ix1][it];
        }     
 
        for(ip=0;ip<nx;ip++)
           {
             for(ix=0;ix<nx;ix++)
               {
                 base[ix].real()=0.0;
                 base[ix].imag()=kx[ip]*x[ix];
        	 base[ix]=exp(base[ix]);
        	 base1[ix].real()=base[ix].real();
        	 base1[ix].imag()=-base[ix].imag();
               }
//vector multipied by matrix
              for(ix=0;ix<nx;ix++)
                base_acs[ix]=(0.0,0.0);

              for(ix=0;ix<nx;ix++)
                {
                  for(ix1=0;ix1<nx;ix1++)  
                    base_acs[ix]+=base[ix1]*acs[ix1][ix];
                }

              for(ix=0;ix<nx;ix++)
                 psd[ip][it]+=base_acs[ix]*base1[ix];
           }
/*
         for(ip=0;ip<nx;ip++)
            cout<<ip+1<<"    "<<psd[ip][it]<<endl;
*/
         for(ip=0;ip<nx;ip++)
           psd_amp_tmp[ip]=sqrt(pow(psd[ip][it].real(),2)+pow(psd[ip][it].imag(),2));

          for(ip=0;ip<nx;ip++)
               swq88.write((char*)&psd_amp_tmp[ip],sizeof(psd_amp_tmp[ip])); 

     }  

   for(ip=0;ip<nx;ip++)
     for(it=lt/2+1;it<lt;it++)
       {
          psd[ip][it].real()=psd[ip][lt-it].real();
          psd[ip][it].imag()=-psd[ip][lt-it].imag();
       }
   return 0;

}

int CLEAN(complex<float> **upf_new, float rou, int ite_max, float **psd_amp, complex<float> **upf, complex <float> *base2,  complex <float> *base22, complex <float> *base3, complex <float> *base4, float *amp_tmp, complex<float> **dxf, complex<float> **psd, complex<float> **dxf1,complex<float> **acsk,complex<float> *base_acsk,  float *omega, float *p, float *x, int nx, int np, int lt, int ifmin, int ifmax)
{
    int ir,ix,ix1,ip,it, ite,ipmax,sign;
    int ite_tmp;

    for(ip=0;ip<np;ip++)
      for(it=0;it<lt;it++)
        {
          upf_new[ip][it]=(0.0,0.0);
          psd_amp[ip][it]=0.0;
        }

    for(ix=0;ix<nx;ix++)
      for(it=0;it<lt;it++)
       {
          dxf1[ix][it].real()=dxf[ix][it].real();
          dxf1[ix][it].imag()=-dxf[ix][it].imag();
       }

    for(it=ifmin;it<ifmax;it++)
//    for(it=79;it<80;it++)
       {
         cout<<it+1<<" frequency slices  begins..."<<endl;

        for(ix=0;ix<nx;ix++)
         {
           for(ix1=0;ix1<nx;ix1++)
            acsk[ix][ix1]=dxf[ix][it]*dxf1[ix1][it];
         }



           for(ite=0;ite<ite_max;ite++)
            {
              for(ip=0;ip<np;ip++)
                psd_amp[ip][it]=sqrt(pow(psd[ip][it].real(),2)+pow(psd[ip][it].imag(),2));  

//              for(ip=0;ip<np;ip++)
//                 cout<<ip<<"    "<<psd_amp[ip][it]<<endl;

              ite_tmp=ite;

                ipmax=0;
                for(ip=0;ip<np-1;ip++)
                  {
                     if(psd_amp[ip+1][it]>psd_amp[ipmax][it])
                       ipmax=ip+1;                        
                  }

                cout<<"   ===="<<ite+1<<" ,  "<<ipmax+1<<endl;

                upf_new[ipmax][it]=upf[ipmax][it];
   
                for(ix=0;ix<nx;ix++)
                 {
                   base2[ix].real()=0.0;
                   base2[ix].imag()=-p[ipmax]*omega[it]*x[ix];
                   base2[ix]=exp(base2[ix]);
                 }
               
                for(ip=0;ip<np;ip++)
                  {
                    base4[ip]=(0.0,0.0); 

                    for(ix=0;ix<nx;ix++)
                     {
                       base3[ix].real()=0.0;
                       base3[ix].imag()=p[ip]*omega[it]*x[ix];
                       base3[ix]=exp(base3[ix]);
                     }
                    for(ix=0;ix<nx;ix++)
                       base4[ip]+=base2[ix]*base3[ix];

                     amp_tmp[ip]=rou*1.0*sqrt(pow(psd[ipmax][it].real(),2)+pow(psd[ipmax][it].imag(),2))/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2));
                  }


                    for(ip=0;ip<np;ip++)
                      { 
//                        cout<<psd[ip][it].real()<<endl;
//                        cout<<psd[ip][it].imag()<<endl;


                        psd[ip][it].real()-=1.0*psd[ipmax][it].real()*rou/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2));
                        psd[ip][it].imag()-=1.0*psd[ipmax][it].imag()*rou/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2));

                        
//                        cout<<2.0*psd[ipmax][it].real()*rou/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2))<<endl;
//                        cout<<2.0*psd[ipmax][it].imag()*rou/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2))<<endl;
                        

                      }



/*
               sign=0;
               for(ip=0;ip<np;ip++)
                 {
                    if(psd_amp[ip][it]>amp_tmp[ip])
                        sign+=1;
                    else
                        sign-=1;  
                 }
                
               cout<<"   "<<ite+1<<", sign===="<<sign<<endl;

//               if(sign==np)
               if(sign>=np/2)
                 {
                    cout<<"   ===="<<ite+1<<" ,  "<<ipmax+1<<endl;

                    for(ip=0;ip<np;ip++)
                      { 
                        psd[ip][it].real()-=2.0*psd[ipmax][it].real()*rou/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2));
                        psd[ip][it].imag()-=2.0*psd[ipmax][it].imag()*rou/nx/nx*(pow(base4[ip].real(),2)+pow(base4[ip].imag(),2));
                      }
                 }
               else
                ite+=ite_max+999;
*/



            }
           cout<<"   ===== iteration time is===="<<ite_tmp+1<<endl;

       }

   for(ip=0;ip<np;ip++)
      for(it=lt/2+1;it<lt;it++)
       {
          psd[ip][it].real()=psd[ip][lt-it].real();
          psd[ip][it].imag()=-psd[ip][lt-it].imag();
       }
   for(ip=0;ip<np;ip++)
      for(it=lt/2+1;it<lt;it++)
        {
          upf_new[ip][it].real()=upf_new[ip][lt-it].real();
          upf_new[ip][it].imag()=-upf_new[ip][lt-it].imag();
        }

   return 0;

}

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256];
   int ns, nr, lt;
   float dx, dt;
   int nx_h, np_h;
   float fmin, fmax,theta_max,theta_min,v, rou, thre;
   int ite_max;

   int is,ir,ix,ip,it;

   ifstream swq;
   swq.open("xf_local_slant_stack_filter.par");
   if(!swq)
     {
        cout<<"Cannot open xf_local_slant_stack_filter.par"<<endl;
        return 0;
     } 
   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>ns>>nr>>lt>>dx>>dt>>nx_h>>np_h>>fmin>>fmax>>theta_max>>theta_min>>v>>rou>>ite_max>>thre;
   swq.close();

   int np=2*np_h+1;
   int nx=2*nx_h+1;

   float *omega;
   omega=alloc1float(lt);
   for(it=0;it<lt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*lt);
   for(it=lt/2+1;it<lt;it++)
      omega[it]=-2*pai*1000/(2*dt)+2*pai*(it-lt/2)*1000/(dt*lt);

   float *kx;
   kx=alloc1float(nx);
   for(ix=0;ix<nx/2+1;ix++)
      kx[ix]=2*pai*float(ix)/float(dx*nx);
   for(ix=nx/2+1;ix<nx;ix++)
      kx[ix]=-2*pai*1.0/float(2*dx)+2*pai*float(ix-nx/2)/float(dx*nx);

   cout<<"dk,kmax"<<2*pai/float(dx*nx)<<" , "<<2*pai*float(nx/2)/float(dx*nx)<<endl;
 
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
      cout<<ip<<p[ip]<<endl;
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
     cout<<ix<<x[ix]<<endl;
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

   complex<float> **uf1;
   uf1=alloc2complex(lt,nx);

   complex<float> **upf;
   upf=alloc2complex(lt,np);
   complex<float> **upf_new;
   upf_new=alloc2complex(lt,np);

   complex<float> **psd;
   psd=alloc2complex(lt,nx);
       
   float **upf_amp;
   upf_amp=alloc2float(lt,np);

   float **psd_amp;
   psd_amp=alloc2float(lt,nx);

   float *psd_amp_tmp;
   psd_amp_tmp=alloc1float(nx);

   complex<float> *base;
   base=alloc1complex(nx); 

   complex<float> *base1;
   base1=alloc1complex(nx); 

   complex<float> *base2;
   base2=alloc1complex(nx); 
   complex<float> *base22;
   base22=alloc1complex(nx); 

   complex<float> *base3;
   base3=alloc1complex(nx); 

   complex<float> *base4;
   base4=alloc1complex(np); 

   complex<float> *base_acs;
   base_acs=alloc1complex(nx); 

   complex<float> *base_acsk;
   base_acsk=alloc1complex(nx); 

   complex<float> **acs;
   acs=alloc2complex(nx,nx);

   complex<float> **acsk;
   acsk=alloc2complex(nx,nx);

   float *amp_tmp;
   amp_tmp=alloc1float(np);

   float **sem;
   sem=alloc2float(lt,np);

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
   ofstream swq4;
   swq4.open(fn4,ios::binary);
   if(!swq4)
     {
       cout<<"cannot open "<<fn4<<endl;
       return 0;
     }
   ofstream swq5;
   swq5.open(fn5,ios::binary);
   if(!swq5)
     {
       cout<<"cannot open "<<fn4<<endl;
       return 0;
     }
   ofstream swq6;
   swq6.open(fn6,ios::binary);
   if(!swq6)
     {
       cout<<"cannot open "<<fn6<<endl;
       return 0;
     }

   ifstream swq7;
   swq7.open(fn7,ios::binary);
   if(!swq7)
     {
       cout<<"cannot open "<<fn7<<endl;
       return 0;
     }

   ofstream swq8;
   swq8.open(fn8,ios::binary);
   if(!swq8)
     {
       cout<<"cannot open "<<fn8<<endl;
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
        for(ir=100;ir<101;ir++)
          {

            swq7.seekg(0,ios::beg);
            for(ix=0;ix<ir;ix++)
               swq7.seekg(np*lt*4,ios::cur);

            for(ip=0;ip<np;ip++)
              for(it=0;it<lt;it++)
               swq7.read((char*)&sem[ip][it],sizeof(sem[ip][it]));

             for(ix=0;ix<nx;ix++)
                for(it=0;it<lt;it++)
                  uf[ix][it]=uexf[ir+ix][it];           

             xf_local_tau_p(uf, upf, base, omega, p, x,  nx,  np, lt, ifmin, ifmax);

             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                  upf_amp[ip][it]=sqrt(pow(upf[ip][it].real(),2)+pow(upf[ip][it].imag(),2));

             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                swq3.write((char*)&upf_amp[ip][it],sizeof(upf_amp[ip][it]));

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
         
            cout<<ir+1<<" trace tp done!"<<endl;     
 
           cal_space_spectrum(uf, psd, uf1, base,base1, base_acs, acs, omega, kx, x, nx,  lt,  ifmin,  ifmax, psd_amp_tmp);

             for(ip=0;ip<nx;ip++)
               for(it=0;it<lt;it++)
                  psd_amp[ip][it]=sqrt(pow(psd[ip][it].real(),2)+pow(psd[ip][it].imag(),2));

             for(ip=0;ip<nx;ip++)
               for(it=0;it<lt;it++)
                swq4.write((char*)&psd_amp[ip][it],sizeof(psd_amp[ip][it]));

             cout<<ir+1<<" trace psd calculation done!"<<endl;     
/*
           CLEAN(upf_new, rou, ite_max, psd_amp, upf, base2,base22, base3, base4, amp_tmp,uf, psd, uf1,acsk, base_acsk,  omega, p, x, nx,  np,  lt,  ifmin,  ifmax);

             for(ip=0;ip<np;ip++)
               for(it=0;it<lt;it++)
                swq5.write((char*)&psd_amp[ip][it],sizeof(psd_amp[ip][it]));

             if((in2==NULL)||(out2==NULL))
               cout<<"memory insufficient"<<endl;
             else
              { 
               for(ip=0;ip<np;ip++)
                 {
                   for(it=0;it<lt;it++)
                    {
                     in2[it][0]=upf_new[ip][it].real();
                     in2[it][1]=upf_new[ip][it].imag();
                    }

                   fftwf_execute(p2);

                   for(it=0;it<lt;it++)
                      upt[ip][it]=out2[it][0];
                 }
               }

            for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               swq6.write((char*)&upt[ip][it],sizeof(upt[ip][it]));

            for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               {
                  if(sem[ip][it]<=thre)
                    upt[ip][it]=0.0; 
               }

            for(ip=0;ip<np;ip++)
             for(it=0;it<lt;it++)
               swq8.write((char*)&upt[ip][it],sizeof(upt[ip][it]));

             cout<<ir+1<<" trace CLEAN calculation done!"<<endl;     
*/
          }

       cout<<is+1<<" shot done!"<<endl; 
       
     }

   swq1.close();
   swq2.close();
   swq3.close();
   swq4.close();


   cout<<"All Done!"<<endl;

  return 0; 
}














