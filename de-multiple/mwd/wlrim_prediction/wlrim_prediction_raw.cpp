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

int green_function( complex <float> *sfk, complex <float> **gfk,  complex <float> **gf, float **g, float *omega, float *kx, float sdep, float rdep, float ds, float dr, int nr, int lt, int ltt, int boul, int bour,  int nrb, float v, float theta, int ifmin, int ifmax, float wd)
{
    float kz;
    complex <float> ps1;
    complex <float> ps2;
    complex <float> ps3;
    complex <float> ps4;
    complex <float> ps5;
    complex <float> ps6;
    complex <float> ps7;

    int ir,is,it;

//zeroing the arrays
    for(ir=0;ir<nrb;ir++) 
       for(it=0;it<ltt;it++)
         {
            gfk[ir][it]=(0.0,0.0);
         }

     for(ir=0;ir<nr;ir++)
       for(it=0;it<ltt;it++)
          {
            gf[ir][it]=(0.0,0.0);
            g[ir][it]=0.0;
          }

    fftwf_complex *in3,*out3;
    fftwf_plan p3;
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    p3=fftwf_plan_dft_1d(nrb,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);

    fftwf_complex *in4,*out4;
    fftwf_plan p4;
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p4=fftwf_plan_dft_1d(ltt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

     for(ir=0;ir<nrb;ir++)
        {
           for(it=0;it<ltt/2+1;it++)
              {
                 if(it<ifmin)
                   {
                      gfk[ir][it].real()=0.0;
                      gfk[ir][it].imag()=0.0;
                   }
                 else if(it>ifmax&&it<ltt/2+1)
                    {
                      gfk[ir][it].real()=0.0;
                      gfk[ir][it].imag()=0.0;
                   }
                 else
                   {
                      if(pow(omega[it]/v,2)-pow(kx[ir],2)<0.0)
                        {
                           gfk[ir][it].real()=0.0;
                           gfk[ir][it].imag()=0.0;
                        }
                      else if(pow(omega[it]/v,2)-pow(kx[ir],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                        {
                             kz=sqrt(pow(omega[it]/v,2)-pow(kx[ir],2));

                             ps1.real()=0.0;               
                             ps1.imag()=-2*wd*kz;
                             ps1=exp(ps1);

                             gfk[ir][it]=sfk[ir]*ps1; 
                        }
                       else
                        {
                           gfk[ir][it].real()=0.0;
                           gfk[ir][it].imag()=0.0;
                        } 
                   }

                }
        }


//ifft

     cout<<"IFFT Begins..."<<endl;
 
       for(it=0;it<ltt/2+1;it++)
         {   
           for(ir=0;ir<nrb;ir++)
              {   
                 in3[ir][0]=gfk[ir][it].real();
                 in3[ir][1]=gfk[ir][it].imag();
              }   
           fftwf_execute(p3); 

           for(ir=0;ir<nr;ir++)
             {   
               gf[ir][it].real()=out3[ir+boul][0];
               gf[ir][it].imag()=out3[ir+boul][1];
             }   
         }   

     for(ir=0;ir<nr;ir++)
       for(it=ltt/2+1;it<ltt;it++)
          {
             gf[ir][it].real()=gf[ir][ltt-it].real();
             gf[ir][it].imag()=-gf[ir][ltt-it].imag();
          }

        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in4[it][0]=gf[ir][it].real();
                 in4[it][1]=gf[ir][it].imag();
              }
           fftwf_execute(p4);

           for(it=0;it<lt;it++)
             g[ir][it]=out4[it][0];
         }


    return 0;

}

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256],fn5[256],fn6[256],fn7[256],fn8[256],fn9[256],fn10[256],fn100[256];
   int ns,nr,lt,ltt,boul,bour, ifmin, ifmax, flag,nrb, n;
   float ds,dr,dt,sdep, rdep,theta, v, fmax, fmin, lamda,wd, min,max;
 
   ifstream swq;
   swq.open("wlrim_prediction.par");
   swq>>fn1>>fn2>>fn3>>fn4>>fn5>>fn6>>fn7>>fn8>>fn9>>fn10>>fn100>>ns>>nr>>lt>>dt>>ds>>dr>>sdep>>rdep>>v>>theta>>boul>>bour>>fmax>>fmin>>wd>>min>>max;
   swq.close();

   cout<<"Fna of Green's Function in x-f domain Real Parts is===="<<fn1<<endl; 
   cout<<"Fna of Green's Function in x-f domain Imaginary Parts is===="<<fn2<<endl; 
   cout<<"Fna of Green's Function in x-t domain is===="<<fn3<<endl; 
   cout<<"Fna of Input CSG in x-t domain is===="<<fn4<<endl; 
   cout<<"Fna of Input CRG in x-t domain is===="<<fn5<<endl; 
   cout<<"Fna of CSG in x-f domain Real Parts is===="<<fn6<<endl; 
   cout<<"Fna of CSG in x-f domain Imaginary Parts is===="<<fn7<<endl; 
   cout<<"Fna of CRG in x-f domain Real Parts is===="<<fn8<<endl; 
   cout<<"Fna of CRG in x-f domain Imaginary Parts is===="<<fn9<<endl; 
   cout<<"Fna of output WLRIM in x-t domain is===="<<fn10<<endl; 
   cout<<"Fna of output SRM in x-t domain is===="<<fn100<<endl; 
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

   int it,it1,is,ir,is1,ir1,is2,ir2;

   float *omega;
   omega=alloc1float(ltt);
   float *kx;
   kx=alloc1float(nrb);

   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   ifmin=int(fmin*dt*ltt/1000);
   ifmax=int(fmax*dt*ltt/1000)+1;

   cout<<"Totally "<<ifmax-ifmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

   for(is=0;is<nrb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dr*nrb);
   for(is=nrb/2+1;is<nrb;is++)
      kx[is]=-2*pai*1.0/float(2*dr)+2*pai*float(is-nrb/2)/float(dr*nrb);
   
   complex<float> *sf;
   sf=alloc1complex(nrb); 

   complex<float> *sfk;
   sfk=alloc1complex(nrb);

   complex<float> **gfk;
   gfk=alloc2complex(ltt,nrb);

   complex<float> **gf;
   gf=alloc2complex(ltt,nr);

   float **g;
   g=alloc2float(ltt,nr);

   float **usraw;
   usraw=alloc2float(lt,nr);
   float **urraw;
   urraw=alloc2float(lt,ns);

   float **us;
   us=alloc2float(ltt,nr);
   float **ur;
   ur=alloc2float(ltt,ns);

   complex<float> **usf;
   usf=alloc2complex(ltt,nr);
   complex<float> **urf;
   urf=alloc2complex(ltt,ns);
    
   complex<float> **srmf;
   srmf=alloc2complex(ltt,nr);
   float **srm;
   srm=alloc2float(ltt,nr);

   complex<float> **imf;
   imf=alloc2complex(ltt,nr);
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

  ofstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
       cout<<"cannot open "<<fn1<<endl;

  ofstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
       cout<<"cannot open "<<fn2<<endl;

  ofstream swq3;
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

  ofstream swq6;
  swq6.open(fn6,ios::binary);
  if(!swq6)
       cout<<"cannot open "<<fn6<<endl;

  ofstream swq7;
  swq7.open(fn7,ios::binary);
  if(!swq7)
       cout<<"cannot open "<<fn7<<endl;

  ofstream swq8;
  swq8.open(fn8,ios::binary);
  if(!swq8)
       cout<<"cannot open "<<fn8<<endl;

  ofstream swq9;
  swq9.open(fn9,ios::binary);
  if(!swq9)
       cout<<"cannot open "<<fn9<<endl;

  for(is=0;is<ns;is++)
    {
       for(ir=0;ir<nrb;ir++)
          sf[ir]=(0.0,0.0);

       sf[is+boul].real()=1.0;
       sf[is+boul].imag()=0.0;
        
       for(ir=0;ir<nrb;ir++)
          {
            in2[ir][0]=sf[ir].real();
            in2[ir][1]=sf[ir].imag();
          }
       
        fftwf_execute(p2);

        for(ir=0;ir<nrb;ir++)
         {
           sfk[ir].real()=out2[ir][0]/nrb;
           sfk[ir].imag()=out2[ir][1]/nrb;
         }

       green_function(sfk, gfk, gf, g, omega, kx, sdep, rdep, ds, dr, nr, lt, ltt, boul, bour, nrb, v, theta, ifmin, ifmax, wd);

       for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
          { 
           swq1.write((char*)&gf[ir][it].real(),sizeof(gf[ir][it].real()));
           swq2.write((char*)&gf[ir][it].imag(),sizeof(gf[ir][it].imag()));
          }
  
       for(ir=0;ir<nr;ir++)
        for(it=0;it<lt;it++)
           swq3.write((char*)&g[ir][it],sizeof(g[ir][it]));

        cout<<is+1<<" shot Green's Function calculation done!"<<endl;
    }
 
   swq1.close();
   swq2.close();
   swq3.close();
 
   cout<<"Green's Function Done!"<<endl;

   for(is=0;is<ns;is++)
     {
        for(ir=0;ir<nr;ir++)
          for(it=0;it<ltt;it++)
            us[ir][it]=0.0;

        for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++)
            swq4.read((char*)&usraw[ir][it],sizeof(usraw[ir][it]));

         for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++)
            us[ir][it]=usraw[ir][it];

        for(ir=0;ir<nr;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in5[it][0]=us[ir][it];
                 in5[it][1]=0.0;
              }
           fftwf_execute(p5);

           for(it=0;it<lt;it++)
             {
               usf[ir][it].real()=out5[it][0];
               usf[ir][it].imag()=out5[it][1];
             }
         }

       for(ir=0;ir<nr;ir++)
        for(it=0;it<ltt;it++)
          { 
           swq6.write((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real()));
           swq7.write((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag()));
          }

      }

   for(is=0;is<nr;is++)
     {
        for(ir=0;ir<ns;ir++)
          for(it=0;it<ltt;it++)
            ur[ir][it]=0.0;

        for(ir=0;ir<ns;ir++)
          for(it=0;it<lt;it++)
            swq5.read((char*)&urraw[ir][it],sizeof(urraw[ir][it]));

         for(ir=0;ir<ns;ir++)
          for(it=0;it<lt;it++)
            ur[ir][it]=urraw[ir][it];

        for(ir=0;ir<ns;ir++)
         {
           for(it=0;it<ltt;it++)
              {
                 in5[it][0]=ur[ir][it];
                 in5[it][1]=0.0;
              }

           fftwf_execute(p5);

           for(it=0;it<lt;it++)
             {
               urf[ir][it].real()=out5[it][0];
               urf[ir][it].imag()=out5[it][1];
             }
         }

       for(ir=0;ir<ns;ir++)
        for(it=0;it<ltt;it++)
          { 
           swq8.write((char*)&urf[ir][it].real(),sizeof(urf[ir][it].real()));
           swq9.write((char*)&urf[ir][it].imag(),sizeof(urf[ir][it].imag()));
          }

      }

   swq4.close();
   swq5.close();
   swq6.close();
   swq7.close();
   swq8.close();
   swq9.close();



  ifstream swq11;
  swq11.open(fn1,ios::binary);
  if(!swq11)
       cout<<"cannot open "<<fn1<<endl;

  ifstream swq22;
  swq22.open(fn2,ios::binary);
  if(!swq22)
       cout<<"cannot open "<<fn2<<endl;

  ifstream swq66;
  swq66.open(fn6,ios::binary);
  if(!swq66)
       cout<<"cannot open "<<fn6<<endl;

  ifstream swq77;
  swq77.open(fn7,ios::binary);
  if(!swq77)
       cout<<"cannot open "<<fn7<<endl;
 
  ifstream swq88;
  swq88.open(fn8,ios::binary);
  if(!swq88)
       cout<<"cannot open "<<fn8<<endl;

  ifstream swq99;
  swq99.open(fn9,ios::binary);
  if(!swq99)
       cout<<"cannot open "<<fn9<<endl;

  ofstream swq10;
  swq10.open(fn10,ios::binary);
  if(!swq10)
       cout<<"cannot open "<<fn10<<endl;

  ofstream swq100;
  swq100.open(fn100,ios::binary);
  if(!swq100)
       cout<<"cannot open "<<fn100<<endl;

  cout<<"Surface-Related Multiple Prediction Begins..."<<endl;

  swq66.seekg(0,ios::beg);
  swq77.seekg(0,ios::beg);

//  for(is=0;is<ns;is++)  //A Loop
  for(is=0;is<1;is++)  //A Loop
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
           srmf[ir][it]=(0.0,0.0);

       cout<<"   "<<is+1<<" shot Surface-Related Multiple Prediction Begins..."<<endl;

       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
            {
              swq66.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real()));
              swq77.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag()));
            }        

       swq88.seekg(0,ios::beg);
       swq99.seekg(0,ios::beg);
        
       for(ir=0;ir<nr;ir++)  //E Loop
         {
           cout<<"      "<<is+1<<" shot, "<<ir+1<<" trace, Surface-Related Multiple Prediction Begins..."<<endl;

           for(is1=0;is1<ns;is1++)
            for(it=0;it<ltt;it++)
             {
               swq88.read((char*)&urf[is1][it].real(),sizeof(urf[is1][it].real()));
               swq99.read((char*)&urf[is1][it].imag(),sizeof(urf[is1][it].imag()));
             }    
    
             if(is<ir)
               {
                 for(it=ifmin;it<ifmax;it++)
                  {
                    if(ir-is<10)
                     { 
                       for(ir1=is;ir1<ir+1;ir1++) //F Loop
                         srmf[ir][it]+=usf[ir1][it]*urf[ir1][it];                           
                     }
                    else
                     { 
                       for(ir1=is+int((ir-is)*min);ir1<is+int((ir-is)*max);ir1++) //F Loop
                         srmf[ir][it]+=usf[ir1][it]*urf[ir1][it];                           
                     }
                  }                
               }             

              else
               {
                 for(it=ifmin;it<ifmax;it++)
                  {
                    if(is-ir<10)
                     {
                      for(ir1=ir;ir1<is+1;ir1++) //F Loop
                       srmf[ir][it]+=usf[ir1][it]*urf[ir1][it];    
                     }
                    else
                     {
                       for(ir1=ir+int((is-ir)*min);ir1<ir+int((is-ir)*max);ir1++) //F Loop
                         srmf[ir][it]+=usf[ir1][it]*urf[ir1][it];    
                     } 
                  }                
               }             

           for(it=ltt/2+1;it<ltt;it++)
            {
              srmf[ir][it].real()=srmf[ir][ltt-it].real();               
              srmf[ir][it].imag()=-srmf[ir][ltt-it].imag();               
            }  

           for(it=0;it<ltt;it++)
              {
                 in6[it][0]=srmf[ir][it].real();
                 in6[it][1]=srmf[ir][it].imag();
              }
           fftwf_execute(p6);

           for(it=0;it<ltt;it++)
             srm[ir][it]=out6[it][0];

           for(it=0;it<lt;it++)
             swq100.write((char*)&srm[ir][it],sizeof(srm[ir][it]));

           cout<<"      "<<is+1<<" shot, "<<ir+1<<" trace Surface-Related Multiple Prediction Done!"<<endl;

         } //E Loop

       cout<<"   "<<is+1<<" shot Surface-Related Multiple Prediction Done!"<<endl;

    } //A Loop

   cout<<"Surface-Related Multiple Prediction Done!"<<endl;

//   return 0;

  cout<<"Water-Bottom-Related Internal Multiple Prediction Begins..."<<endl;

  swq66.seekg(0,ios::beg);
  swq77.seekg(0,ios::beg);

//  for(is=0;is<ns;is++)  //A Loop
  for(is=0;is<1;is++)  //A Loop
    {
       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
           imf[ir][it]=(0.0,0.0);

       cout<<"   "<<is+1<<" shot Internal Multiple Prediction Begins..."<<endl;

       for(ir=0;ir<nr;ir++)
         for(it=0;it<ltt;it++)
            {
              swq66.read((char*)&usf[ir][it].real(),sizeof(usf[ir][it].real()));
              swq77.read((char*)&usf[ir][it].imag(),sizeof(usf[ir][it].imag()));
            }        

       swq88.seekg(0,ios::beg);
       swq99.seekg(0,ios::beg);
        
       for(ir=0;ir<nr;ir++)  //E Loop
         {
           cout<<"      "<<is+1<<" shot, "<<ir+1<<" trace,  Internal Multiple Prediction Begins..."<<endl;

           for(is1=0;is1<ns;is1++)
            for(it=0;it<ltt;it++)
             {
               swq88.read((char*)&urf[is1][it].real(),sizeof(urf[is1][it].real()));
               swq99.read((char*)&urf[is1][it].imag(),sizeof(urf[is1][it].imag()));
             }    
    
             if(is<ir)
               {
//                for(ir1=0;ir1<nr;ir1++) //F Loop
//                for(ir1=is;ir1<ir+1;ir1++) //F Loop
                for(ir1=is+int(ir-is+1)*min;ir1<is+int(ir-is+1)*max;ir1++) //F Loop
                  {
                     swq11.seekg(0,ios::beg);
                     swq22.seekg(0,ios::beg);

                     for(ir2=0;ir2<ir1;ir2++)
                       {
                         swq11.seekg(ns*ltt*4,ios::cur);
                         swq22.seekg(ns*ltt*4,ios::cur);
                       }

                     for(ir2=0;ir2<ns;ir2++)
                       for(it1=0;it1<ltt;it1++)
                         {
                           swq11.read((char*)&gf[ir2][it1].real(),sizeof(gf[ir2][it1].real()));
                           swq22.read((char*)&gf[ir2][it1].imag(),sizeof(gf[ir2][it1].imag()));
                         }

                      for(ir2=0;ir2<ns;ir2++)
                       for(it1=0;it1<ltt;it1++)
                          gf[ir2][it1].imag()*=-1.0;                    

                     for(it=ifmin;it<ifmax;it++)
                       {//1111 
//                         for(ir2=0;ir2<nr;ir2++) //G Loop
//                          for(ir2=ir1;ir2<ir+1;ir2++) //G Loop
                          for(ir2=ir1+int(ir-ir1+1)*min;ir2<ir1+int(ir-ir1+1)*max;ir2++) //G Loop
                            imf[ir][it]+=usf[ir2][it]*urf[ir1][it]*gf[ir2][it];                           

//                            imf[ir][it].real()/=nr;
//                            imf[ir][it].imag()/=nr;

//                            imf[ir][it].real()/=(ir-ir1+1);
//                            imf[ir][it].imag()/=(ir-ir1+1);

                            imf[ir][it].real()/=(int(ir-ir1+1)*max-int(ir-ir1+1)*min);
                            imf[ir][it].imag()/=(int(ir-ir1+1)*max-int(ir-ir1+1)*min);
                       }//1111
                  }                
               }             

              else
               {
//                for(ir1=0;ir1<nr;ir1++) //F Loop
//                for(ir1=ir;ir1<is+1;ir1++) //F Loop
                for(ir1=ir+int(is-ir+1)*min;ir1<ir+int(is-ir+1)*max;ir1++) //F Loop
                  {
                     swq11.seekg(0,ios::beg);
                     swq22.seekg(0,ios::beg);

                     for(ir2=0;ir2<ir1;ir2++)
                       {
                         swq11.seekg(ns*ltt*4,ios::cur);
                         swq22.seekg(ns*ltt*4,ios::cur);
                       }

                     for(ir2=0;ir2<ns;ir2++)
                       for(it1=0;it1<ltt;it1++)
                         {
                           swq11.read((char*)&gf[ir2][it1].real(),sizeof(gf[ir2][it1].real()));
                           swq22.read((char*)&gf[ir2][it1].imag(),sizeof(gf[ir2][it1].imag()));
                         }

                      for(ir2=0;ir2<ns;ir2++)
                       for(it1=0;it1<ltt;it1++)
                          gf[ir2][it1].imag()*=-1.0;                    

                     for(it=ifmin;it<ifmax;it++)
                       {//1111 
//                         for(ir2=0;ir2<nr;ir2++) //G Loop
//                         for(ir2=ir1;ir2<is+1;ir2++) //G Loop
                          for(ir2=ir1+int(is-ir1+1)*min;ir2<ir1+int(is-ir1+1)*max;ir2++) //G Loop
                            imf[ir][it]+=usf[ir2][it]*urf[ir1][it]*gf[ir2][it];                           

//                            imf[ir][it].real()/=nr;
//                            imf[ir][it].imag()/=nr;

//                           imf[ir][it].real()/=(is-ir1+1);
//                           imf[ir][it].imag()/=(is-ir1+1);

                            imf[ir][it].real()/=(int(is-ir1+1)*max-int(is-ir1+1)*min);
                            imf[ir][it].imag()/=(int(is-ir1+1)*max-int(is-ir1+1)*min);
                       }//1111
                  }                
               }             

           for(it=ltt/2+1;it<ltt;it++)
            {
              imf[ir][it].real()=imf[ir][ltt-it].real();               
              imf[ir][it].imag()=-imf[ir][ltt-it].imag();               
            }  

           for(it=0;it<ltt;it++)
              {
                 in6[it][0]=imf[ir][it].real();
                 in6[it][1]=imf[ir][it].imag();
              }
           fftwf_execute(p6);

           for(it=0;it<ltt;it++)
             im[ir][it]=out6[it][0];

           for(it=0;it<lt;it++)
             swq10.write((char*)&im[ir][it],sizeof(im[ir][it]));

           cout<<"      "<<is+1<<" shot, "<<ir+1<<" trace Internal Multiple Prediction Done!"<<endl;

         } //E Loop

       cout<<"   "<<is+1<<" shot Internal Multiple Prediction Done!"<<endl;

    } //A Loop

   cout<<"Internal Multiple Prediction Done!"<<endl;

   swq11.close();
   swq22.close();
   swq66.close();
   swq77.close();
   swq88.close();
   swq99.close();
   swq10.close();

   return 0; 

}
































