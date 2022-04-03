#include "iostream.h"
#include "fstream.h"
#include "math.h"
#include "complex.h"
#include "alloc.c"
#include "stdlib.h"
#include "fftw3.h"
#define pai 3.14159265

using namespace std;

int main()
{
  char fn1[256],fn2[256],fn3[256];
  int ns,nr,lt,itmin,itmax;
  float dt,sz,rz,dz,dsx,drx,ds,dr,lamda,fmin,fmax,v,kz,gh_amp;
  int it,is,ir,ix,irtmp,tmp;

  ifstream swq;
  swq.open("fk_deghosting.par");
  if(!swq)
     {
        cout<<"cannot open fk_deghosting.par"<<endl;
        abort();
     }

  swq>>fn1>>fn2>>fn3>>ns>>nr>>lt>>dt>>sz>>rz>>dz>>dsx>>drx>>ds>>dr>>lamda>>fmin>>fmax>>v;
  swq.close();

  cout<<"Fna of input data with Ghosts is===="<<fn1<<endl; 
  cout<<"Fna of output data without Ghosts is===="<<fn2<<endl;
  cout<<"No. of shots is===="<<ns<<endl;
  cout<<"No. of receivers per shot is===="<<nr<<endl;
  cout<<"No. of temperal samplings is===="<<lt<<endl;
  cout<<"Temperal Interval (ms) is===="<<dt<<endl;
  cout<<"Depth of shots and receiver are===="<<sz<<" , "<<rz<<" respectively"<<endl;
  cout<<"Spatial Interval in z direction is===="<<dz<<endl;
  cout<<"Spatial Intervals of shots and receivers are===="<<dsx<<" , "<<drx<<" respectively"<<endl;
  cout<<"Spatial Interval of shots and receivers for calculation are===="<<ds<<" , "<<dr<<" respectively"<<endl;
  cout<<"Minimum and Maximum frequency to be calculated are==== "<<fmin<<" , "<<fmax<<endl;
  cout<<"Velocity of water is===="<<v<<endl;

  fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4,*in5,*out5,*in6,*out6;
  
  in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
  in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);
  out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * lt);

  fftwf_plan p1,p2,p3,p4;

  p1=fftwf_plan_dft_1d(lt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
  p2=fftwf_plan_dft_1d(nr,in2,out2,FFTW_FORWARD,FFTW_MEASURE);
  p3=fftwf_plan_dft_1d(nr,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
  p4=fftwf_plan_dft_1d(lt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

  ifstream swq1;
  swq1.open(fn1,ios::binary);
  if(!swq1)
    {
       cout<<"cannot open "<<fn1<<endl;
       abort();
    } 

  ofstream swq2;
  swq2.open(fn2,ios::binary);
  if(!swq2)
    {
       cout<<"cannot open "<<fn2<<endl;
       abort();
    }
   ofstream swq22;
  swq22.open(fn3,ios::binary);
  if(!swq22)
    {
       cout<<"cannot open "<<fn3<<endl;
       abort();
    } 


 ofstream swq3;
 swq3.open("gh_amp.dat",ios::binary);

 ofstream swq4;
 swq4.open("orig_amp.dat",ios::binary);

  float **u;
  u=alloc2float(lt,nr);


  complex<float> **uf;
  uf=alloc2complex(lt,nr);

  complex<float> **ufk;
  ufk=alloc2complex(lt,nr);

  complex<float> *ps;
  ps=alloc1complex(nr);
 
  complex<float> *dg;
  dg=alloc1complex(nr); 

  complex<float> *ps2;
  ps2=alloc1complex(nr);

  complex<float> *dg2;
  dg2=alloc1complex(nr);

  complex<float> **ufk1;
  ufk1=alloc2complex(lt,nr);
  complex<float> **ufk2;
  ufk2=alloc2complex(lt,nr);

  complex<float> **uf1;
  uf1=alloc2complex(lt,nr);

  float **u1;
  u1=alloc2float(lt,nr);

  float **degh_amp_tmp;
  degh_amp_tmp=alloc2float(nr,lt/2+1);

  float **degh_amp;
  degh_amp=alloc2float(lt,nr);

  float **amp_orig;
  amp_orig=alloc2float(lt,nr);

  float **filt;
  filt=alloc2float(lt/2+1,nr);

  float *omega;
  omega=alloc1float(lt);

  float *kx;
  kx=alloc1float(nr);

//calculate omega and kx
  for(it=0;it<lt/2+1;it++)
     omega[it]=2*pai*it*1000/(dt*lt);
  for(it=lt/2+1;it<lt;it++)
     omega[it]=2*pai*(-1000/(2*dt)+(it-lt/2)*1000/(dt*lt));

  for(ix=0;ix<nr/2+1;ix++)
     kx[ix]=2*pai*(float)ix/(float)(dr*nr);    
  for(ix=nr/2+1;ix<nr;ix++) 
     kx[ix]=-2*pai/(2*dr)+2*pai*(float)(ix-nr/2)/(float)(dr*nr);

  itmin=int(fmin*dt*lt/1000);
  itmax=int(fmax*dt*lt/1000)+1; 
   
  cout<<"Totally "<<itmax-itmin+1<<" frequency slices needed to be calculated..."<<endl; 

  cout<<"Deghosting Starts..."<<endl;

  for(is=0;is<ns;is++)
    {
       cout<<is+1<<" shot deghosting starts..."<<endl;

       for(ir=0;ir<nr;ir++)
         for(it=0;it<lt;it++)
           swq1.read((char*)&u[ir][it],sizeof(u[ir][it]));

       for(ir=0;ir<nr;ir++)
          {
             for(it=0;it<lt;it++)
                {
                  in1[it][0]=u[ir][it];
                  in1[it][1]=0.0;
                }

             fftwf_execute(p1);

             for(it=0;it<lt;it++)
                {
                   uf[ir][it].real()=out1[it][0]/lt;
                   uf[ir][it].imag()=out1[it][1]/lt;
                } 
          }

       for(it=0;it<lt;it++)
          {
             for(ix=0;ix<nr;ix++)
               {
                  in2[ix][0]=uf[ix][it].real();
                  in2[ix][1]=uf[ix][it].imag();
               }
             
             fftwf_execute(p2);

             for(ix=0;ix<nr;ix++)
               {
                  ufk[ix][it].real()=out2[ix][0]/nr;
                  ufk[ix][it].imag()=out2[ix][1]/nr;
               }
          }

       for(ix=0;ix<nr;ix++)
          for(it=0;it<lt;it++)
            {
              amp_orig[ix][it]=sqrt(pow(ufk[ix][it].real(),2)+pow(ufk[ix][it].imag(),2));
              swq4.write((char*)&amp_orig[ix][it],sizeof(amp_orig[ix][it])); 
            } 

       for(it=0;it<lt/2+1;it++)
          {
            if((it+1)%1000==0)  
              cout<<it+1<<" frequency slices done... "<<endl;
              
              if(it<itmin)
                {
                   for(ix=0;ix<nr;ix++)
                     {
                        ufk1[ix][it].real()=0.0;
                        ufk1[ix][it].imag()=0.0;

                        ufk2[ix][it].real()=0.0;
                        ufk2[ix][it].imag()=0.0;

                        gh_amp=0.0;
                        swq3.write((char*)&gh_amp,sizeof(gh_amp));
                     }   
                }       
              else if(it>itmax&&it<lt/2+1)
                {
                   for(ix=0;ix<nr;ix++)
                     {
                        ufk1[ix][it].real()=0.0;
                        ufk1[ix][it].imag()=0.0;

                        ufk2[ix][it].real()=0.0;
                        ufk2[ix][it].imag()=0.0;

                        gh_amp=0.0;
                        swq3.write((char*)&gh_amp,sizeof(gh_amp));
                     }  
                }
              else
                {
                   for(ix=0;ix<nr;ix++)
                     {
                       if(pow(omega[it]/v,2)-pow(kx[ix],2)<0)
                          {
                            gh_amp=0.0;
                            swq3.write((char*)&gh_amp,sizeof(gh_amp));
                           
                            ufk1[ix][it].real()=0.0;
                            ufk1[ix][it].imag()=0.0;

                            ufk2[ix][it].real()=0.0;
                            ufk2[ix][it].imag()=0.0;
 
                          } 
                       else
                          {
                             ps[ix].real()=0.0;
                             ps[ix].imag()=2*sqrt(pow(omega[it]/v,2)-pow(kx[ix],2))*rz;
                             ps[ix]=exp(ps[ix]);

                             dg[ix].real()=1.0-ps[ix].real()+lamda;
                             dg[ix].imag()=ps[ix].imag();
                             
                             gh_amp=sqrt(pow(dg[ix].real()-lamda,2)+pow(dg[ix].imag(),2));
                             swq3.write((char*)&gh_amp,sizeof(gh_amp));
                          
                             ufk1[ix][it]=ufk[ix][it]/dg[ix];

                             ps2[ix].real()=0.0;
                             ps2[ix].imag()=-2*sqrt(pow(omega[it]/v,2)-pow(kx[ix],2))*rz;
                             ps2[ix]=exp(ps2[ix]);

                             dg2[ix].real()=1.0-ps2[ix].real()+lamda;
                             dg2[ix].imag()=ps2[ix].imag();
                             ufk2[ix][it]=ufk[ix][it]/dg2[ix];

                          }
                     }
                }
            }  
        swq3.close();

        ifstream swq5;
        swq5.open("gh_amp.dat",ios::binary);
        ofstream swq6;
        swq6.open("gh_amp_all.dat",ios::binary);
        ofstream swq7;
        swq7.open("filter.dat",ios::binary);

        for(it=0;it<lt/2+1;it++)
           for(ix=0;ix<nr;ix++)
              swq5.read((char*)&degh_amp_tmp[it][ix],sizeof(degh_amp_tmp[it][ix]));        
        swq5.close();

        for(ix=0;ix<nr;ix++)
           for(it=0;it<lt/2+1;it++)
               degh_amp[ix][it]=degh_amp_tmp[it][ix];

        for(ix=0;ix<nr;ix++)
           for(it=lt/2+1;it<lt;it++)
              degh_amp[ix][it]=degh_amp[ix][lt-it];     

        for(ix=0;ix<nr;ix++)
           for(it=0;it<lt;it++)
             swq6.write((char*)&degh_amp[ix][it],sizeof(degh_amp[ix][it]));
        swq6.close();

        for(ix=0;ix<nr;ix++)
          { 
            for(it=0;it<lt/2+1;it++)
               {
                  filt[ix][it]=sqrt(pow(ufk1[ix][it].real(),2)+pow(ufk1[ix][it].imag(),2));
                  swq7.write((char*)&filt[ix][it],sizeof(filt[ix][it]));          
               }
          }
        swq7.close();


        for(it=0;it<lt/2+1;it++)
           {
              for(ix=0;ix<nr;ix++)
                {
                   in3[ix][0]=ufk1[ix][it].real();
                   in3[ix][1]=ufk1[ix][it].imag();
                }

              fftwf_execute(p3);

              for(ix=0;ix<nr;ix++)
                {
                   uf1[ix][it].real()=out3[ix][0];
                   uf1[ix][it].imag()=out3[ix][1];
                } 
           }

        for(ix=0;ix<nr;ix++)
          for(it=lt/2+1;it<lt;it++)
            {
               uf1[ix][it].real()=uf1[ix][lt-it].real();
               uf1[ix][it].imag()=-uf1[ix][lt-it].imag();
            }      
        
        for(ix=0;ix<nr;ix++)
           {
              for(it=0;it<lt;it++)
                {
                   in4[it][0]=uf1[ix][it].real();
                   in4[it][1]=uf1[ix][it].imag();
                }
              
              fftwf_execute(p4);   

              for(it=0;it<lt;it++)
                 u1[ix][it]=out4[it][0];    

           }  
        
        for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++)
            swq2.write((char*)&u1[ir][it],sizeof(u1[ir][it]));

        for(it=0;it<lt/2+1;it++)
           {
              for(ix=0;ix<nr;ix++)
                {
                   in3[ix][0]=ufk2[ix][it].real();
                   in3[ix][1]=ufk2[ix][it].imag();
                }

              fftwf_execute(p3);

              for(ix=0;ix<nr;ix++)
                {
                   uf1[ix][it].real()=out3[ix][0];
                   uf1[ix][it].imag()=out3[ix][1];
                }
           }

        for(ix=0;ix<nr;ix++)
          for(it=lt/2+1;it<lt;it++)
            {
               uf1[ix][it].real()=uf1[ix][lt-it].real();
               uf1[ix][it].imag()=-uf1[ix][lt-it].imag();
            }

        for(ix=0;ix<nr;ix++)
           {
              for(it=0;it<lt;it++)
                {
                   in4[it][0]=uf1[ix][it].real();
                   in4[it][1]=uf1[ix][it].imag();
                }

              fftwf_execute(p4);

              for(it=0;it<lt;it++)
                 u1[ix][it]=out4[it][0];

           }

        for(ir=0;ir<nr;ir++)
          for(it=0;it<lt;it++)
            swq22.write((char*)&u1[ir][it],sizeof(u1[ir][it]));
 
        cout<<is+1<<" shot deghosting done!"<<endl;       
    }

   cout<<"ALL DONE!"<<endl;

   return 0;


}






































