#include <iostream>
#include "math.h"
#include <fstream>
#include "stdlib.h"
#include <complex>

using namespace std;

#include "alloc.c"
#include "fftw3.h"
#include "mpi.h"
#define pai 3.14159265

int main(int argc,char **argv)
{

  char fn1[256],fn2[256],fn3[256],fn4[256];
  int ns_all,ns,nr,nsb,nrb,lt,ltt,bou_x,bou_y,d;
  float dx,dy,dt,fmin,fmax,v,theta,dep,dep_s,lamda;
  float kz;
  complex<float> ps;
  complex<float> ps1;
  complex<float> ps2;
  complex<float> ps3;
  complex<float> ps4;
 
  int np,myid;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);


  ifstream swq;
  swq.open("3d_mwgd.par");
  if(!swq)
    { 
      cout<<"Cannot open 3d_mwgd.par"<<endl;
      abort();
    }
  swq>>fn1>>fn2>>fn3>>ns_all>>ns>>nr>>dx>>dy>>lt>>dt>>fmin>>fmax>>v>>bou_x>>bou_y>>theta>>d>>dep>>dep_s>>lamda;
  swq.close();

  cout<<fn1<<endl;
  cout<<fn2<<endl;
  cout<<fn3<<endl;
  cout<<"ns_all===="<<ns_all<<endl;
  cout<<"ns per line===="<<ns<<endl;
  cout<<"nr===="<<nr<<endl;
  cout<<"lt===="<<lt<<endl;
  cout<<"dt===="<<dt<<endl;
  cout<<"dx===="<<dx<<endl;
  cout<<"dy===="<<dy<<endl;
  cout<<"fmin===="<<fmin<<endl;
  cout<<"fmax===="<<fmax<<endl;
  cout<<"v===="<<v<<endl;
  cout<<"bou_x===="<<bou_x<<endl;
  cout<<"bou_y===="<<bou_y<<endl;
  cout<<"theta===="<<theta<<endl;
  cout<<"distance===="<<d<<endl;
  cout<<"depth of water bottom is===="<<dep<<endl;
  cout<<"depth of streamers is===="<<dep_s<<endl;
  cout<<"regulization parameter is===="<<lamda<<endl;


//  ltt=2*lt;

  ltt=int(1.3*lt);
  if(ltt%2!=0)
    ltt+=1;     

  nsb=ns+2*bou_x;
  nrb=nr+2*bou_y;
  int ix,iy,it,is,tmp;
  int itmin,itmax;

  
//  cout<<nsb*nrb*ltt/1024/1024/4<<endl;
//  cout<<nsb<<" , "<<nrb<<" , "<<ltt<<endl;
//  return 0;

  float *omega;
  omega=alloc1float(ltt);
  float *kx;
  kx=alloc1float(nsb);
  float *ky;
  ky=alloc1float(nrb);

  for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
  for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

  itmin=int(fmin*dt*ltt/1000);
  itmax=int(fmax*dt*ltt/1000)+1;

  cout<<"Totally "<<itmax-itmin+1<<" Frequency Slices Needed to be Calculated..."<<endl;

  for(is=0;is<nsb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dx*nsb);
   for(is=nsb/2+1;is<nsb;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nsb/2)/float(dx*nsb);

   for(is=0;is<nrb/2+1;is++)
      ky[is]=2*pai*float(is)/float(dy*nrb);
   for(is=nrb/2+1;is<nrb;is++)
      ky[is]=-2*pai*1.0/float(2*dy)+2*pai*float(is-nrb/2)/float(dy*nrb);

  float **u2t;
  u2t=alloc2float(ltt,nrb);

  complex <float> **u2f;
  u2f=alloc2complex(ltt,nrb);
  complex <float> **u2fky;
  u2fky=alloc2complex(ltt,nrb); 

  complex <float> ***u3fft;
  u3fft=alloc3complex(ltt,nrb,nsb);

  complex <float> ***upkf;
  upkf=alloc3complex(ltt,nrb,nsb);


  float mem;
  mem=0.0;
  mem+=(ltt+nsb+nrb+ltt*nrb*5+ltt*nrb*nsb*2*2);
  mem=mem*4/1024.0/1024.0;
  cout<<"Memory Needed to be Allocated is===="<<mem<<"MB"<<endl;

//  return 0;

  fstream swq1;
  swq1.open(fn1,ios::binary|ios::in);
  if(!swq1)
    {
       cout<<"cannot open "<<fn1<<endl;
       abort();
    } 

  fstream swq2;
  swq2.open(fn2,ios::binary|ios::out);
  if(!swq2)
    {
       cout<<"cannot create "<<fn2<<endl;
       abort();
    }

  fstream swq3;
  swq3.open(fn3,ios::binary|ios::out);
  if(!swq3)
    {
       cout<<"cannot create "<<fn3<<endl;
       abort();
    }

    fftwf_complex *in1,*out1;
    fftwf_plan p1;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    p2=fftwf_plan_dft_1d(nrb,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in3,*out3;
    fftwf_plan p3;
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nsb);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nsb);
    p3=fftwf_plan_dft_1d(nsb,in3,out3,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in4,*out4;
    fftwf_plan p4;
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p4=fftwf_plan_dft_1d(ltt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

    fftwf_complex *in5,*out5;
    fftwf_plan p5;
    in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nrb);
    p5=fftwf_plan_dft_1d(nrb,in5,out5,FFTW_BACKWARD,FFTW_MEASURE);

    fftwf_complex *in6,*out6;
    fftwf_plan p6;
    in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nsb);
    out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nsb);
    p6=fftwf_plan_dft_1d(nsb,in6,out6,FFTW_BACKWARD,FFTW_MEASURE);

    for(is=myid;is<ns_all;is+=np)
//   for(is=0;is<1;is++) 
     { 
       for(iy=0;iy<nrb;iy++)
         for(it=0;it<ltt;it++)
            u2t[iy][it]=0.0;
 
       swq1.seekg(0,ios::beg);
       for(ix=0;ix<is;ix++)
          swq1.seekg(nr*lt*4,ios::cur);

       for(iy=bou_y;iy<nr+bou_y;iy++)
         for(it=0;it<lt;it++)
           swq1.read((char*)&u2t[iy][it],sizeof(u2t[iy][it]));        

       for(ix=0;ix<nsb;ix++)
        {
          for(iy=0;iy<nrb;iy++)
           {
            for(it=0;it<ltt;it++)
            {
              u3fft[ix][iy][it].real()=0.0;
              u3fft[ix][iy][it].imag()=0.0;
            }
          }
        } 

        if((in1==NULL)||(out1==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
             for(iy=0;iy<nrb;iy++)
              {
               for(it=0;it<ltt;it++)
                 {
                    in1[it][0]=u2t[iy][it];
                    in1[it][1]=0.0;
                 }
               
               fftwf_execute(p1);

               for(it=0;it<ltt;it++)
                 {
                   u2f[iy][it].real()=out1[it][0]/ltt;
                   u2f[iy][it].imag()=out1[it][1]/ltt;
                 }

              }
           }

         for(it=0;it<ltt;it++)
              {
               for(iy=0;iy<nrb;iy++)
                 {
                    in2[iy][0]=u2f[iy][it].real();
                    in2[iy][1]=u2f[iy][it].imag();
                 }
               
               fftwf_execute(p2);

               for(iy=0;iy<nrb;iy++)
                 {
                   u2fky[iy][it].real()=out2[iy][0]/nrb;
                   u2fky[iy][it].imag()=out2[iy][1]/nrb;
                 }
              }

        cout<<is+1<<" shot 2d fft done..."<<endl; 

          for(iy=0;iy<nrb;iy++)
           { 
             for(it=0;it<ltt;it++)
              { 
                u3fft[d+bou_x][iy][it].real()=u2fky[iy][it].real();
                u3fft[d+bou_x][iy][it].imag()=u2fky[iy][it].imag();
              }
           }   
          
          for(it=0;it<ltt;it++)
           {
             for(iy=0;iy<nrb;iy++)
               {
                  if((in1==NULL)||(out1==NULL))
                    cout<<"memory insufficient"<<endl;      

                  else
                   {
                     for(ix=0;ix<nsb;ix++)
                        {
                           in3[ix][0]=u3fft[ix][iy][it].real();
                           in3[ix][1]=u3fft[ix][iy][it].imag();
                        }
                     fftwf_execute(p3);

                     for(ix=0;ix<nsb;ix++)
                       {   
                         u3fft[ix][iy][it].real()=out3[ix][0]/nsb;
                         u3fft[ix][iy][it].imag()=out3[ix][1]/nsb;
                       }   
                   }
               }
            if(it%100==0)         
             cout<<is+1<<" shot,  "<<it+1<<"  frequency slice 3dfft done!"<<endl;  
           }


       cout<<is+1<<" shot 3d_fft done!"<<endl;
  
       cout<<is+1<<" shot 3d_mwgd begins..."<<endl;

       for(ix=0;ix<nsb;ix++)
          {
            for(iy=0;iy<nrb;iy++)
              {
                  for(it=0;it<ltt/2+1;it++)
                     {
                       if(it<itmin)
                          {
                            upkf[ix][iy][it].real()=0.0;//upward wavefield
                            upkf[ix][iy][it].imag()=0.0;
                            u3fft[ix][iy][it].real()=0.0;//downward wavefield
                            u3fft[ix][iy][it].imag()=0.0;
                          }
                       else if(it>itmax&&it<ltt/2+1)
                          {
                            upkf[ix][iy][it].real()=0.0;
                            upkf[ix][iy][it].imag()=0.0;
                            u3fft[ix][iy][it].real()=0.0;
                            u3fft[ix][iy][it].imag()=0.0;
                          }
                       else
                          {
                             if(pow(omega[it]/v,2)-pow(kx[ix],2)-pow(ky[iy],2)<0.0) 
                               {
                                  upkf[ix][iy][it].real()=0.0;
                                  upkf[ix][iy][it].imag()=0.0;
                                  u3fft[ix][iy][it].real()=0.0;
                                  u3fft[ix][iy][it].imag()=0.0;
                               }      
                             else if(pow(omega[it]/v,2)-pow(kx[ix],2)-pow(ky[iy],2)>pow(cos(2*pai*theta/360.0)*omega[it]/v,2))
                               {
                                  kz=sqrt(pow(omega[it]/v,2)-pow(kx[ix],2)-pow(ky[iy],2));

/*
                                  ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                                  ps1.imag()=-2*dep_s*kz;
                                  ps1=exp(ps1);
                                  ps2.real()=1.0-ps1.real();
                                  ps2.imag()=-ps1.imag();  

                                  ps3.real()=0.0;               //phase-shift extrapolating operator to achieve downward field
                                  ps3.imag()=2*dep_s*kz;
                                  ps3=exp(ps3);
                                  ps4.real()=1.0-ps3.real();
                                  ps4.imag()=-ps3.imag();  
                                                               
            
                                  upkf[ix][iy][it].real()=(u3fft[ix][iy][it].real()*ps2.real()+u3fft[ix][iy][it].imag()*ps2.imag())/(ps2.real()*ps2.real()-ps2.imag()*ps2.imag()+lamda);
                                  upkf[ix][iy][it].imag()=(u3fft[ix][iy][it].imag()*ps2.real()-u3fft[ix][iy][it].real()*ps2.imag())/(ps2.real()*ps2.real()-ps2.imag()*ps2.imag()+lamda);

                                  u3fft[ix][iy][it].real()=(u3fft[ix][iy][it].real()*ps4.real()+u3fft[ix][iy][it].imag()*ps4.imag())/(ps4.real()*ps4.real()-ps4.imag()*ps4.imag()+lamda);
                                  u3fft[ix][iy][it].imag()=(u3fft[ix][iy][it].imag()*ps4.real()-u3fft[ix][iy][it].real()*ps4.imag())/(ps4.real()*ps4.real()-ps4.imag()*ps4.imag()+lamda);
*/

                                  ps1.real()=0.0;               //phase-shift extrapolating operator to achieve upward field
                                  ps1.imag()=-2*dep_s*kz;
                                  ps1=exp(ps1);
                                  ps2.real()=1.0-ps1.real()+lamda;
                                  ps2.imag()=-ps1.imag();  

                                  ps3.real()=0.0;               //phase-shift extrapolating operator to achieve downward field
                                  ps3.imag()=2*dep_s*kz;
                                  ps3=exp(ps3);
                                  ps4.real()=1.0-ps3.real()+lamda;
                                  ps4.imag()=-ps3.imag();  
 
                                  upkf[ix][iy][it]=u3fft[ix][iy][it]/ps2;   //upward field
                                  u3fft[ix][iy][it]=u3fft[ix][iy][it]/ps4;  //downward field



                                  ps.real()=0.0;               //phase-shift extrapolating operator to achieve wlrm
                                  ps.imag()=-2*(dep-dep_s)*kz;
                                  ps=exp(ps);

                                  u3fft[ix][iy][it]=u3fft[ix][iy][it]*ps;
                               }
                            else
                               {
                                  upkf[ix][iy][it].real()=0.0;
                                  upkf[ix][iy][it].imag()=0.0;  

                                  u3fft[ix][iy][it].real()=0.0;
                                  u3fft[ix][iy][it].imag()=0.0;
                               }
                          }
                     }
                 for(it=ltt/2+1;it<ltt;it++)
                    {
                      upkf[ix][iy][it].real()=upkf[ix][iy][ltt-it].real();
                      upkf[ix][iy][it].imag()=-upkf[ix][iy][ltt-it].imag();

                      u3fft[ix][iy][it].real()=u3fft[ix][iy][ltt-it].real();
                      u3fft[ix][iy][it].imag()=-u3fft[ix][iy][ltt-it].imag();
                    }

              }
          }

       cout<<is+1<<" shot 3d_mwgd  in k-f domain done!"<<endl;

       cout<<is+1<<" shot 3d_ifft of upward wavefilds  begins..."<<endl;
      
       cout<<is+1<<" shot 3d_ifft kx domain to x domain begins..."<<endl;

       for(it=0;it<ltt;it++)
           {
             for(iy=0;iy<nrb;iy++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(ix=0;ix<nsb;ix++)
                        {
                           in6[ix][0]=upkf[ix][iy][it].real();
                           in6[ix][1]=upkf[ix][iy][it].imag();
                        }
                     fftwf_execute(p6);

                     for(ix=0;ix<nsb;ix++)
                       {
                         upkf[ix][iy][it].real()=out6[ix][0];
                         upkf[ix][iy][it].imag()=out6[ix][1];
                       }
                   }

               }
    
             if(it%100==0)
                cout<<is+1<<" shot, "<<it<<" frequency slices, kx to x  done..."<<endl;
           }

        cout<<is+1<<" shot 3d_ifft ky domain to y domain begins..."<<endl;
 
        for(it=0;it<ltt;it++)
           {
             for(ix=0;ix<nsb;ix++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(iy=0;iy<nrb;iy++)
                        {
                           in5[iy][0]=upkf[ix][iy][it].real();
                           in5[iy][1]=upkf[ix][iy][it].imag();
                        }

                     fftwf_execute(p5);

                     for(iy=0;iy<nrb;iy++)
                       {
                         upkf[ix][iy][it].real()=out5[iy][0];
                         upkf[ix][iy][it].imag()=out5[iy][1];
                       }
                   }
               }

             if(it%100==0)
                cout<<is+1<<" shot, "<<it<<" frequency slices, ky to y  done..."<<endl;
           }

        for(ix=0;ix<nsb;ix++)
           {
             for(iy=0;iy<nrb;iy++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(it=0;it<ltt;it++)
                        {
                           in4[it][0]=upkf[ix][iy][it].real();
                           in4[it][1]=upkf[ix][iy][it].imag();
                        }

                     fftwf_execute(p4);
 
                     for(it=0;it<lt;it++)
                       upkf[ix][iy][it].real()=out4[it][0];
 
                   } 
               }
          }

         swq3.seekg(0,ios::beg);
         for(ix=0;ix<is;ix++)
           swq3.seekg(nr*lt*4,ios::cur);

             for(iy=0;iy<nr;iy++)
              {
                for(it=0;it<lt;it++)
                  swq3.write((char*)&upkf[d+bou_x][iy+bou_y][it].real(),sizeof(upkf[d+bou_x][iy+bou_y][it].real()));
              }

              cout<<is+1<<" shot 3d_deghosting done!"<<endl;

       cout<<is+1<<" shot 3d_ifft of wlrm  begins..."<<endl;
      
       cout<<is+1<<" shot wlrm 3d_ifft kx domain to x domain begins..."<<endl;

       for(it=0;it<ltt;it++)
           {
             for(iy=0;iy<nrb;iy++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(ix=0;ix<nsb;ix++)
                        {
                           in6[ix][0]=u3fft[ix][iy][it].real();
                           in6[ix][1]=u3fft[ix][iy][it].imag();
                        }
                     fftwf_execute(p6);

                     for(ix=0;ix<nsb;ix++)
                       {
                         u3fft[ix][iy][it].real()=out6[ix][0];
                         u3fft[ix][iy][it].imag()=out6[ix][1];
                       }
                   }

               }
    
             if(it%100==0)
                cout<<is+1<<" shot, "<<it<<" frequency slices, kx to x  done..."<<endl;

           }

        cout<<is+1<<" shot wlrm 3d_ifft ky domain to y domain begins..."<<endl;
 
        for(it=0;it<ltt;it++)
           {
             for(ix=0;ix<nsb;ix++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(iy=0;iy<nrb;iy++)
                        {
                           in5[iy][0]=u3fft[ix][iy][it].real();
                           in5[iy][1]=u3fft[ix][iy][it].imag();
                        }

                     fftwf_execute(p5);

                     for(iy=0;iy<nrb;iy++)
                       {
                         u3fft[ix][iy][it].real()=out5[iy][0];
                         u3fft[ix][iy][it].imag()=out5[iy][1];
                       }
                   }
               }

             if(it%100==0)
                cout<<is+1<<" shot, "<<it<<" frequency slices, ky to y  done..."<<endl;
           }

        for(ix=0;ix<nsb;ix++)
           {
             for(iy=0;iy<nrb;iy++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(it=0;it<ltt;it++)
                        {
                           in4[it][0]=u3fft[ix][iy][it].real();
                           in4[it][1]=u3fft[ix][iy][it].imag();
                        }

                     fftwf_execute(p4);
 
                     for(it=0;it<lt;it++)
                       u3fft[ix][iy][it].real()=out4[it][0];
 
                   } 
               }
          }

         swq2.seekg(0,ios::beg);
         for(ix=0;ix<is;ix++)
          swq2.seekg(nr*lt*4,ios::cur);

             for(iy=0;iy<nr;iy++)
              {
                for(it=0;it<lt;it++)
                  swq2.write((char*)&u3fft[d+bou_x][iy+bou_y][it].real(),sizeof(u3fft[d+bou_x][iy+bou_y][it].real()));
              }

              cout<<is+1<<" shot 3d_mwd done!"<<endl;

    }

        cout<<"all done!"<<endl;

     swq1.close();
     swq2.close();
     swq3.close();

     return 0;

}






