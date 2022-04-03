#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"

int main()
{
/*
  int ixx,iyy,izz;
  float ***a;
  a=alloc3float(1,2,3);
  for(ixx=0;ixx<1;ixx++)
    for(iyy=0;iyy<2;iyy++)
      for(izz=0;izz<3;izz++)
          a[ixx][iyy][izz]=ixx+iyy+izz;
  for(ixx=0;ixx<1;ixx++)
    for(iyy=0;iyy<2;iyy++)
      for(izz=0;izz<3;izz++)
         cout<<a[ixx][iyy][izz]<<endl;
  return 0;
*/

  char fn1[256],fn2[256],fn3[256],fn4[256];
  int ns_all,ns,nr,nsr,lt,ltt,d;

  ifstream swq;
  swq.open("3d_fftw.par");
  if(!swq)
    { 
      cout<<"Cannot open 3d_fftw.par"<<endl;
      abort();
    }
  swq>>fn1>>fn2>>fn3>>fn4>>ns_all>>ns>>nr>>lt>>d;
  swq.close();

  cout<<fn1<<endl;
  cout<<fn2<<endl;
  cout<<fn3<<endl;
  cout<<fn4<<endl;
  cout<<"ns_all===="<<ns_all<<endl;
  cout<<"ns per line===="<<ns<<endl;
  cout<<"nr===="<<nr<<endl;
  cout<<"lt===="<<lt<<endl;
  cout<<"distance===="<<d<<endl;


  ltt=2*lt;
  int ix,iy,it,is,tmp;

  float **u2t;
  u2t=alloc2float(ltt,nr);

  complex <float> **u2f;
  u2f=alloc2complex(ltt,nr);
  complex <float> **u2fky;
  u2fky=alloc2complex(ltt,nr); 


  complex <float> **u2fky_tmp;
  u2fky_tmp=alloc2complex(ltt,nr);
  float **u2t_tmp;
  u2t_tmp=alloc2float(ltt,nr);



  complex <float> ***u3fft;
  u3fft=alloc3complex(ltt,nr,ns);
  complex <float> ***u3fft_tmp1;
  u3fft_tmp1=alloc3complex(ltt,nr,ns);  
  complex <float> ***u3fft_tmp2;
  u3fft_tmp2=alloc3complex(ltt,nr,ns);  
  complex <float> ***u3fft_tmp3;
  u3fft_tmp3=alloc3complex(ltt,nr,ns);  

  float ***u3t;
  u3t=alloc3float(ltt,nr,ns);


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
       cout<<"cannot create "<<fn2<<endl;
       abort();
    } 

  ofstream swq3;
  swq3.open(fn3,ios::binary);
  if(!swq3)
    {
       cout<<"cannot create "<<fn3<<endl;
       abort();
    }

  ofstream swq4;
  swq4.open(fn4,ios::binary);
  if(!swq4)
    {
       cout<<"cannot create "<<fn4<<endl;
       abort();
    }

    fftwf_complex *in1,*out1;
    fftwf_plan p1;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
    p2=fftwf_plan_dft_1d(nr,in2,out2,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in3,*out3;
    fftwf_plan p3;
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    p3=fftwf_plan_dft_1d(ns,in3,out3,FFTW_FORWARD,FFTW_MEASURE);

    fftwf_complex *in4,*out4;
    fftwf_plan p4;
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    p4=fftwf_plan_dft_1d(ltt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

    fftwf_complex *in5,*out5;
    fftwf_plan p5;
    in5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
    out5=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nr);
    p5=fftwf_plan_dft_1d(nr,in5,out5,FFTW_BACKWARD,FFTW_MEASURE);

    fftwf_complex *in6,*out6;
    fftwf_plan p6;
    in6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out6=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    p6=fftwf_plan_dft_1d(ns,in6,out6,FFTW_BACKWARD,FFTW_MEASURE);
 

//   for(is=0;is<ns_all;is++) 
   for(is=0;is<1;is++) 
     { 
       for(iy=0;iy<nr;iy++)
         for(it=0;it<ltt;it++)
            u2t[iy][it]=0.0;

       for(iy=0;iy<nr;iy++)
         for(it=0;it<lt;it++)
           swq1.read((char*)&u2t[iy][it],sizeof(u2t[iy][it]));        
/*
       ofstream swq88;
       swq88.open("check.dat",ios::binary);
       for(iy=0;iy<nr;iy++)
         for(it=0;it<lt;it++)
           swq88.write((char*)&u2t[iy][it],sizeof(u2t[iy][it]));
       return 0;
*/

       for(ix=0;ix<ns;ix++)
        {
          for(iy=0;iy<nr;iy++)
           {
            for(it=0;it<ltt;it++)
            {
              u3t[ix][iy][it]=0.0;
              u3fft[ix][iy][it].real()=0.0;
              u3fft[ix][iy][it].imag()=0.0;
            }
          }
        } 


        if((in1==NULL)||(out1==NULL))
            cout<<"memory insufficient"<<endl;
         else
           {
             for(iy=0;iy<nr;iy++)
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
               for(iy=0;iy<nr;iy++)
                 {
                    in2[iy][0]=u2f[iy][it].real();
                    in2[iy][1]=u2f[iy][it].imag();
                 }
               
               fftwf_execute(p2);

               for(iy=0;iy<nr;iy++)
                 {
                   u2fky[iy][it].real()=out2[iy][0]/nr;
                   u2fky[iy][it].imag()=out2[iy][1]/nr;
                 }
              }

        cout<<is+1<<" shot 2d fft done..."<<endl; 

/*
        for(it=0;it<ltt;it++)
              {
               for(iy=0;iy<nr;iy++)
                 {
                    in5[iy][0]=u2fky[iy][it].real();
                    in5[iy][1]=u2fky[iy][it].imag();
                 }

               fftwf_execute(p5);

               for(iy=0;iy<nr;iy++)
                 {
                   u2fky_tmp[iy][it].real()=out5[iy][0];
                   u2fky_tmp[iy][it].imag()=out5[iy][1];
                 }
              }

          for(iy=0;iy<nr;iy++)
              {
               for(it=0;it<ltt;it++)
                 {
                    in4[it][0]=u2fky_tmp[iy][it].real();
                    in4[it][1]=u2fky_tmp[iy][it].imag();
                 }

               fftwf_execute(p4);

               for(it=0;it<ltt;it++)
                   u2t_tmp[iy][it]=out4[it][0];

              }
        ofstream swq88;
        swq88.open("check1.dat",ios::binary);
        for(iy=0;iy<nr;iy++)
          for(it=0;it<ltt;it++)
           swq88.write((char*)&u2fky[iy][it].real(),sizeof(u2fky[iy][it].real()));        
        return 0;
*/


          for(iy=0;iy<nr;iy++)
           { 
             for(it=0;it<ltt;it++)
              { 
//                u3fft[d][iy][it]=u2fky[iy][it];
                u3fft[d][iy][it].real()=u2fky[iy][it].real();
                u3fft[d][iy][it].imag()=u2fky[iy][it].imag();
              }
           }   
/*
        ofstream swq888;
        swq888.open("check.dat",ios::binary);
        for(ix=0;ix<ns;ix++)
         {
          for(iy=0;iy<nr;iy++)
           {
            for(it=0;it<ltt;it++)
              swq888.write((char*)&u3fft[ix][iy][it].real(),sizeof(u3fft[ix][iy][it].real()));        
           }
         }
        cout<<"222222"<<endl;
        return 0;
*/
          
          for(it=0;it<ltt;it++)
           {
             for(iy=0;iy<nr;iy++)
               {
                  if((in1==NULL)||(out1==NULL))
                    cout<<"memory insufficient"<<endl;      

                  else
                   {
                     for(ix=0;ix<ns;ix++)
                        {
                           in3[ix][0]=u3fft[ix][iy][it].real();
                           in3[ix][1]=u3fft[ix][iy][it].imag();
                        }
                     fftwf_execute(p3);

                     for(ix=0;ix<ns;ix++)
                       {   
                         u3fft_tmp1[ix][iy][it].real()=out3[ix][0]/ns;
                         u3fft_tmp1[ix][iy][it].imag()=out3[ix][1]/ns;
                       }   
                   }
                     
               }
            if(it%100==0)         
             cout<<is+1<<" shot,  "<<it+1<<"  frequency slice 3dfft done!"<<endl;  
           }

        for(ix=0;ix<ns;ix++)
          {
            for(iy=0;iy<nr;iy++)
             {
              for(it=0;it<ltt;it++)
               {
                 swq2.write((char*)&u3fft_tmp1[ix][iy][it].real(),sizeof(u3fft_tmp1[ix][iy][it].real()));        
                 swq3.write((char*)&u3fft_tmp1[ix][iy][it].imag(),sizeof(u3fft_tmp1[ix][iy][it].imag()));        
               }
            }
          }

       cout<<is+1<<" shot 3d_fft done!"<<endl;

       cout<<is+1<<" shot 3d_ifft begins..."<<endl;
      
       cout<<is+1<<" shot 3d_ifft kx domain to x domain begins..."<<endl;

       for(it=0;it<ltt;it++)
           {
             for(iy=0;iy<nr;iy++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(ix=0;ix<ns;ix++)
                        {
                           in6[ix][0]=u3fft_tmp1[ix][iy][it].real();
                           in6[ix][1]=u3fft_tmp1[ix][iy][it].imag();
                        }
                     fftwf_execute(p6);

                     for(ix=0;ix<ns;ix++)
                       {
                         u3fft_tmp2[ix][iy][it].real()=out6[ix][0];
                         u3fft_tmp2[ix][iy][it].imag()=out6[ix][1];
                       }
                   }

               }
    
             if(it%100==0)
                cout<<is+1<<" shot, "<<it<<" frequency slices, kx to x  done..."<<endl;

           }

        cout<<is+1<<" shot 3d_ifft ky domain to y domain begins..."<<endl;
 
        for(it=0;it<ltt;it++)
           {
             for(ix=0;ix<ns;ix++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(iy=0;iy<nr;iy++)
                        {
                           in5[iy][0]=u3fft_tmp2[ix][iy][it].real();
                           in5[iy][1]=u3fft_tmp2[ix][iy][it].imag();
                        }

                     fftwf_execute(p5);

                     for(iy=0;iy<nr;iy++)
                       {
                         u3fft_tmp3[ix][iy][it].real()=out5[iy][0];
                         u3fft_tmp3[ix][iy][it].imag()=out5[iy][1];
                       }
                   }
               }

             if(it%100==0)
                cout<<is+1<<" shot, "<<it<<" frequency slices, ky to y  done..."<<endl;
           }

        for(ix=0;ix<ns;ix++)
           {
             for(iy=0;iy<nr;iy++)
               {
                  if((in2==NULL)||(out2==NULL))
                    cout<<"memory insufficient"<<endl;

                  else
                   {
                     for(it=0;it<ltt;it++)
                        {
                           in4[it][0]=u3fft_tmp3[ix][iy][it].real();
                           in4[it][1]=u3fft_tmp3[ix][iy][it].imag();
                        }

                     fftwf_execute(p4);
 
                     for(it=0;it<lt;it++)
                       u3t[ix][iy][it]=out4[it][0];
 
                   } 
               }
          }

        for(ix=0;ix<ns;ix++)
          {
             for(iy=0;iy<nr;iy++)
              {
                for(it=0;it<lt;it++)
                  swq4.write((char*)&u3t[ix][iy][it],sizeof(u3t[ix][iy][it]));
              }
          }
     }

        cout<<"all done!"<<endl;

     swq1.close();
     swq2.close();
     swq3.close();
     swq4.close();

     return 0;

}






