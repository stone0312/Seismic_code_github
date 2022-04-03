#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"

int main()
{
   int ns,lt,ltt;
   ns=200;
   lt=4000;
   ltt=2*lt;

   int is,it;

   float **u;
   u=alloc2float(lt,ns);

   float **amp;
   amp=alloc2float(ltt,ns);


   float **u1;
   u1=alloc2float(ltt,ns);

   complex<float> **u1f;
   u1f=alloc2complex(ltt,ns);

   complex<float> **u1fk;
   u1fk=alloc2complex(ltt,ns);

   complex<float> **u1fkf;
   u1fkf=alloc2complex(ltt,ns);

   complex<float> **u1f1;
   u1f1=alloc2complex(ltt,ns);
 
   complex<float> **u1f2;
   u1f2=alloc2complex(ltt,ns);

   for(is=0;is<ns;is++)
      for(it=0;it<ltt;it++)
        u1[is][it]=0.0;

    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3,*in4,*out4;
    fftwf_plan p1,p2,p3,p4;

    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ns);
    in4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out4=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p1=fftwf_plan_dft_1d(ltt,in1,out1,FFTW_FORWARD,FFTW_MEASURE);
    p2=fftwf_plan_dft_1d(ns,in2,out2,FFTW_FORWARD,FFTW_MEASURE);
    p3=fftwf_plan_dft_1d(ns,in3,out3,FFTW_BACKWARD,FFTW_MEASURE);
    p4=fftwf_plan_dft_1d(ltt,in4,out4,FFTW_BACKWARD,FFTW_MEASURE);

    ifstream swq1;
    swq1.open("/data1/swq/srme/orig_shots.dat",ios::binary);

    ofstream swq2;
    swq2.open("/data1/swq/srme/test_2dfft.dat",ios::binary);
    
    for(is=0;is<ns;is++)
       {
          for(it=0;it<lt;it++)
            swq1.read((char*)&u[is][it],sizeof(u[is][it]));
       }
    swq1.close();


    for(is=0;is<ns;is++)
       for(it=0;it<lt;it++)
          u1[is][it]=u[is][it];
/*
    ofstream swq22;
    swq22.open("/data1/swq/srme/test.dat",ios::binary);
    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq22.write((char*)&u1[is][it],sizeof(u1[is][it]));      
    swq22.close();
    return 0;
*/

    for(is=0;is<ns;is++)
       {
          if((in1==NULL)||(out1==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(it=0;it<ltt;it++)
                  {
                     in1[it][0]=u1[is][it];
                     in1[it][1]=0.0;
                  }
               fftwf_execute(p1);

               for(it=0;it<ltt;it++)
                  {
                     u1f[is][it].real()=out1[it][0]/ltt;
                     u1f[is][it].imag()=out1[it][1]/ltt;
                  }
            } 
       }  

/*
     ofstream swq22;
    swq22.open("/data1/swq/srme/test.dat",ios::binary);
    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq22.write((char*)&u1f[is][it].real(),sizeof(u1f[is][it].real()));      
    swq22.close();
    return 0;
*/

 
    for(is=0;is<ns;is++)
      for(it=ltt/2+1;it<ltt;it++)
         {
             u1f[is][it].real()=0.0;
             u1f[is][it].imag()=0.0;
         }  

    for(is=0;is<ns;is++)
      for(it=ltt/2+1;it<ltt;it++)
         {
            u1f[is][it].real()=u1f[is][ltt-it].real();
            u1f[is][it].imag()=-u1f[is][ltt-it].imag();
         }


/*
     for(is=0;is<ns;is++)
       {
          if((in4==NULL)||(out4==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(it=0;it<ltt;it++)
                  {
                     in4[it][0]=u1f[is][it].real();
                     in4[it][1]=u1f[is][it].imag();
                  }
               fftwf_execute(p4);

               for(it=0;it<ltt;it++)
                  {
                     u1f2[is][it].real()=out4[it][0];
                     u1f2[is][it].imag()=out4[it][1];
                  }
            }
        }
     ofstream swq22;
    swq22.open("/data1/swq/srme/test.dat",ios::binary);
    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq22.write((char*)&u1f2[is][it].real(),sizeof(u1f2[is][it].real()));      
    swq22.close();
    return 0;
*/
   
    for(it=0;it<ltt;it++)
      {
          if((in2==NULL)||(out2==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(is=0;is<ns;is++)
                  {
                     in2[is][0]=u1f[is][it].real();
                     in2[is][1]=u1f[is][it].imag();
                  }
               fftwf_execute(p2);

               for(is=0;is<ns;is++)
                  {
                     u1fk[is][it].real()=out2[is][0]/ns;
                     u1fk[is][it].imag()=out2[is][1]/ns;
                  } 
            }
      }
/*
   cout<<"(0,1)==="<<u1fk[100][200]<<endl;
   cout<<"("<<ns-1<<","<<ltt-1<<")==="<<u1fk[ns-100][ltt-200]<<endl;
   return 0;
*/

/*
    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
          amp[is][it]=sqrt(u1fk[is][it].real()*u1fk[is][it].real()+u1fk[is][it].imag()*u1fk[is][it].imag());
*/
     for(it=ltt/2+1;it<ltt;it++)
         {
            u1fk[0][it].real()=u1fk[0][ltt-it].real();
            u1fk[0][it].imag()=-u1fk[0][ltt-it].imag();
         }

     for(it=ltt/2+1;it<ltt;it++)
       {
         for(is=1;is<ns;is++)
           {
              u1fk[is][it].real()=u1fk[ns-is][ltt-it].real();
              u1fk[is][it].imag()=-u1fk[ns-is][ltt-it].imag();
           } 
       } 


    for(it=0;it<ltt;it++)
      {
          if((in3==NULL)||(out3==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(is=0;is<ns;is++)
                  {
                     in3[is][0]=u1fk[is][it].real();
                     in3[is][1]=u1fk[is][it].imag();
                  }
               fftwf_execute(p3);
    
               for(is=0;is<ns;is++)
                  {
                     u1f1[is][it].real()=out3[is][0];
                     u1f1[is][it].imag()=out3[is][1];
                  } 
            }
      }
/*
   for(it=0;it<ltt;it++)
      cout<<it<<"  "<<u1f1[100][it].imag()<<endl;
   return 0; 


    ofstream swq222;
    swq222.open("/data1/swq/srme/test_fkif_real.dat",ios::binary);
    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq222.write((char*)&u1f1[is][it].real(),sizeof(u1f1[is][it].real()));
    swq222.close();

    ofstream swq22;
    swq22.open("/data1/swq/srme/test_fkif_imag.dat",ios::binary);
    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq22.write((char*)&u1f1[is][it].imag(),sizeof(u1f1[is][it].imag()));      
    swq22.close();
*/

     for(is=0;is<ns;is++)
       {
          if((in4==NULL)||(out4==NULL))
            cout<<"memory insufficient"<<endl;
          else
            {
               for(it=0;it<ltt;it++)
                  {
                     in4[it][0]=u1f1[is][it].real();
                     in4[it][1]=u1f1[is][it].imag();
                  }
               fftwf_execute(p4);

               for(it=0;it<ltt;it++)
                  {
                     u1f2[is][it].real()=out4[it][0];
                     u1f2[is][it].imag()=out4[it][1];
                  }
            }
       }

    for(is=0;is<ns;is++)
       for(it=0;it<ltt;it++)
         swq2.write((char*)&u1f2[is][it].real(),sizeof(u1f2[is][it].real()));
    swq2.close();

    cout<<"ALL DONE!"<<endl;
   
    return 0;


}



















































