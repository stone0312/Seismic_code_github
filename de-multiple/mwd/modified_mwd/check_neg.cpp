#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

int main()
 { 

   int is,it,nr,lt,ltt;
   nr=500;
   lt=5000;
   ltt=2*lt; 

   complex<float> *f_tmp;
    f_tmp=alloc1complex(ltt);

    float *m;
    m=alloc1float(ltt);

    fftwf_complex *in2,*out2;
    fftwf_plan p2;
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);

    p2=fftwf_plan_dft_1d(ltt,in2,out2,FFTW_BACKWARD,FFTW_MEASURE);

    ofstream swq2;
    swq2.open("/data1/swq/srme_data/check_1st_shot_water_layer_related_multiple_v11_boundary1000.dat",ios::binary);
    if(!swq2)
       {
          cout<<"cannot open "<<"/data1/swq/srme_data/check_1st_shot_water_layer_related_multiple_v11_boundary1000.dat"<<endl;
          abort();
       }
 
       ofstream swq666;
       swq666.open("/data1/swq/srme_data/check_fft_predicted_model_real_v11_all.dat",ios::binary);
       if(!swq666)
         {
          cout<<"cannot open "<<"/data1/swq/srme_data/check_fft_predicted_model_real_v11_all.dat"<<endl;
          abort();
         }
    
       ofstream swq777;
       swq777.open("/data1/swq/srme_data/check_fft_predicted_model_imag_v11_all.dat",ios::binary);
       if(!swq777)
         {
            cout<<"cannot open "<<"/data1/swq/srme_data/check_fft_predicted_model_imag_v11_all.dat"<<endl;
            abort();
         }
       
       ifstream swq66;
       swq66.open("/data1/swq/srme_data/fft_predicted_model_real_v11_positive.dat",ios::binary);
       if(!swq66)
         {
           cout<<"cannot open "<<"/data1/swq/srme_data/fft_predicted_model_real_v11_positive.dat"<<endl;
           abort();
         }

       ifstream swq77;
       swq77.open("/data1/swq/srme_data/fft_predicted_model_imag_v11_positive.dat",ios::binary);
       if(!swq77)
         {
           cout<<"cannot open "<<"/data1/swq/srme_data/fft_predicted_model_imag_v11_positive.dat"<<endl;
           abort();
         }

//start to compute the predicticted model of the negative frequency
        cout<<"Negative frequcies prediction Done!"<<endl;
        cout<<"IFFT Starts..."<<endl;

        for(is=0;is<nr;is++)
          {
           for(it=0;it<lt+1;it++)
             {
               swq66.seekg((it*nr+is)*4,ios::beg);
               swq77.seekg((it*nr+is)*4,ios::beg);
               swq66.read((char*)&f_tmp[it].real(),sizeof(f_tmp[it].real()));
               swq77.read((char*)&f_tmp[it].imag(),sizeof(f_tmp[it].imag()));
             }
           for(it=lt+1;it<ltt;it++)
               {
                 f_tmp[it].real()=f_tmp[ltt-it].real();
                 f_tmp[it].imag()=-f_tmp[ltt-it].imag();
               }
           for(it=0;it<ltt;it++)
               {
                 swq666.write((char*)&f_tmp[it].real(),sizeof(f_tmp[it].real()));
                 swq777.write((char*)&f_tmp[it].imag(),sizeof(f_tmp[it].imag()));
               }

           if((in2==NULL)||(out2==NULL))
              cout<<"memory insufficient"<<endl;
           else
              {
                for(it=0;it<ltt;it++)
                  {
                     in2[it][0]=f_tmp[it].real();
                     in2[it][1]=f_tmp[it].imag();
                  }
              }
           fftwf_execute(p2);

              for(it=0;it<ltt;it++)
                 m[it]=out2[it][0];

              for(it=0;it<ltt;it++)
                swq2.write((char*)&m[it],sizeof(m[it]));

          }

        swq2.close();
        swq66.close();
        swq77.close();
        swq666.close();
        swq777.close();

        return 0;
}
