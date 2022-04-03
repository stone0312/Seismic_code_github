#include "iostream.h"
#include "fstream.h"

using namespace std;

int main()
{
   int i,j,k,is,it;
   int nshot,nt,lt,ltt;
   char input[256],inputcp[256],output[256];
   
   ifstream swq;
   swq.open("convolution_srm.par");
   if(!swq)
     {
        cout<<"cannot open convolution_srm.par"<<endl;
        abort();
     }
   swq>>input>>inputcp>>output>>nshot>>nt>>lt;
   swq.close();

   cout<<"fna of original shot gather is===="<<input<<endl;
   cout<<"fna of output multiple model is===="<<output<<endl;
   cout<<"No. of shots is===="<<nshot<<endl;
   cout<<"No. of traces per shot is===="<<nt<<endl;
   cout<<"No. of temperal samples is===="<<lt<<endl;

   ltt=2*lt;

   float **us=new float *[nt];
   for(i=0;i<nt;i++)
     us[i]=new float [ltt];

   float **ur=new float *[nshot];
   for(i=0;i<nshot;i++)
     ur[i]=new float [ltt];

  ifstream swq1;
  swq1.open(input,ios::binary);
  if(!swq1)
    {
       cout<<"cannot open "<<input<<endl;
       abort();
    } 
  ifstream swq2;
  swq2.open(inputcp,ios::binary);
  if(!swq2)
    {
       cout<<"cannot open "<<inputcp<<endl;
       abort();
    } 
  
  ofstream swq3;
  swq3.open("check_read_common_shot.dat",ios::binary);
  if(!swq3)
    {
       cout<<"cannot open "<<output<<endl;
       abort();
    }

  ofstream swq4;
  swq4.open("check_read_common_receiver.dat",ios::binary);
  if(!swq4)
    {
       cout<<"cannot open "<<output<<endl;
       abort();
    } 
 
  for(is=0;is<nshot;is++)
     {
        for(i=0;i<nt;i++)
          for(j=0;j<ltt;j++)
             us[i][j]=0.0; 
             
        for(i=0;i<nt;i++)
          for(j=0;j<lt;j++)
             swq1.read((char*)&us[i][j],sizeof(us[i][j]));//read common shot gather 
        for(i=0;i<nt;i++)
          for(j=0;j<ltt;j++) 
              swq3.write((char*)&us[i][j],sizeof(us[i][j]));
         cout<<is+1<<" shot read and write done"<<endl;
      }
      cout<<"shot gather done"<<endl;

   for(it=0;it<nt;it++)
     {
       for(i=1;i<nshot;i++)
          for(j=0;j<ltt;j++)
             ur[i][j]=0.0; 

       swq2.seekg(it*lt*4,ios::beg);
       for(j=0;j<lt;j++)
          swq2.read((char*)&ur[0][j],sizeof(ur[0][j])); 

        for(i=1;i<nshot;i++)//read common receiver gather
           {
              swq2.seekg((nt-1)*lt*4,ios::cur);
              for(j=0;j<lt;j++)
                swq2.read((char*)&ur[i][j],sizeof(ur[i][j]));
           }
         
        for(i=0;i<nshot;i++)
           for(j=0;j<ltt;j++)
               swq4.write((char*)&ur[i][j],sizeof(ur[i][j]));
         cout<<it+1<<" receiver read and write done"<<endl;      
      }

     cout<<"ALL DONE!"<<endl;
     swq1.close();
     swq2.close();
     swq3.close();
     swq4.close();
 return 0;

}











