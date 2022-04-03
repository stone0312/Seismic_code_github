#include "iostream.h"
#include "fstream.h"
#include "iomanip.h"
#include "math.h"

using namespace std;

int main()
{
   int i,j;  
   int ntrace,lt;
	 char input1[256],input2[256],output[256];
	 
   ifstream swq88;
	 swq88.open("substract.par");
	 swq88>>input1>>input2>>output>>ntrace>>lt;
	 swq88.close();

	 cout<<"no. of traces to be substracted for comparision===="<<ntrace<<endl;
	 
   float **u1=new float *[ntrace];
   for(i=0;i<ntrace;i++)
         u1[i]=new float [lt];

   float **u2=new float *[ntrace];
   for(i=0;i<ntrace;i++)
         u2[i]=new float [lt];

   float **u3=new float *[ntrace];
   for(i=0;i<ntrace;i++)
         u3[i]=new float [lt];

   ifstream swq1;
   swq1.open(input1,ios::binary);
   if(!swq1)
   {
        cout<<"cannot open"<<input1<<endl;
        abort();
   }

     for(i=0;i<ntrace;i++)
           for(j=0;j<lt;j++)
                swq1.read((char*)&u1[i][j],sizeof(u1[i][j]));
   swq1.close();
 
   ifstream swq2;
   swq2.open(input2,ios::binary);
   if(!swq2)
   {
        cout<<"cannot open"<<input2<<endl;
        abort();
   }

     for(i=0;i<ntrace;i++)
           for(j=0;j<lt;j++)
                swq2.read((char*)&u2[i][j],sizeof(u2[i][j]));
   swq2.close();

   for(i=0;i<ntrace;i++)
      for(j=0;j<lt;j++)
         u3[i][j]=u1[i][j]-u2[i][j];

    ofstream swq3;
    swq3.open(output,ios::binary);
    if(!swq3)
    {
        cout<<"cannot open"<<output<<endl;
        abort();
    }
     
    for(i=0;i<ntrace;i++)
         for(j=0;j<lt;j++)
            swq3.write((char*)&u3[i][j],sizeof(u3[i][j]));
    swq3.close();

    return 0;
}
