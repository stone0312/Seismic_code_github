#include "iostream.h"
#include "fstream.h"
#include "iomanip.h"
#include "vector.h"
using namespace std;
int main()
{
  int i,j,ntrace,nt1,nt2;
  char fna1[256],fna2[256],fna3[256];

  ifstream swq1;
  swq1.open("cut_direct.par");
  if(!swq1)
    {
		cout<<"cannot open"<<endl;
        abort();
	}
   swq1>>fna1>>fna2>>fna3>>ntrace>>nt1>>nt2;
   swq1.close();   
											
   vector< vector<float> > u1(nt1); 
   vector< vector<float> > u2(nt1); 
   vector< vector<float> > u3(nt2);
   vector< vector<float> > u4(nt1);

   for (i=0;i<nt1;i++)
        {
		    u1[i].resize(ntrace);
		    u2[i].resize(ntrace); 	 
        u4[i].resize(ntrace);
		}  
    for (i=0;i<nt2;i++)
	        u3[i].resize(ntrace);
   
    ifstream swq2;
	swq2.open(fna2,ios::binary);
	if(!swq2)
	   {
	        cout<<"cannot open"<<endl;
	        abort();
	   }

    for(j=0;j<ntrace;j++)
	   for(i=0;i<nt1;i++)
			swq2.read((char*)&u2[i][j],sizeof(u2[i][j]));
    swq2.close();
																	 
    ifstream swq3;
	swq3.open(fna3,ios::binary);
    if(!swq3)
	   {
	       cout<<"cannot open"<<endl;
		   abort();
       }

	for(j=0;j<ntrace;j++)
		for(i=0;i<nt2;i++)
			swq3.read((char*)&u3[i][j],sizeof(u3[i][j]));
	swq3.close();
	
    for(j=0;j<ntrace;j++)
	     for(i=0;i<nt1;i++)
		     {
                u1[i][j]=0.0;
				u4[i][j]=0.0;
			 }

     for(j=0;j<ntrace;j++)
	       for(i=0;i<nt2;i++)
			    u4[i][j]=u3[i][j];

	 for(j=0;j<ntrace;j++)
	       for(i=0;i<nt1;i++)
			    u2[i][j]-=u4[i][j];
					
	 ofstream swq4;
     swq4.open(fna1,ios::binary);

     for(j=0;j<ntrace;j++)
	    for(i=0;i<nt1;i++)
				swq4.write((char*)&u2[i][j],sizeof(u2[i][j]));
	swq4.close();
										   
	return 0;
}














