#include "iostream.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

int main()
{
   srand((int)time(0));
   for(int ix=0;ix<=10;ix++)
      cout<<(float)(rand()/float(RAND_MAX))*0.05-0.025<<endl;

   cout<<endl;

   srand((int)time(0));
   for(int ix=0;ix<=10;ix++)
//      cout<<(int)(rand()/int(RAND_MAX))*50<<endl;
      cout<<rand()%50+25<<endl;
   cout<<"====================="<<endl;   
   
   for(int ix=0;ix<=10;ix++)
      cout<<rand()%50+25<<endl;

   return 0;

}
