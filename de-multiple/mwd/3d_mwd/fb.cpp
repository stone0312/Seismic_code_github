#include <iostream>
#include <stdlib.h>
#include <complex>

using namespace std;

//prototypes.
complex <float> *alloc1complex(size_t n1);
void free1complex( complex <float>*p);

complex <float> *alloc1complex(size_t n1) 
{
        return (complex <float>*)calloc(n1,sizeof(complex <float>));
}

/*
// free a 1-d array of complexs 
void free1complex( complex <float>*p)
{
        free1(p);
}
*/
//#include "alloc1.c"

//using namespace std;

int main()
{

int	n1=1000;
int	n2=1000;
complex <float> *u1 =  (complex <float>*)calloc(n1,sizeof(complex <float>));
complex <float> *u2 =  (complex <float>*)alloc1complex(n1);

return 0;
}

