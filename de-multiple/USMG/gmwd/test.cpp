#include "iostream.h"
#include "complex.h"

int main()
{
  complex<float> b;

  b.real()=1.0;
  b.imag()=2.0;

  complex<float> c;
  c=1+b;

  cout<<c<<endl;

  return 0;
}
