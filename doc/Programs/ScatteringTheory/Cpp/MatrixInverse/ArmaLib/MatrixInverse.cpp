// Using armadillo as library

#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

int main()
{
   mat A = randu<mat>(5,5);

  A.print("A before inversion=");
  // find inverse
  mat B = inv(A);
  B.print("A after inversion=");
  // Test of inverse 
  (A*B).print("Is the matrix diagonal?");
  return 0;
}

