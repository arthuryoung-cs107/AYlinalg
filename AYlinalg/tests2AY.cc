#include "AYlinalg.hh"
#include "CGIHT.hh"

extern "C"
{
  #include "AYaux.h"
}

void CGIHT_vs_Cholesky_test()
{
  int M = 5, N = 4;
  AYmat sym05(M, N); sym05.init_randuni();
  AYsym sym1(N); sym1.init_sqrmat(&sym05);
  AY_Choleskyspace space(N);

  space.Cholesky_decomp(&sym1, L_);

  printf("L_:\n");
  L_->print_mat();
  printf("r_:\n");
  r_->print_mat();



}

int main()
{
  CGIHT_vs_Cholesky_test();
  return 0;
}
