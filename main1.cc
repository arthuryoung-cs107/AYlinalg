#include <sys/stat.h>
#include <sstream>
#include <string>
#include <fstream>

#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

void preliminary_test1()
{
  int W=5;
  int M=4;
  int N=3;
  AYtens T1(W, M, N); T1.init_123();
  AYtens T2(W, M, N); T2.init_mats123();
  AYmat M1(M, N); M1.init_123();

  printf("T1\n");
  T1.print_tens();

  printf("T2\n");
  T2.print_tens();

  printf("M1\n");
  M1.print_mat();

  char name1[50]; name_gen(name1, 50, "./dat_dir/tens1");
  char name2[50]; name_gen(name2, 50, "./dat_dir/mat1");
  M1.fprintf_mat(name2);
  T1.fprintf_tens(name1);
}

void preliminary_test2()
{
  char name1_back[50]; name_gen(name1_back, 50, "./dat_dir/tens1_back");
  char name2_back[50]; name_gen(name2_back, 50, "./dat_dir/mat1_back");

  printf("T1 back\n");
  AYtens * T1_back = aysml_read_tens(name1_back);
  T1_back->print_tens();
  printf("M1 back\n");
  AYmat * M1_back = aysml_read(name2_back);
  M1_back->print_mat();
}

void AYsym_test()
{
  int N = 4;
  int M = 5;
  AYsym eye(N); eye.init_eye();
  AYsym sym2(N); sym2.init_123();
  AYsym sym3(N);
  AYvec v1(N); v1.init_123();
  AYvec v2(N); v2.init_0();
  AYvec v3(N); v3.init_123();
  AYvec v4(N); v4.init_123();
  AYvec v5(N); v5.init_0();
  AYmat m1(M, N); m1.init_123();

  printf("v1\n");
  v1.print_vec();
  printf("v2\n");
  v2.print_vec();

  printf("sym2\n");
  sym2.print_mat();

  eye.mult_vec(&v1, &v2);
  printf("v2\n");
  v2.print_vec();

  sym2.mult_vec(&v1, &v3);
  printf("v3\n");
  v3.print_vec();

  sym2.mult_vec(&v1, &v4, true);
  printf("v4\n");
  v4.print_vec();

  printf("v^T A v = %f\n", sym2.vT_A_v(&v1, &v5));
  printf("v5\n");
  v5.print_vec();

  printf("m1\n");
  m1.print_mat();

  sym3.init_sqrmat(&m1);
  printf("sym3\n");
  sym3.print_mat();

  char name3[50]; name_gen(name3, 50, "./dat_dir/sym3");
  sym3.fprintf_sym(name3);
}

void Cholesky_test()
{
  int M = 5, N = 3;
  AYmat sym05(M, N); sym05.init_randuni();
  AYsym sym1(N); sym1.init_sqrmat(&sym05);
  AYvec x(N);
  AYvec b(N); b.init_123();
  sym05.print_mat();
  sym1.print_mat();
  b.print_vec();
  AY_Choleskyspace space(N); space.alloc_workspace();
  space.load_mat(&sym1);

  space.solve_system(&x, &b);
  x.print_mat();

  AYsym L(N);
  AY_Choleskyspace space2(N);
  space2.Cholesky_decomp(&sym1, &L);
  L.print_mat();

}

void Cholesky_solve_test()
{
  int i,j;
  int M = 5, N = 4;
  AYmat sym05(M, N); sym05.init_randuni();
  AYsym sym1(N); sym1.init_sqrmat(&sym05);
  AY_Choleskyspace space(N);
  AYsym * L_= new AYsym(N);
  AYvec * r_= new AYvec(N); r_->init_123();
  AYvec * z_= new AYvec(N);

  space.Cholesky_decomp(&sym1, L_);

  printf("L_:\n");
  L_->print_mat();
  printf("r_:\n");
  r_->print_mat();


  for ( i = 0; i < N; i++)
  {
    double sum = r_->A_ptr[i];
    for ( j = 0; j < i; j++) sum -= z_->A_ptr[j]*L_->A[j][i-j];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
  for ( i = N-1; i >= 0; i--)
  {
    double sum = z_->A_ptr[i];
    for ( j = N-1; j > i; j--) sum -= z_->A_ptr[j]*L_->A[i][j-i];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
  printf("z_:\n");
  z_->print_mat();

  char name1[50]; name_gen(name1, 50, "./dat_dir/sym1");
  sym1.fprintf_sym(name1);

}

int main()
{
  // preliminary_test1();
  // preliminary_test2();
  // AYsym_test();
  // Cholesky_test();
  Cholesky_solve_test();
  return 0;
}
