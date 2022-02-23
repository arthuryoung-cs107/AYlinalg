#include "AYlinalg.hh"
#include "CGIHT.hh"

extern "C"
{
  #include "AYaux.h"
}
void iCholesky_test()
{
  int M = 20, N = 10;
  double low = -10.0, high = 10.0, thresh = 1e-1;
  AYmat sym05(M, N); sym05.init_randuni(low, high);
  AYsym sym1(N); sym1.init_sqrmat(&sym05);
  AY_Choleskyspace space1(N), space2(N);
  AYsym L1(N), L2(N);

  sym1.print_mat("sym1");

  space1.Cholesky_decomp(&sym1, &L1); L1.print_mat("L1");
  space2.iCholesky_decomp(&sym1, &L2, thresh); L2.print_mat("L1");

}
void CGIHT_vs_Cholesky_test()
{
  int N_sparse=3; int s_indices[] = {2, 4, 7};

  int M = 20, N = 10;
  double low = -10.0, high = 10.0;
  AYmat sym05(M, N); sym05.init_randuni(low, high);
  AYsym sym1(N); sym1.init_sqrmat(&sym05);
  AYvec xrand(N_sparse), xtrue(N), x0(N), x1(N), x2(N), x3(N), x4(N), x5(N), b1(N);
  AY_Choleskyspace space1(N), space2(N), space3(N);

  xrand.init_randuni(low, high); xtrue.init_0(); x0.init_0(); b1.init_0();

  for (int i = 0; i < N_sparse; i++) xtrue.A_ptr[s_indices[i]] = xrand.A_ptr[i];

  sym1.mult_vec(&xtrue, &b1);

  sym05.fprintf_mat("./dat_dir/sym05", true); sym1.fprintf_sym("./dat_dir/sym1", true); xtrue.fprintf_vec("./dat_dir/xtrue", true); b1.fprintf_vec("./dat_dir/b1", true);

  xtrue.print_vec("xtrue");

  space1.solve_system(&sym1, &x1, &b1); x1.print_vec("x1"); printf("x1 r2: %e\n", x1.norm_2(&xtrue));

  CGIHT cg1(&sym1, &b1); cg1.solve_CG(&x0, &x2, false); x2.print_vec("x2"); printf("x2 r2: %e\n", x2.norm_2(&xtrue));

  CGIHT cg2(&sym1, &b1); cg2.solve_CG(&space2, &x0, &x3, false); x3.print_vec("x3"); printf("x3 r2: %e\n", x3.norm_2(&xtrue));

  CGIHT cg3(&sym1, &b1); cg3.solve_CGIHT(N_sparse, &x0, &x4, false); x4.print_vec("x4"); printf("x4 r2: %e\n", x4.norm_2(&xtrue));

  CGIHT cg4(&sym1, &b1); cg4.solve_CGIHT(N_sparse, &space3, &x0, &x5, false); x5.print_vec("x5"); printf("x5 r2: %e\n", x5.norm_2(&xtrue));


}

int main()
{
  // iCholesky_test();
  CGIHT_vs_Cholesky_test();
  return 0;
}
