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
}

void AYtens_test()
{
  int W=5;
  int M=4;
  int N=3;
  AYtens T1(W, M, N); T1.init_123();
  double *** T1_tpr = AYd3tensor(W, M, N);

  for (int i = 0; i < W*M*N; i++)
  {
    printf("%f\n", T1.T_AT[0][0][i]);
    T1_tpr[0][0][i] = (double) i+1;
  }

  printf("\n");

  int count = 0;
  for (int k = 0; k < W; k++)
  {
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < N; j++)
      {
        printf("%f ", T1_tpr[k][i][j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  printf("\n" );

  for (int i = 0; i < W; i++)
  {
    T1.mat[i].print_mat();
  }
  free_AYd3tensor(T1_tpr);
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

  char name3[] = "./dat_dir/sym3";
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

void AYdiag_test()
{
  int N = 4;
  int M = 5;
  AYmat mat1(M, N); mat1.init_randuni();
  AYdiag diag1(M); diag1.init_123();
  AYsym sym1(N); sym1.init_XTWX(&mat1, &diag1);

  printf("mat1:\n"); mat1.print_mat();
  printf("diag1:\n"); diag1.print_mat();
  printf("sym1:\n"); sym1.print_mat();

  diag1.fprintf_diag("./dat_dir/diag1", true);
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

  char name1[] = "./dat_dir/sym1";
  sym1.fprintf_sym(name1);
}

void sorting_test()
{
  int N = 10;
  int Ntop = 10;
  AYvec vec1(N); vec1.init_randuni();
  AYvec top1(Ntop);
  int * ind1 = new int[Ntop];

  for (int i = 0; i < N; i++) vec1.A_ptr[i] = 0.001*(vec1.A_ptr[i]-0.5);
  printf("vec 1:\n");
  vec1.print_vec(false);
  printf("norm 1: %e\n", vec1.norm_1());


  vec1.max_mag_elements_ordered(&top1, ind1);

  printf("\nvec 1, ordered:\n");
  for (int i = 0; i < Ntop; i++) printf("%d %f\n", ind1[i], top1.A_ptr[i]);

  printf("\nvec 1, projected\n ");
  vec1.Proj_1(&top1, ind1);
  top1.print_vec();
  printf("norm 1: %e\n", top1.norm_1());

  delete ind1 ;
}

void write_test()
{
  int M=5,N=3,W=4;

  AYmat mat1(M, N); mat1.init_123();
  AYvec vec1(M); vec1.init_123();
  AYsym sym1(M); sym1.init_123();
  AYtens tens1(W, M, N); tens1.init_123();

  AYmat mat2(M, N); mat2.init_randuni();
  AYvec vec2(M); vec2.init_randuni();
  AYsym sym2(N); sym2.init_sqrmat(&mat2);
  AYtens tens2(W, M, N); tens2.init_mats123();

  printf("mat1:\n"); mat1.print_mat();
  printf("vec1:\n"); vec1.print_vec();
  printf("sym1:\n"); sym1.print_mat();
  printf("tens1:\n"); tens1.print_tens();
  printf("mat2:\n"); mat2.print_mat();
  printf("vec2:\n"); vec2.print_vec();
  printf("sym2:\n"); sym2.print_mat();
  printf("tens2:\n"); tens2.print_tens();

  mat1.fprintf_mat("./dat_dir/mat1", true);
  vec1.fprintf_vec("./dat_dir/vec1", true);
  sym1.fprintf_sym("./dat_dir/sym1", true);
  tens1.fprintf_tens("./dat_dir/tens1", 1, true);
  mat2.fprintf_mat("./dat_dir/mat2", true);
  vec2.fprintf_vec("./dat_dir/vec2", true);
  sym2.fprintf_sym("./dat_dir/sym2", true);
  tens2.fprintf_tens("./dat_dir/tens2", 1, true);
  tens1.fprintf_tens("./dat_dir/tens1_split", 0, true);
  tens2.fprintf_tens("./dat_dir/tens2_split", 0, true);
}

void read_test()
{
  AYmat mat1("./dat_dir/mat1_b");
  AYvec vec1("./dat_dir/vec1_b");
  AYtens tens1("./dat_dir/tens1_b");
  AYmat mat2("./dat_dir/mat2_b");
  AYvec vec2("./dat_dir/vec2_b");
  AYtens tens2("./dat_dir/tens2_b");
  AYtens tens1_split("./dat_dir/tens1_split_b");
  AYtens tens2_split("./dat_dir/tens2_split_b");

  printf("mat1:\n"); mat1.print_mat();
  printf("vec1:\n"); vec1.print_vec();
  printf("tens1:\n"); tens1.print_tens();
  printf("mat2:\n"); mat2.print_mat();
  printf("vec2:\n"); vec2.print_vec();
  printf("tens2:\n"); tens2.print_tens();
  printf("tens1_split:\n"); tens1_split.print_tens();
  printf("tens2_split:\n"); tens2_split.print_tens();

}

int main()
{
  // preliminary_test1();
  // AYtens_test();
  // AYsym_test();
  // AYdiag_test();
  // Cholesky_test();
  Cholesky_solve_test();
  // sorting_test();
  // write_test();
  // read_test();

  return 0;
}
