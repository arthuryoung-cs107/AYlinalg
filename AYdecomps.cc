#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

AY_SVDspace::AY_SVDspace(AYmat * mat_): M_in(mat_->M), N_in(mat_->N), U(gsl_matrix_alloc(mat_->M, mat_->N)), V(gsl_matrix_alloc(mat_->N, mat_->N)), s(gsl_vector_alloc(mat_->N)), work(gsl_vector_alloc(mat_->N)) // assume long thin matrix for now
{}
AY_SVDspace::AY_SVDspace(int M_, int N_): M_in(M_), N_in(N_), U(gsl_matrix_alloc(M_, N_)), V(gsl_matrix_alloc(N_, N_)), s(gsl_vector_alloc(N_)), work(gsl_vector_alloc(N_)) // assume long thin matrix for now
{}
AY_SVDspace::~AY_SVDspace() { gsl_matrix_free(U); gsl_matrix_free(V); gsl_vector_free(s); gsl_vector_free(work); }

void AY_SVDspace::load_U(AYmat * X_) {for (int i = 0; i < M_in; i++) {for (int j = 0; j < N_in; j++) {gsl_matrix_set(U, i, j, X_->AT[j][i]);}}}
void AY_SVDspace::svd() {gsl_linalg_SV_decomp(U, V, s, work);}
void AY_SVDspace::unpack(AYmat * U_, AYmat * S_, AYmat * V_)
{
  for (int j = 0; j < N_in; j++)
  {
    for (int i = 0; i < M_in; i++)
    {
      U_->set(i, j, gsl_matrix_get(U, i, j));
      if (i < N_in)
      {
        V_->set(i, j, gsl_matrix_get(V, i, j));
        S_->set(i, j, 0.0);
      }
    }
    S_->set(j, j, gsl_vector_get(s, j));
  }
}

AY_Choleskyspace::AY_Choleskyspace(AYsym * mat_): N_in(mat_->N), mat_gsl(gsl_matrix_alloc(mat_->N, mat_->N))
{load_mat(mat_);}

AY_Choleskyspace::AY_Choleskyspace(int N_): N_in(N_), mat_gsl(gsl_matrix_alloc(N_, N_))
{}

AY_Choleskyspace::~AY_Choleskyspace()
{
  gsl_matrix_free(mat_gsl);
  if (workspace_alloc)
  {
    gsl_vector_free(x_gsl);
  }
}

void AY_Choleskyspace::load_mat(AYsym * mat_)
{
  int i,j;
  for ( i = 0; i < N_in; i++)
  {
    //lower diagonal components and diagonal only, exploiting gsl technique
    for ( j = 0; j < i; j++) gsl_matrix_set(mat_gsl, i, j, mat_->A[j][i-j]);
    //possibly redundant?
    for ( j = i; j < N_in; j++) gsl_matrix_set(mat_gsl, i, j, mat_->A[i][j-i]);
  }
}
void AY_Choleskyspace::load_mat(AYsym * mat_, double scal_)
{
  for (int i = 0; i < N_in; i++)
  {
    //lower diagonal components and diagonal only, exploiting gsl technique
    for (int j = 0; j < i; j++) gsl_matrix_set(mat_gsl, i, j, scal_*mat_->A[j][i-j]);
    //possibly redundant?
    for (int j = i; j < N_in; j++) gsl_matrix_set(mat_gsl, i, j, scal_*mat_->A[i][j-i]);
  }
}

void AY_Choleskyspace::Cholesky_decomp()
{gsl_linalg_cholesky_decomp(mat_gsl);}

void AY_Choleskyspace::alloc_workspace()
{
  workspace_alloc = true;
  x_gsl = gsl_vector_alloc(N_in);
}
void AY_Choleskyspace::solve_system(AYvec * x_in)
{
  gsl_linalg_cholesky_decomp(mat_gsl);
  gsl_linalg_cholesky_svx(mat_gsl, x_gsl);
  x_in->GSL_2_AYvec_copy(x_gsl);
}
void AY_Choleskyspace::solve_system(AYvec * x_in, AYvec * b_in)
{
  gsl_linalg_cholesky_decomp(mat_gsl);
  b_in->AYvec_2_GSL_copy(x_gsl);
  gsl_linalg_cholesky_svx(mat_gsl, x_gsl);
  x_in->GSL_2_AYvec_copy(x_gsl);
}

void AYlinalg_svd(AYmat * mat_, AY_SVDspace * space_) // assume long thin matrix for now
{
  mat_->AYmat_2_GSL_copy(space_->U);
  gsl_linalg_SV_decomp(space_->U, space_->V, space_->s, space_->work);
}
