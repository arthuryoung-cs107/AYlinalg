#include "CGIHT.hh"

extern "C"
{
  #include "AYaux.h"
}

CGIHT::CGIHT(AYsym * A_, AYvec * b_): A(A_), b(b_), N(A_->N)
{}
CGIHT::~CGIHT()
{}
void CGIHT::init(double CG_tol_, int CG_max_mult_)
{CG_tol = CG_tol_; CG_max_mult = CG_max_mult_;}
void CGIHT::solve_CG(AYvec * x0, AYvec * x_final, bool verbose_)
{
  verbose = verbose_; precond_flag = false;
  CG_max = (CG_max==0) ? CG_max_mult*N : CG_max;
  if (verbose_capped) CG_verbose_div  = (int)(CG_max/verbose_capped);
  bool CG_cont=true;
  double inner_rk, inner_rk_old, l2_rk, inner_pTAp, CG_cond0;

  AYvec pk(N), work(N), xk(N), rk(N);

  for (int i = 0; i < N; i++)
  {
    pk.A_ptr[i] = 0.0;
    xk.A_ptr[i] = x0->A_ptr[i];
    rk.A_ptr[i] = b->A_ptr[i];
  }
  A->mult_vec(&xk, &rk, true); // initialize residual
  CG_cond0 = inner_rk = l2_rk = rk.norm_2(); inner_rk *= l2_rk; beta = 0.0; CG_count = 0;
  do
  {
    // compute orthogonalization weight
    if (CG_count++ != 0) beta = inner_rk/inner_rk_old;

    //pk = rk + beta pk
    for (int i = 0; i < N; i++)
    {
      pk.A_ptr[i] = beta*(pk.A_ptr[i]) + rk.A_ptr[i];
      work.A_ptr[i] = 0.0;
    }

    // compute alpha:
    inner_pTAp = A->vT_A_v(&pk, &work); // work = Apk
    alpha = (inner_rk)/(inner_pTAp);

    inner_rk_old = inner_rk; inner_rk = 0.0;
    for (int i = 0; i < N; i++)
    {
      rk.A_ptr[i] -= alpha*(work.A_ptr[i]);
      inner_rk += (rk.A_ptr[i])*(rk.A_ptr[i]);
      xk.A_ptr[i] += alpha*(pk.A_ptr[i]);
    }
    l2_rk = sqrt(inner_rk); CG_cond = l2_rk; // error criterion

    verbose_CG();
    if ((CG_cond < CG_tol) || (CG_count >= CG_max)) CG_cont = false;
  } while (CG_cont);
  xk.copy_set(x_final);
  verbose_done();
}
void CGIHT::solve_CG(AY_Choleskyspace * space, AYvec * x0, AYvec * x_final, bool verbose_)
{
  verbose = verbose_; precond_flag = true;
  CG_max = (CG_max==0) ? CG_max_mult*N : CG_max;
  if (verbose_capped) CG_verbose_div = (int)(CG_max/verbose_capped);
  bool CG_cont=true;
  double l2_rk, inner_pTAp, CG_cond0;
  double inner_rz, inner_rz_old;

  AYvec pk(N), work(N), xk(N), rk(N);
  AYvec zk(N);
  AYsym L(N);

  space->iCholesky_decomp(A, &L, precond_thresh);

  for (int i = 0; i < N; i++)
  {
    xk.A_ptr[i] = x0->A_ptr[i];
    rk.A_ptr[i] = b->A_ptr[i];
  }
  A->mult_vec(&xk, &rk, true); // initialize residual
  AYlinalg_Cholesky_solve(&L, &zk, &rk); // initialize zk
  inner_rz=rk.inner(&zk);
  CG_cond0 = l2_rk = rk.norm_2(); beta = 0.0; CG_count = 0;
  do
  {
    // compute orthogonalization weight
    if (CG_count++ != 0) beta = inner_rz/inner_rz_old;

    //pk = rk + beta pk
    for (int i = 0; i < N; i++)
    {
      pk.A_ptr[i] = beta*(pk.A_ptr[i]) + zk.A_ptr[i];
      work.A_ptr[i] = 0.0;
    }

    // compute alpha:
    inner_pTAp = A->vT_A_v(&pk, &work); // work = Apk
    alpha = (inner_rz)/(inner_pTAp);
    inner_rz_old = inner_rz;

    // this step will simultaneously update rk, xk, zk, inner_rz, and CG_cond
    inner_rz = M_inv(&L, &zk, &rk, &work, &xk, &pk, alpha);
    l2_rk = CG_cond;

    verbose_CG();
    if ((CG_cond < CG_tol) || (CG_count >= CG_max)) CG_cont = false;
  } while (CG_cont);
  xk.copy_set(x_final);
  verbose_done();
}

void CGIHT::solve_CGIHT(int k_, AYvec * x0, AYvec * x_final, bool verbose_)
{
  verbose = verbose_; precond_flag = false;
  CG_max = (CG_max==0) ? CG_max_mult*N : CG_max;
  if (verbose_capped) CG_verbose_div = (int)(CG_max/verbose_capped);
  bool CG_cont=true;
  double inner_rk, inner_rkproj, inner_rkproj_old, l2_rk, inner_pTAp, CG_cond0;
  int k = k_;

  AYvec pk(N), work(N), xk(N), rk(N), rk_old(N);
  AYvec tk(k), wk(N);
  int *Tk = new int[k], *Tk_old = new int[k];
  int i_small;

  for (int i = 0; i < k; i++) Tk_old[i] = -1;
  for (int i = 0; i < N; i++)
  {
    pk.A_ptr[i] = xk.A_ptr[i] = 0.0;
    wk.A_ptr[i] = x0->A_ptr[i];
    rk.A_ptr[i] = rk_old.A_ptr[i] = b->A_ptr[i];
  }
  i_small = wk.max_mag_elements(&tk, Tk);
  proj_T(&xk, &tk, Tk, false);
  A->mult_vec(&xk, &rk, true); // initialize residual
  CG_cond0 = l2_rk = rk.norm_2(); CG_count = 0;
  do
  {
    inner_rkproj = inner_rkproj_old = 0.0;
    for (int i = 0; i < k; i++)
    {
      inner_rkproj += (rk.A_ptr[Tk[i]])*(rk.A_ptr[Tk[i]]);
      inner_rkproj_old += (rk_old.A_ptr[Tk[i]])*(rk_old.A_ptr[Tk[i]]);
    }
    // compute orthogonalization weight, evaluating restart condition
    if (restart_flag) beta = (T_compare(Tk, Tk_old, k)) ? inner_rkproj/inner_rkproj_old : 0.0;
    else beta = (CG_count!=0) ? inner_rkproj/inner_rkproj_old : 0.0;

    //pk = rk + beta pk
    for (int i = 0; i < N; i++)
    {
      pk.A_ptr[i] = beta*(pk.A_ptr[i]) + rk.A_ptr[i];
      work.A_ptr[i] = wk.A_ptr[i] = 0.0;
    }
    proj_T(&work, &pk, &tk, Tk, false);
    inner_pTAp = A->vT_A_v(&work, &wk); // work = Apk
    // compute alpha:
    alpha = (inner_rkproj)/(inner_pTAp);
    for (int i = 0; i < N; i++)
    {
      wk.A_ptr[i] = xk.A_ptr[i] + alpha*(pk.A_ptr[i]);
      xk.A_ptr[i] = 0.0;
      rk_old.A_ptr[i] = rk.A_ptr[i];
      rk.A_ptr[i] = b->A_ptr[i];
    }
    for (int i = 0; i < k; i++) Tk_old[i] = Tk[i];
    i_small = wk.max_mag_elements(&tk, Tk);
    proj_T(&xk, &tk, Tk, false);
    A->mult_vec(&xk, &rk, true); // new residual
    CG_cond = l2_rk = rk.norm_2(); // error criterion

    if ((CG_cond < CG_tol) || (++CG_count >= CG_max)) CG_cont = false;
    verbose_CG();
  } while (CG_cont);
  xk.copy_set(x_final);
  verbose_done();
  delete Tk;
  delete Tk_old;
}

void CGIHT::solve_CGIHT(int k_, AY_Choleskyspace * space, AYvec * x0, AYvec * x_final, bool verbose_)
{
  verbose = verbose_; precond_flag = true;
  CG_max = (CG_max==0) ? CG_max_mult*N : CG_max;
  if (verbose_capped) CG_verbose_div = (int)(CG_max/verbose_capped);
  bool CG_cont=true;
  double inner_rk, inner_rkproj, inner_rkproj_old, l2_rk, inner_pTAp, CG_cond0;
  double inner_rz, inner_rz_old;
  int k = k_;

  AYvec pk(N), work(N), xk(N), rk(N), rk_old(N);
  AYvec tk(k), wk(N);
  AYvec zk(N), zk_old(N);
  AYsym L(N);

  space->iCholesky_decomp(A, &L, precond_thresh);

  int *Tk = new int[k], *Tk_old = new int[k];
  int i_small;

  for (int i = 0; i < k; i++) Tk_old[i] = -1;

  for (int i = 0; i < N; i++)
  {
    pk.A_ptr[i] = xk.A_ptr[i] = 0.0;
    wk.A_ptr[i] = x0->A_ptr[i];
    rk.A_ptr[i] = rk_old.A_ptr[i] = zk_old.A_ptr[i] = b->A_ptr[i];
  }
  i_small = wk.max_mag_elements(&tk, Tk);
  proj_T(&xk, &tk, Tk, false);
  A->mult_vec(&xk, &rk, true); // initialize residual
  AYlinalg_Cholesky_solve(&L, &zk, &rk); // initialize zk

  CG_cond0 = l2_rk = rk.norm_2(); CG_count = 0;
  do
  {
    inner_rz = inner_rz_old = 0.0;
    for (int i = 0; i < k; i++)
    {
      inner_rz += (rk.A_ptr[Tk[i]])*(zk.A_ptr[Tk[i]]);
      inner_rz_old += (rk_old.A_ptr[Tk[i]])*(zk_old.A_ptr[Tk[i]]);
    }
    // compute orthogonalization weight, evaluating restart condition
    if (restart_flag) beta = (T_compare(Tk, Tk_old, k)) ? inner_rz/inner_rz_old : 0.0;
    else beta = (CG_count!=0) ? inner_rz/inner_rz_old : 0.0;


    //pk = rk + beta pk
    for (int i = 0; i < N; i++)
    {
      pk.A_ptr[i] = beta*(pk.A_ptr[i]) + zk.A_ptr[i];
      work.A_ptr[i] = wk.A_ptr[i] = 0.0; // initializing
    }
    proj_T(&work, &pk, &tk, Tk, false);
    inner_pTAp = A->vT_A_v(&work, &wk); // work = Apk
    // compute alpha:
    alpha = (inner_rz)/(inner_pTAp);
    for (int i = 0; i < N; i++)
    {
      wk.A_ptr[i] = xk.A_ptr[i] + alpha*(pk.A_ptr[i]);
      xk.A_ptr[i] = 0.0;
      rk_old.A_ptr[i] = rk.A_ptr[i];
      rk.A_ptr[i] = b->A_ptr[i];
    }
    for (int i = 0; i < k; i++) Tk_old[i] = Tk[i];
    i_small = wk.max_mag_elements(&tk, Tk);
    proj_T(&xk, &tk, Tk, false);
    A->mult_vec(&xk, &rk, true); // new residual
    CG_cond = l2_rk = M_inv(&L, &zk, &rk, &zk_old);

    if ((CG_cond < CG_tol) || (++CG_count >= CG_max)) CG_cont = false;
    verbose_CG();
  } while (CG_cont);
  xk.copy_set(x_final);
  verbose_done();
  delete Tk;
  delete Tk_old;
}

double CGIHT::M_inv(AYsym * L_, AYvec * z_, AYvec *r_, AYvec *w_, AYvec *x_, AYvec *p_, double alpha_) // assumes L_ is already the decomposition
{
  CG_cond = 0.0;
  for (int i = 0; i < N; i++)
  {
    x_->A_ptr[i] += alpha_*(p_->A_ptr[i]); // updating x
    double sum = r_->A_ptr[i] -= alpha_*(w_->A_ptr[i]); // initialized to residual value, which we are updating in the same loop
    CG_cond += sum*sum; // updating the residual norm in the same step
    for ( int j = 0; j < i; j++) sum -= z_->A_ptr[j]*L_->A[j][i-j];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
  CG_cond = sqrt(CG_cond); // residual updated
  double inner_rz_out = 0.0; // next, determine rz inner product with known z
  for (int i = N-1; i >= 0; i--)
  {
    double sum = z_->A_ptr[i];
    for (int j = N-1; j > i; j--) sum -= z_->A_ptr[j]*L_->A[i][j-i];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
    inner_rz_out += (z_->A_ptr[i])*(r_->A_ptr[i]);
  }
  return inner_rz_out;
}

double CGIHT::M_inv(AYsym * L_, AYvec * z_, AYvec *r_, AYvec * z_old_) // assumes L_ is already the decomposition
{
  double res_out = 0.0;
  int i, j, N = L_->N;
  for ( i = 0; i < N; i++)
  {
    z_old_->A_ptr[i] = z_->A_ptr[i];
    double sum = r_->A_ptr[i];
    res_out += sum*sum;
    for ( j = 0; j < i; j++) sum -= z_->A_ptr[j]*L_->A[j][i-j];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
  for ( i = N-1; i >= 0; i--)
  {
    double sum = z_->A_ptr[i];
    for ( j = N-1; j > i; j--) sum -= z_->A_ptr[j]*L_->A[i][j-i];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
  return sqrt(res_out);
}

void CGIHT::proj_T(AYvec *out_, AYvec * t_, int * T_, bool init_0_)
{
  if (init_0_) out_->init_0();
  for (int i = 0; i < t_->M; i++) out_->A_ptr[T_[i]] = t_->A_ptr[i];
}

void CGIHT::proj_T(AYvec *out_, AYvec * in_, AYvec * t_, int * T_, bool init_0_)
{
  if (init_0_) out_->init_0();
  for (int i = 0; i < t_->M; i++) out_->A_ptr[T_[i]] = in_->A_ptr[T_[i]];
}

bool CGIHT::T_compare(int * T1_, int * T2_, int k_)
{
  int sum1=0, sum2=0;
  for (int i = 0; i < k_; i++)
  {
    sum1+=T1_[i]; sum2+=T2_[i];
  }
  if (sum1==sum2) // they might be the same
  {
    bool return_val=true;
    for (int i = 0; i < k_; i++)
    {
      int check = T1_[i];
      bool found = false;
      for (int j = 0; j < k_; j++) if (check == T2_[j]) {found = true; break;}
      if (!found) {return_val=false; break;}
    }
    return return_val;
  }
  else return false;
}

void CGIHT::verbose_done()
{
  if (verbose_done_flag) printf("      (thread %d of %d) CG: Complete. Ran %d of %d CG_it. Condition: %e, convergence completion: %f%% \n", omp_get_thread_num(), omp_get_max_threads(), CG_count, CG_max, CG_cond, 100.0*CG_tol/CG_cond);
}

void CGIHT::verbose_CG()
{ if (verbose && ((CG_count % CG_verbose_div) == 0) ) printf("        (thread %d of %d) CG: CG_count: %d of %d, beta = %e, err = %e, convergence: %f%% \n", omp_get_thread_num(), omp_get_max_threads(), CG_count, CG_max, beta, CG_cond, 100.0*CG_tol/CG_cond); }

void CGIHT::verbose_CGIHT()
{ if (verbose && ((CG_count % CG_verbose_div) == 0) ) printf("        (thread %d of %d) CGIHT: CG_count: %d of %d, beta = %e, err = %e, convergence: %f%% \n", omp_get_thread_num(), omp_get_max_threads(), CG_count, CG_max, beta, CG_cond, 100.0*CG_tol/CG_cond); }
