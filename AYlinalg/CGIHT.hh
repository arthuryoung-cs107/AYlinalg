#ifndef CGIHT_HH
#define CGIHT_HH

#include "omp.h"

#include "AYlinalg.hh"

class CGIHT
{
  public:
    CGIHT(AYsym * A_, AYvec * b_);
    ~CGIHT();

    int N, CG_max=0, verbose_capped=0, CG_max_mult=100, CG_count, CG_verbose_div = 1;
    double CG_cond, CG_conv, CG_tol=1.0e-11, precond_thresh=1e-1, alpha, beta;

    AYsym * A;
    AYvec * b;
    bool verbose=false, verbose_done_flag = true, precond_flag=false, restart_flag = true;
    virtual void init(double CG_tol_, int CG_max_mult_);

    virtual void solve_CG(AYvec * x0, AYvec * x_final, bool verbose_ = false);
    virtual void solve_CG(AY_Choleskyspace * space_, AYvec * x0, AYvec * x_final, bool verbose_ = false);

    virtual void solve_CGIHT(int k_, AYvec * x0, AYvec * x_final, bool verbose_ = false);
    virtual void solve_CGIHT(int k_, AY_Choleskyspace * space_, AYvec * x0, AYvec * x_final, bool verbose_ = false);
  protected:
    virtual void verbose_CG();
    virtual void verbose_CGIHT();
    virtual void verbose_done();

    double max(double a, double b) { if (a > b) {return a;} else {return b;} }
    double min(double a, double b) { if (a < b) {return a;} else {return b;} }

    double M_inv(AYsym * L_, AYvec * z_, AYvec *r_, AYvec *w_, AYvec *x_, AYvec *p_, double alpha_); // assumes L_ is already the decomposition
    double M_inv(AYsym * L_, AYvec * z_, AYvec *r_, AYvec * z_old_);
    void proj_T(AYvec * out_, AYvec * t_, int * T_, bool init_0_=true);
    void proj_T(AYvec * out_, AYvec * in_, AYvec * t_, int * T_, bool init_0_=true);
    bool T_compare(int * T1, int * T2, int k);
};


#endif
