#ifndef STO_FUNC_H
#define STO_FUNC_H

template <typename REAL>
class Sto_Func
{
  public:
    Sto_Func();
    ~Sto_Func(){};
    REAL tem; // temperature
    REAL mu;  // chemical potential
    REAL* Emin = nullptr;
    REAL* Emax = nullptr;
    void set_E_range(REAL* Emin_in, REAL* Emax_in);

  public:
    REAL root_fd(REAL e) const;
    REAL fd(REAL e) const;
    REAL nroot_fd(REAL e) const;
    REAL nfd(REAL e) const;
    REAL nxfd(REAL e) const;
    REAL fdlnfd(REAL e) const;
    REAL nfdlnfd(REAL e) const;
    REAL n_root_fdlnfd(REAL e) const;
    REAL nroot_mfd(REAL e) const;

  public:
    REAL t;
    REAL ncos(REAL e) const;
    REAL nsin(REAL e) const;
    REAL n_sin(REAL e) const;

  public:
    REAL sigma;
    REAL targ_e;
    REAL gauss(REAL e) const;
    REAL ngauss(REAL e) const;
    REAL nroot_gauss(REAL e) const;
};

#endif