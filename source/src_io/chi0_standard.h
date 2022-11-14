#ifndef CHI0_STANDARD_H
#define CHI0_STANDARD_H

#include "../src_pw/wavefunc.h"
#include "../src_parallel/parallel_global.h"
#include "../module_pw/pw_basis.h"
#include "module_elecstate/elecstate.h"

class Chi0_standard
{
public:

      Chi0_standard();
	  ~Chi0_standard();
	  
	  bool epsilon;
	  std::string system;
	  double eta;
	  double domega;
	  int nomega;
	  int dim;
	  int oband;
	  double q_start[3];
	  double direct[3];
	  int start_q;
	  int interval_q;
	  int nq;
	  bool out_epsilon;
	  
	  void Chi(const elecstate::ElecState* pelec);
	  void Parallel_G();
	  void Init();
	  void Delete();
	  void Cal_Psi(int iq, std::complex<double> **psi_r);
	  void Cal_b(int iq, int ik, int iqk,  ModulePW::PW_Basis *rho_basis);
	  void Cal_weight(int iq, int ik, double omega, const elecstate::ElecState* pelec);
	  void Cal_last();
	  void Cal_first();
	  void Cal_chi0(int iq, double omega, const elecstate::ElecState* pelec);
	  void Cal_rpa(int iq);
	  void Cal_chi();
	  double qg2( int iq, int g0);
	  int Cinv(int n, std::complex<double>** a);
	  void CMatrixMul(int m, int n, int l, std::complex<double>** A, std::complex<double>** B, std::complex<double>** C);
	  int Cal_iq(int ik, int iq, int a, int b, int c);
	  int parallel_g();
	  void chi0_para_g();
	  
	  ModuleBase::Vector3<double> *all_gcar;	  
	  std::complex<double> **chi0;
	  std::complex<double> **chi;
	  
	  
private:
	  
	  bool init_finish;
	  std::complex<double> **psi_r1;
	  std::complex<double> **psi_r2;
	  std::complex<double> **b;
	  std::complex<double> **A;
	  std::complex<double> **B;
	  std::complex<double> *weight;
	  std::complex<double> **rpa;
	  std::complex<double> *b_core;
	  int *num_G_core;
	  int *num_G_dis;
	  std::complex<double> *b_summary;
	  std::complex<double> *b_order;
	  double *G_r_core;
	  int *num_Gvector_core;
	  int *num_Gvector_dis;
	  int *flag;
	  double *G_r;
	  double *Gvec_core;
	  double *Gvec;
	  int *flag1;
	  double **para_g;
	  int dim_para;
	  std::complex<double> **chi0_para;
	  	
};

namespace GlobalC
{
extern Chi0_standard chi0_standard;
}

#endif
