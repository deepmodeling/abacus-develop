#ifndef TO_W90_H
#define TO_W90_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <unordered_set>

#include "source_base/complexmatrix.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_cell/klist.h"

class toW90
{
  public:
    toW90();

    toW90(
      const bool &out_wannier_mmn, 
      const bool &out_wannier_amn, 
      const bool &out_wannier_unk, 
      const bool &out_wannier_eig,
      const bool &out_wannier_wvfn_formatted, 
      const std::string &nnkpfile,
      const std::string &wannier_spin,
      const int &nspin,
      const int &nbands,
      const int &nqx,
      const double &dq,
      const int &npol
    );
    ~toW90();

    void calculate();
    void read_nnkp(const UnitCell& ucell, const K_Vectors& kv);

    void out_eig(const ModuleBase::matrix& ekb);
    void cal_Amn();
    void cal_Mmn();
    void out_unk();

  protected:
    bool try_read_nnkp(const UnitCell& ucell, const K_Vectors& kv);

    // Parameters related to k point
    int num_kpts=0;
    int cal_num_kpts=0;
    std::vector<std::vector<int>> nnlist;
    std::vector<std::vector<ModuleBase::Vector3<double>>> nncell;
    int nntot = 0;
    int start_k_index = 0;

    // Parameters related to trial orbitals
    int num_wannier=0; // Number of Wannier orbits
    std::vector<ModuleBase::Vector3<double>> R_centre;
    std::vector<int> L;
    std::vector<int> m;
    std::vector<int> rvalue;
    std::vector<ModuleBase::Vector3<double>> z_axis;
    std::vector<ModuleBase::Vector3<double>> x_axis;
    std::vector<double> alfa;
    std::vector<int> spin_eig; // 'up' state is 1, 'down' state is -1
    std::vector<ModuleBase::Vector3<double>> spin_qaxis; // spin quantisation axis
    std::vector<std::complex<double>> up_con;
    std::vector<std::complex<double>> dn_con;

    // Wannier control parameters
    bool out_wannier_mmn = true;
    bool out_wannier_amn = true;
    bool out_wannier_unk = true;
    bool out_wannier_eig = true;
    bool out_wannier_wvfn_formatted = true;

    std::string nnkpfile = "";
    std::string wannier_file_name = "seedname";
    std::string wannier_spin = "up";

    int nspin_ = 1; // number of spin states, 1, 2 or 4
    int nbands_ = 0; // number of bands from input
    int nqx_ = 0; // number of points in q-space integration
    double dq_ = 0.0; // step of dq for q-space integration
    int npol_ = 1; // number of polar states

    int num_exclude_bands = 0;
    std::unordered_set<int> exclude_bands;
    int num_bands = 0;
    std::vector<int> cal_band_index;
    bool gamma_only_wannier = false;
    

};

#endif
