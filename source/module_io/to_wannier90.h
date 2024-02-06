#ifndef TOWannier90_H
#define TOWannier90_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <unordered_set>

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_psi/psi.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_reduce.h"

class toWannier90
{
  public:
    toWannier90();

    toWannier90(
      const bool &out_wannier_mmn, 
      const bool &out_wannier_amn, 
      const bool &out_wannier_unk, 
      const bool &out_wannier_eig,
      const bool &out_wannier_wvfn_formatted, 
      const std::string &nnkpfile,
      const std::string &wannier_spin
    );
    ~toWannier90();

    void calculate();
    void read_nnkp(const K_Vectors& kv);

    void out_eig(const ModuleBase::matrix& ekb);
    void cal_Amn();
    void cal_Mmn();
    void out_unk();

  protected:
    bool try_read_nnkp(const K_Vectors& kv);

    // Parameters related to k point
    int num_kpts;
    int cal_num_kpts;
    std::vector<std::vector<int>> nnlist;
    std::vector<std::vector<ModuleBase::Vector3<double>>> nncell;
    int nntot = 0;
    int start_k_index = 0;

    // Parameters related to trial orbitals
    int num_wannier; // Number of Wannier orbits
    ModuleBase::Vector3<double> *R_centre = nullptr;
    int *L = nullptr;
    int *m = nullptr;
    int *rvalue = nullptr;
    ModuleBase::Vector3<double> *z_axis = nullptr;
    ModuleBase::Vector3<double> *x_axis = nullptr;
    double *alfa = nullptr;
    int *spin_eig = nullptr; // 'up' state is 1, 'down' state is -1
    ModuleBase::Vector3<double> *spin_qaxis = nullptr; // spin quantisation axis
    std::complex<double> *up_con = nullptr;
    std::complex<double> *dn_con = nullptr;

    // Wannier control parameters
    bool out_wannier_mmn = true;
    bool out_wannier_amn = true;
    bool out_wannier_unk = true;
    bool out_wannier_eig = true;
    bool out_wannier_wvfn_formatted = true;

    std::string nnkpfile = "";
    std::string wannier_file_name = "seedname";
    std::string wannier_spin = "up";

    int num_exclude_bands = 0;
    // int *exclude_bands = nullptr;
    std::unordered_set<int> exclude_bands;
    // bool *tag_cal_band = nullptr;
    int num_bands = 0;
    int *cal_band_index = nullptr;
    bool gamma_only_wannier = false;
    

};

#endif
