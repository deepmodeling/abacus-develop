#include "source_lcao/module_lr/lr_core.h"

#include "source_base/tool_title.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_lr/lr_gs.h"

#include <cassert>
#include <complex>

namespace LR
{

template <typename T>
void Core<T>::read_wfc(const std::string& in_dir, const Input_para& inp, const K_Vectors& kv, std::ostream& log)
{
    assert(this->gs.psi != nullptr);
    GsWfc<T> wfc = read_gs_wfc<T>(*this->gs.psi, this->gs.eig, this->gs.wg, this->gs.para_mat, this->gs.para_mat_all,
                                  in_dir, inp, this->dim.nocc_in, this->dim.nvirt_in, this->dim.nocc_max,
                                  this->dim.nbands, this->dim.nbasis, kv, inp.cal_force, log);
    this->gs.psi_all_own = std::move(wfc.all_own);
    this->gs.psi_all = wfc.all;
    this->gs.eig_all = std::move(wfc.eig_all);
    this->gs.wg_all = std::move(wfc.wg_all);
}

template <typename T>
void Core<T>::fill_z(const int* desc_src, const K_Vectors& kv, const bool cal_force, const bool first,
                     std::ostream& log)
{
    ModuleBase::TITLE("LR::Core", "fill_z");
    ZWin<T> z = build_z_win<T>(this->gs.psi_all, this->gs.eig_all, this->dim.nocc, this->dim.nvirt, this->dim.nocc_max,
                               this->dim.nbands, this->dim.nbasis, this->dim.nk, this->dim.openshell, kv,
                               this->amp.para_c, this->gs.para_mat, desc_src, cal_force, first, log);
    this->zw.psi = std::move(z.psi);
    this->zw.eig = std::move(z.eig);
    this->zw.para_c = std::move(z.para_c);
    this->zw.para_x = std::move(z.para_x);
    this->zw.nvirt = std::move(z.nvirt);
    this->zw.nbands = z.nbands;
    this->zw.nloc = z.nloc;
}

template void Core<double>::read_wfc(const std::string&, const Input_para&, const K_Vectors&, std::ostream&);
template void Core<std::complex<double>>::read_wfc(const std::string&, const Input_para&, const K_Vectors&,
                                                   std::ostream&);
template void Core<double>::fill_z(const int*, const K_Vectors&, bool, bool, std::ostream&);
template void Core<std::complex<double>>::fill_z(const int*, const K_Vectors&, bool, bool, std::ostream&);

} // namespace LR
