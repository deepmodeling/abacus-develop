#include "dfpt_pw.h"
#include "dfpt_pw_impl.h"

#include "dfpt_hamilt_shift.h"
#include "dfpt_kq_basis.h"
#include "dfpt_metal.h"
#include "dfpt_pert.h"
#include "dfpt_phon.h"
#include "dfpt_pw_data.h"
#include "dfpt_q0.h"
#include "dfpt_rho.h"
#include "dfpt_stern.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/qlist.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace ModuleDFPT
{

DFPT_PW::Impl::Impl()
{
}

DFPT_PW::Impl::~Impl()
{
}

bool DFPT_PW::Impl::wired() const
{
    return pw_rho_ != nullptr && pw_wfc_ != nullptr;
}

DFPT_PW::DFPT_PW() : pimpl_(std::unique_ptr<Impl>(new Impl()))
{
}

DFPT_PW::~DFPT_PW()
{
}

bool DFPT_PW::get_with_u() const
{
    return pimpl_->data_.with_u();
}

bool DFPT_PW::get_u_active() const
{
    return pimpl_->data_.u_active();
}

int DFPT_PW::get_nq() const
{
    return pimpl_->qlist_.get_nq();
}

ModuleBase::Vector3<double> DFPT_PW::get_qvec(int q_idx) const
{
    return pimpl_->data_.get_qvec(q_idx);
}

std::vector<double> DFPT_PW::get_phonon_freq(int q_idx) const
{
    return pimpl_->data_.get_phon_freq(q_idx);
}

std::vector<double> DFPT_PW::get_phon_freq_loto() const
{
    return pimpl_->data_.get_phon_freq_loto();
}

ModuleBase::Vector3<double> DFPT_PW::get_loto_dir() const
{
    return pimpl_->data_.get_loto_dir();
}

std::string DFPT_PW::format_q_report(int q_idx) const
{
    return pimpl_->phon_.format_q_report(q_idx, pimpl_->data_);
}

std::string DFPT_PW::format_loto_report() const
{
    return pimpl_->phon_.format_loto_report(pimpl_->data_);
}

ModuleBase::matrix DFPT_PW::get_dielectric_tensor() const
{
    return pimpl_->data_.get_dielectric();
}

ModuleBase::matrix DFPT_PW::get_born_charges(int atom_idx) const
{
    return pimpl_->data_.get_born(atom_idx);
}

void DFPT_PW::set_qfile(const std::string& filename)
{
    pimpl_->qfile_ = filename;
}

void DFPT_PW::set_qmesh(int nqx, int nqy, int nqz)
{
    pimpl_->nqx_ = nqx;
    pimpl_->nqy_ = nqy;
    pimpl_->nqz_ = nqz;
}

void DFPT_PW::set_conv_thr(double thr)
{
    pimpl_->conv_thr_ = thr;
    pimpl_->data_.set_conv_thr(thr);
}

void DFPT_PW::set_max_iter(int max_iter)
{
    pimpl_->max_iter_ = max_iter;
    pimpl_->data_.set_max_iter(max_iter);
}

void DFPT_PW::set_mix_beta(double beta)
{
    if (beta > 0.0 && beta <= 1.0)
    {
        pimpl_->mix_beta_ = beta;
    }
}

void DFPT_PW::set_compute_q0(bool flag)
{
    pimpl_->data_.set_compute_q0(flag);
}

void DFPT_PW::set_loto(bool flag)
{
    pimpl_->data_.set_loto(flag);
}

void DFPT_PW::set_loto_dir(const ModuleBase::Vector3<double>& dir)
{
    pimpl_->data_.set_loto_dir(dir);
}

} // namespace ModuleDFPT
