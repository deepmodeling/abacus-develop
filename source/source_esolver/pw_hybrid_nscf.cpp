#include "pw_hybrid_nscf.h"
#include "esolver_ks_pw.h"
#include "source_io/module_wf/exx_source_io.h"
#include "source_io/module_wf/read_wfc_pw.h"
#include "source_base/module_out/filename.h"
#include "source_hamilt/module_xc/general_exx_info.h"
#include "source_pw/module_pwdft/exx_helper.h"

#include <cerrno>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ModuleESolver
{
namespace
{
std::string exx_restart_configuration(const Input_para& inp, const General_Exx_Info& info)
{
    std::ostringstream out;
    out << std::setprecision(17) << inp.dft_functional << ':' << inp.ecutwfc << ':'
        << info.ecut_exx << ':' << info.hybrid_alpha << ':' << info.gamma_extrapolation;
    for (const auto& kernel : info.coulomb_param)
    {
        out << ':' << static_cast<int>(kernel.first);
        for (const auto& term : kernel.second)
        {
            for (const auto& parameter : term)
            {
                out << ':' << parameter.first << '=' << parameter.second;
            }
        }
    }
    return out.str();
}
}

namespace
{
bool checkpoint_output_enabled(const Input_para& inp)
{
    return inp.out_wfc_pw == 2 && inp.out_freq_ion == 0 && inp.kpar == 1
           && inp.bndpar == 1 && inp.nspin != 4 && inp.symmetry != "1";
}
}

void validate_hybrid_nscf(const Input_para& inp, const General_Exx_Info& info)
{
    if (inp.calculation != "nscf" || !info.cal_exx)
    {
        return;
    }
    const auto fock = info.coulomb_param.find(Conv_Coulomb_Pot_K::Coulomb_Type::Fock);
    const bool has_fock = fock != info.coulomb_param.end() && !fock->second.empty();
    const bool compatible_execution = inp.device == "cpu" && inp.kpar == 1 && inp.bndpar == 1 && inp.nspin != 4;
    const bool compatible_exchange = !inp.exxace && !inp.exx_gamma_extrapolation && !has_fock;
    const bool compatible_targets = inp.symmetry != "1" && inp.init_wfc != "file" && !inp.cal_force && !inp.cal_stress;
    if (!compatible_execution || !compatible_exchange || !compatible_targets)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW",
            "Hybrid NSCF currently requires CPU, kpar=bndpar=1, nspin=1/2, symmetry=-1/0, "
            "exxace=false, exx_gamma_extrapolation=false, screened exchange, "
            "cal_force=cal_stress=false and fresh target wavefunctions");
    }
}

void invalidate_exx_source(const Input_para& inp,
                           const General_Exx_Info& info,
                           const ModulePW::PW_Basis_K& basis,
                           const std::string& out_dir)
{
    if (inp.calculation != "scf" || !info.cal_exx || inp.out_wfc_pw != 2
        || inp.kpar != 1 || inp.bndpar != 1 || basis.poolrank != 0)
    {
        return;
    }
    const std::string checkpoint = out_dir + "EXX_SOURCE";
    if (std::remove(checkpoint.c_str()) != 0 && errno != ENOENT)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "Cannot invalidate the previous EXX source checkpoint");
    }
}

void save_exx_source(const Input_para& inp,
                     const General_Exx_Info& info,
                     const K_Vectors& points,
                     const ModuleBase::matrix& weights,
                     const ModulePW::PW_Basis_K& basis,
                     const std::string& out_dir)
{
    const bool output_enabled = checkpoint_output_enabled(inp);
    if (inp.calculation != "scf" || !info.cal_exx || !output_enabled || basis.poolrank != 0)
    {
        return;
    }
    const std::string filename = out_dir + "EXX_SOURCE";
    const std::string configuration = exx_restart_configuration(inp, info);
    std::ofstream out(filename);
    try
    {
        ModuleIO::write_exx_source(out, points, weights, inp.nspin, basis.nx, basis.ny, basis.nz, configuration);
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", error.what());
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::prepare_exx_nscf(const UnitCell& ucell, const std::string& readin_dir)
{
    ModuleBase::timer::start("ESolver_KS_PW", "prepare_exx_nscf");
    const Input_para& inp = *this->inp_;
    const std::string configuration = exx_restart_configuration(inp, this->general_exx_info_);
    const std::string filename = readin_dir + "EXX_SOURCE";
    std::ifstream in(filename);
    try
    {
        ModuleIO::read_exx_source(in, exx_source_points_, exx_source_weights_, inp.nspin,
                                 this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz, configuration);
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", error.what());
    }
    const int nqs = exx_source_points_.get_nks();
    int nq = exx_source_points_.get_nkstot_nospin();
    exx_source_points_.para_k.kinfo(nq, 1, 0, this->pw_wfc->poolrank, this->pw_wfc->poolnproc, inp.nspin);
    exx_source_basis_.reset(new ModulePW::PW_Basis_K(inp.device, inp.precision));
#ifdef __MPI
    exx_source_basis_->initmpi(this->pw_wfc->poolnproc, this->pw_wfc->poolrank, this->pw_wfc->pool_world);
#endif
    exx_source_basis_->initgrids(ucell.lat0, ucell.latvec,
                               this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz);
    exx_source_basis_->initparameters(false, inp.ecutwfc, nqs, exx_source_points_.kvec_d.data());
    exx_source_basis_->fft_bundle.initfftmode(inp.fft_mode);
    exx_source_basis_->setuptransform();
    exx_source_basis_->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);
    for (int iq = 0; iq < nqs; ++iq)
    {
        exx_source_points_.ngk[iq] = exx_source_basis_->npwk[iq];
        exx_source_points_.kvec_c[iq] = exx_source_points_.kvec_d[iq] * ucell.G;
    }
    const int nbands = exx_source_weights_.nc;
    const int nbasis = exx_source_basis_->npwk_max;
    exx_source_psi_.reset(new psi::Psi<T, Device>(nqs, nbands, nbasis, exx_source_points_.ngk, true));
    for (int iq = 0; iq < nqs; ++iq)
    {
        const std::string wfc_filename = ModuleIO::filename_output(readin_dir, "wf", "pw", iq,
            exx_source_points_.ik2iktot, inp.nspin, nqs, 2, false, false, -1);
        ModuleBase::ComplexMatrix wfc(nbands, nbasis);
        ModuleIO::read_wfc_pw(wfc_filename, exx_source_basis_.get(), this->pw_wfc->poolrank,
                             this->pw_wfc->poolnproc, nbands, 1, iq, iq, nqs, wfc);
        exx_source_psi_->fix_k(iq);
        std::vector<T> coefficients(nbands * nbasis);
        for (int ib = 0; ib < nbands; ++ib)
        {
            for (int ig = 0; ig < nbasis; ++ig)
            {
                coefficients[ib * nbasis + ig] = static_cast<T>(wfc(ib, ig));
            }
        }
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(
            exx_source_psi_->get_pointer(), coefficients.data(), coefficients.size());
    }
    auto* helper = static_cast<Exx_Helper<T, Device>*>(this->exx_helper);
    helper->op_exx->set_source(*exx_source_basis_, exx_source_points_, *exx_source_psi_, exx_source_weights_);
    ModuleBase::timer::end("ESolver_KS_PW", "prepare_exx_nscf");
}


template void ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>::prepare_exx_nscf(
    const UnitCell&, const std::string&);
template void ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::prepare_exx_nscf(
    const UnitCell&, const std::string&);
#if defined(__CUDA) || defined(__ROCM)
template void ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>::prepare_exx_nscf(
    const UnitCell&, const std::string&);
template void ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>::prepare_exx_nscf(
    const UnitCell&, const std::string&);
#endif
} // namespace ModuleESolver
