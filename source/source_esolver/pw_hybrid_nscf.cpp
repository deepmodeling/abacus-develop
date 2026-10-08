#include "pw_hybrid_nscf.h"
#include "esolver_ks_pw.h"
#include "pw_hybrid_source.h"
#include "source_io/module_wf/read_wfc_pw.h"
#include "source_base/module_out/filename.h"
#include "source_hamilt/module_xc/general_exx_info.h"
#include "source_pw/module_pwdft/exx_helper.h"

#include <fstream>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace ModuleESolver
{
namespace
{
// Check the stored Miller mapping before delegating coefficient loading to the unchanged IO reader.
void validate_source_grid(const std::string& filename, const ModulePW::PW_Basis_K& basis, int iq)
{
    std::ifstream in(filename, std::ios::binary);
    const auto header = read_exx_source_header(in);
    in.seekg(160);
    int record = 0;
    in.read(reinterpret_cast<char*>(&record), sizeof(record));
    if (record != 12LL * header.npw)
    {
        throw std::runtime_error("Invalid EXX source Miller record");
    }
    const int grid[3] = {basis.nx, basis.ny, basis.nz};
    std::set<long long> indices;
    for (int ig = 0; ig < header.npw; ++ig)
    {
        int miller[3] = {};
        ModuleBase::Vector3<double> direct_g;
        for (int axis = 0; axis < 3; ++axis)
        {
            in.read(reinterpret_cast<char*>(&miller[axis]), sizeof(int));
            if (!in || miller[axis] < 0 || miller[axis] >= grid[axis])
            {
                throw std::runtime_error("EXX source Miller index is incompatible with FFT grid");
            }
            direct_g[axis] = miller[axis] <= grid[axis] / 2 ? miller[axis] : miller[axis] - grid[axis];
        }
        const auto gk = direct_g * basis.G + basis.kvec_c[iq];
        const long long index = (static_cast<long long>(miller[0]) * basis.ny + miller[1]) * basis.nz + miller[2];
        if (gk.norm2() > basis.gk_ecut * (1.0 + 1e-10) || !indices.insert(index).second)
        {
            throw std::runtime_error("EXX source Miller mapping is incompatible with FFT basis");
        }
    }
}

void warn_source_configuration(const Input_para& inp, const std::string& readin_dir)
{
    std::map<std::string, std::string> expected;
    expected["dft_functional"] = inp.dft_functional;
    expected["exx_singularity_correction"] = inp.exx_singularity_correction;
    const std::vector<std::pair<std::string, std::vector<std::string>>> terms = {
        {"exx_erfc_alpha", inp.exx_erfc_alpha},
        {"exx_erfc_omega", inp.exx_erfc_omega},
        {"exx_fock_alpha", inp.exx_fock_alpha}};
    for (const auto& term : terms)
    {
        std::ostringstream value;
        for (const auto& number : term.second)
        {
            value << number << ' ';
        }
        expected[term.first] = value.str();
    }
    std::ostringstream cutoff;
    cutoff << inp.ecutexx;
    expected["ecutexx"] = cutoff.str();
    expected["exx_gamma_extra"] = inp.exx_gamma_extra ? "1" : "0";
    std::ifstream input(readin_dir + "INPUT.info");
    std::string line;
    std::map<std::string, std::string> saved;
    while (std::getline(input, line))
    {
        line = line.substr(0, line.find('#'));
        std::istringstream record(line);
        std::string key;
        record >> key;
        if (expected.count(key))
        {
            std::getline(record, saved[key]);
        }
    }
    for (const auto& item : expected)
    {
        if (!saved.count(item.first))
        {
            ModuleBase::WARNING("ESolver_KS_PW", "EXX source configuration is incomplete or unavailable in INPUT.info; "
                "continuing with existing wavefunctions and eig_occ.txt. Verify the source SCF settings.");
            return;
        }
        std::istringstream original(saved[item.first]);
        std::istringstream current(item.second);
        std::string source_value;
        std::string target_value;
        bool mismatch = false;
        while (original >> source_value)
        {
            if (!(current >> target_value))
            {
                mismatch = true;
                break;
            }
            if (source_value != target_value)
            {
                mismatch = true;
            }
        }
        if (current >> target_value)
        {
            mismatch = true;
        }
        if (mismatch)
        {
            const std::string message = "EXX source configuration differs for " + item.first
                + "; continuing with the frozen SCF source and current NSCF exchange settings. "
                  "Use matching settings for consistent SCF/NSCF bands.";
            ModuleBase::WARNING("ESolver_KS_PW", message);
        }
    }
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
    bool compatible_device = inp.device == "cpu";
#ifdef __CUDA
    compatible_device = compatible_device || inp.device == "gpu";
#endif
    const bool compatible_execution = compatible_device && inp.kpar == 1 && inp.bndpar == 1 && inp.nspin != 4;
    const bool compatible_exchange = !inp.exxace && !inp.exx_gamma_extra && !has_fock;
    const bool compatible_targets = inp.symmetry != "1" && inp.init_wfc != "file" && !inp.cal_force && !inp.cal_stress;
    if (!compatible_execution || !compatible_exchange || !compatible_targets)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW",
            "Hybrid NSCF currently requires CPU or CUDA GPU, kpar=bndpar=1, nspin=1/2, symmetry=-1/0, "
            "exxace=false, exx_gamma_extra=false, screened exchange, "
            "cal_force=cal_stress=false and fresh target wavefunctions");
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::prepare_exx_nscf(const UnitCell& ucell, const std::string& readin_dir)
{
    ModuleBase::timer::start("ESolver_KS_PW", "prepare_exx_nscf");
    const Input_para& inp = *this->inp_;
    warn_source_configuration(inp, readin_dir);
    const int spin_mult = inp.nspin == 2 ? 2 : 1;
    const std::vector<int> first_index(1, 0);
    const std::string first_filename = ModuleIO::filename_output(readin_dir, "wf", "pw", 0,
        first_index, inp.nspin, spin_mult, 2, false, false, -1);
    std::vector<ExxSourceHeader> headers;
    try
    {
        std::ifstream first(first_filename, std::ios::binary);
        const auto header = read_exx_source_header(first);
        const int nks = header.nks;
        if (nks % spin_mult != 0)
        {
            throw std::runtime_error("Invalid EXX source spin dimensions");
        }
        exx_source_points_.set_nks(nks);
        exx_source_points_.set_nkstot(nks);
        const int nq = nks / spin_mult;
        exx_source_points_.set_nkstot_nospin(nq);
        exx_source_points_.set_spin_mult(spin_mult);
        exx_source_points_.kvec_d.resize(nks);
        exx_source_points_.kvec_c.resize(nks);
        exx_source_points_.wk.resize(nks);
        exx_source_points_.isk.resize(nks);
        exx_source_points_.ngk.resize(nks);
        exx_source_points_.ik2iktot.resize(nks);
        for (int iq = 0; iq < nks; ++iq)
        {
            exx_source_points_.ik2iktot[iq] = iq;
        }
        const ModuleBase::Matrix3 inverse_reciprocal = ucell.G.Inverse();
        for (int iq = 0; iq < nks; ++iq)
        {
            const std::string filename = ModuleIO::filename_output(readin_dir, "wf", "pw", iq,
                exx_source_points_.ik2iktot, inp.nspin, nks, 2, false, false, -1);
            std::ifstream in(filename, std::ios::binary);
            const auto metadata = read_exx_source_header(in);
            if (metadata.ecutwfc != inp.ecutwfc)
            {
                throw std::runtime_error("EXX source wavefunction cutoff differs from the current basis");
            }
            headers.push_back(metadata);
            auto& q = exx_source_points_.kvec_c[iq];
            q.x = metadata.kvec_c[0];
            q.y = metadata.kvec_c[1];
            q.z = metadata.kvec_c[2];
            exx_source_points_.kvec_d[iq] = q * inverse_reciprocal;
            exx_source_points_.wk[iq] = metadata.weight;
            exx_source_points_.isk[iq] = iq / nq;
        }
        std::ifstream occupations(readin_dir + "eig_occ.txt");
        read_exx_source_occupations(occupations, headers, inp.nspin, exx_source_weights_);
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
    // Preserve the exact stored Cartesian coordinates expected by the existing reader.
    for (int iq = 0; iq < nqs; ++iq)
    {
        exx_source_basis_->kvec_c[iq] = exx_source_points_.kvec_c[iq];
    }
    exx_source_basis_->fft_bundle.initfftmode(inp.fft_mode);
    exx_source_basis_->setuptransform();
    exx_source_basis_->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);
    for (int iq = 0; iq < nqs; ++iq)
    {
        exx_source_points_.ngk[iq] = exx_source_basis_->npwk[iq];
    }
    const int nbands = exx_source_weights_.nc;
    const int nbasis = exx_source_basis_->npwk_max;
    exx_source_psi_.reset(new psi::Psi<T, Device>(nqs, nbands, nbasis, exx_source_points_.ngk, true));
    for (int iq = 0; iq < nqs; ++iq)
    {
        const std::string wfc_filename = ModuleIO::filename_output(readin_dir, "wf", "pw", iq,
            exx_source_points_.ik2iktot, inp.nspin, nqs, 2, false, false, -1);
        try
        {
            validate_source_grid(wfc_filename, *exx_source_basis_, iq);
        }
        catch (const std::exception& error)
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", error.what());
        }
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
