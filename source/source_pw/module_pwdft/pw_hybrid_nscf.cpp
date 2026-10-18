#include "pw_hybrid_nscf.h"
#include "pw_hybrid_source.h"
#include "source_io/module_wf/read_wfc_pw.h"
#include "source_base/module_out/filename.h"
#include "source_hamilt/module_xc/general_exx_info.h"
#include "source_pw/module_pwdft/op_pw_exx.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace pw
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

bool source_parameter_matches(const std::string& key, std::string source_value, std::string target_value)
{
    if (key == "dft_functional")
    {
        const auto lowercase = [](unsigned char value) { return std::tolower(value); };
        std::transform(source_value.begin(), source_value.end(), source_value.begin(), lowercase);
        std::transform(target_value.begin(), target_value.end(), target_value.begin(), lowercase);
    }
    if (source_value == target_value)
    {
        return true;
    }
    std::istringstream source(source_value);
    std::istringstream target(target_value);
    double source_number = 0;
    double target_number = 0;
    if (!(source >> source_number) || !(target >> target_number))
    {
        return false;
    }
    std::string remaining;
    if ((source >> remaining) || (target >> remaining))
    {
        return false;
    }
    // INPUT.info may print doubles with the standard six significant digits.
    const double scale = std::max(std::abs(source_number), std::abs(target_number));
    return std::isfinite(source_number) && std::isfinite(target_number)
           && std::abs(source_number - target_number) <= 5e-6 * scale;
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
            // Zero-valued Fock coefficients are equivalent to omitting the term.
            if (term.first == "exx_fock_alpha" && std::stod(number) == 0.0)
            {
                continue;
            }
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
            ModuleBase::WARNING("PW_Hybrid_NSCF", "EXX source configuration is incomplete or unavailable in INPUT.info; "
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
            // Mirror the target-side normalization for saved zero Fock coefficients.
            if (item.first == "exx_fock_alpha" && std::stod(source_value) == 0.0)
            {
                continue;
            }
            if (!(current >> target_value))
            {
                mismatch = true;
                break;
            }
            if (!source_parameter_matches(item.first, source_value, target_value))
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
            ModuleBase::WARNING("PW_Hybrid_NSCF", message);
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
    bool has_fock = false;
    const auto fock = info.coulomb_param.find(Conv_Coulomb_Pot_K::Coulomb_Type::Fock);
    if (fock != info.coulomb_param.end())
    {
        // An explicit zero coefficient (`exx_fock_alpha 0`) retains a Fock entry but is
        // equivalent to omitting the term; only a nonzero coefficient is unscreened exchange.
        for (const auto& term : fock->second)
        {
            if (std::stod(term.at("alpha")) != 0.0)
            {
                has_fock = true;
                break;
            }
        }
    }
    bool compatible_device = inp.device == "cpu";
#ifdef __CUDA
    compatible_device = compatible_device || inp.device == "gpu";
#endif
    const bool compatible_execution = compatible_device && inp.kpar == 1 && inp.bndpar == 1 && inp.nspin != 4;
    const bool compatible_exchange = !inp.exxace && !inp.exx_gamma_extra && !has_fock;
    const bool compatible_targets = inp.symmetry != "1" && inp.init_wfc != "file" && !inp.cal_force && !inp.cal_stress;
    if (!compatible_execution || !compatible_exchange || !compatible_targets)
    {
        ModuleBase::WARNING_QUIT("PW_Hybrid_NSCF",
            "Hybrid NSCF currently requires CPU or CUDA GPU, kpar=bndpar=1, nspin=1/2, symmetry=-1/0, "
            "exxace=false, exx_gamma_extra=false, screened exchange, "
            "cal_force=cal_stress=false and fresh target wavefunctions");
    }
}

template <typename T, typename Device>
void prepare_exx_nscf(const Input_para& inp,
                      const UnitCell& ucell,
                      const ModulePW::PW_Basis_K* pw_wfc,
                      hamilt::OperatorEXXPW<T, Device>& op_exx,
                      const std::string& readin_dir)
{
    ModuleBase::timer::start("PW_Hybrid_NSCF", "prepare_exx_nscf");
    std::unique_ptr<K_Vectors> source_points(new K_Vectors);
    ModuleBase::matrix source_weights;
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
        source_points->set_nks(nks);
        source_points->set_nkstot(nks);
        const int nq = nks / spin_mult;
        source_points->set_nkstot_nospin(nq);
        source_points->set_spin_mult(spin_mult);
        source_points->kvec_d.resize(nks);
        source_points->kvec_c.resize(nks);
        source_points->wk.resize(nks);
        source_points->isk.resize(nks);
        source_points->ngk.resize(nks);
        source_points->ik2iktot.resize(nks);
        for (int iq = 0; iq < nks; ++iq)
        {
            source_points->ik2iktot[iq] = iq;
        }
        const ModuleBase::Matrix3 inverse_reciprocal = ucell.G.Inverse();
        for (int iq = 0; iq < nks; ++iq)
        {
            const std::string filename = ModuleIO::filename_output(readin_dir, "wf", "pw", iq,
                source_points->ik2iktot, inp.nspin, nks, 2, false, false, -1);
            std::ifstream in(filename, std::ios::binary);
            const auto metadata = read_exx_source_header(in);
            if (metadata.ecutwfc != inp.ecutwfc)
            {
                throw std::runtime_error("EXX source wavefunction cutoff differs from the current basis");
            }
            headers.push_back(metadata);
            auto& q = source_points->kvec_c[iq];
            q.x = metadata.kvec_c[0];
            q.y = metadata.kvec_c[1];
            q.z = metadata.kvec_c[2];
            source_points->kvec_d[iq] = q * inverse_reciprocal;
            source_points->wk[iq] = metadata.weight;
            source_points->isk[iq] = iq / nq;
        }
        std::ifstream occupations(readin_dir + "eig_occ.txt");
        read_exx_source_occupations(occupations, headers, inp.nspin, source_weights);
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("PW_Hybrid_NSCF", error.what());
    }
    const int nqs = source_points->get_nks();
    int nq = source_points->get_nkstot_nospin();
    source_points->para_k.kinfo(nq, 1, 0, pw_wfc->poolrank, pw_wfc->poolnproc, inp.nspin);
    std::unique_ptr<ModulePW::PW_Basis_K> source_basis(
        new ModulePW::PW_Basis_K(inp.device, inp.precision));
#ifdef __MPI
    source_basis->initmpi(pw_wfc->poolnproc, pw_wfc->poolrank, pw_wfc->pool_world);
#endif
    source_basis->initgrids(ucell.lat0, ucell.latvec,
                               pw_wfc->nx, pw_wfc->ny, pw_wfc->nz);
    source_basis->initparameters(false, inp.ecutwfc, nqs, source_points->kvec_d.data());
    // Preserve the exact stored Cartesian coordinates expected by the existing reader.
    for (int iq = 0; iq < nqs; ++iq)
    {
        source_basis->kvec_c[iq] = source_points->kvec_c[iq];
    }
    source_basis->fft_bundle.initfftmode(inp.fft_mode);
    source_basis->setuptransform();
    source_basis->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);
    for (int iq = 0; iq < nqs; ++iq)
    {
        source_points->ngk[iq] = source_basis->npwk[iq];
    }
    const int nbands = source_weights.nc;
    const int nbasis = source_basis->npwk_max;
    std::unique_ptr<psi::Psi<T, Device>> source_psi(
        new psi::Psi<T, Device>(nqs, nbands, nbasis, source_points->ngk, true));
    for (int iq = 0; iq < nqs; ++iq)
    {
        const std::string wfc_filename = ModuleIO::filename_output(readin_dir, "wf", "pw", iq,
            source_points->ik2iktot, inp.nspin, nqs, 2, false, false, -1);
        try
        {
            validate_source_grid(wfc_filename, *source_basis, iq);
        }
        catch (const std::exception& error)
        {
            ModuleBase::WARNING_QUIT("PW_Hybrid_NSCF", error.what());
        }
        ModuleBase::ComplexMatrix wfc(nbands, nbasis);
        ModuleIO::read_wfc_pw(wfc_filename, source_basis.get(), pw_wfc->poolrank,
                             pw_wfc->poolnproc, nbands, 1, iq, iq, nqs, wfc);
        source_psi->fix_k(iq);
        std::vector<T> coefficients(nbands * nbasis);
        for (int ib = 0; ib < nbands; ++ib)
        {
            for (int ig = 0; ig < nbasis; ++ig)
            {
                coefficients[ib * nbasis + ig] = static_cast<T>(wfc(ib, ig));
            }
        }
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(
            source_psi->get_pointer(), coefficients.data(), coefficients.size());
    }
    op_exx.set_source(std::move(source_basis), std::move(source_points),
        *source_psi, std::move(source_weights));
    ModuleBase::timer::end("PW_Hybrid_NSCF", "prepare_exx_nscf");
}


template void prepare_exx_nscf<std::complex<float>, base_device::DEVICE_CPU>(
    const Input_para&, const UnitCell&, const ModulePW::PW_Basis_K*,
    hamilt::OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>&, const std::string&);
template void prepare_exx_nscf<std::complex<double>, base_device::DEVICE_CPU>(
    const Input_para&, const UnitCell&, const ModulePW::PW_Basis_K*,
    hamilt::OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>&, const std::string&);
#if defined(__CUDA) || defined(__ROCM)
template void prepare_exx_nscf<std::complex<float>, base_device::DEVICE_GPU>(
    const Input_para&, const UnitCell&, const ModulePW::PW_Basis_K*,
    hamilt::OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>&, const std::string&);
template void prepare_exx_nscf<std::complex<double>, base_device::DEVICE_GPU>(
    const Input_para&, const UnitCell&, const ModulePW::PW_Basis_K*,
    hamilt::OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>&, const std::string&);
#endif
} // namespace pw
