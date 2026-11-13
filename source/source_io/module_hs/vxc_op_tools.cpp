#include "vxc_op_tools.h"

#include "source_base/formatter.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_out/filename.h"
#include "source_io/module_output/provenance_header.h"

#include <fstream>
#include <iomanip>

namespace ModuleIO
{

void set_para2d_MO(const Parallel_Orbitals& pv, const int nbands, Parallel_Blkcyclic& p2d)
{
#ifdef __MPI
    p2d.set(nbands, nbands, pv.para_2d().nb(), pv.para_2d().ictxt());
#else
    p2d.set_serial(nbands, nbands);
#endif
}

template <typename T>
std::vector<T> cVc(const T* V,
                   const T* c,
                   const int nbasis,
                   const int nbands,
                   const Parallel_Orbitals& pv,
                   const Parallel_Blkcyclic& p2d)
{
    std::vector<T> Vc(pv.nloc_wfc, 0.0);
    char transa = 'N';
    char transb = 'N';
    const T alpha = static_cast<T>(1.0);
    const T beta = static_cast<T>(0.0);
#ifdef __MPI
    const int i1 = 1;
    ScalapackConnector::gemm(transa, transb,
        nbasis, nbands, nbasis,
        alpha, V, i1, i1, pv.para_2d().desc(),
        c, i1, i1, pv.desc_wfc,
        beta, Vc.data(), i1, i1, pv.desc_wfc);
#else
    container::BlasConnector::gemm(transa, transb, nbasis, nbands, nbasis, alpha, V, nbasis, c, nbasis, beta, Vc.data(), nbasis);
#endif
    std::vector<T> cVc_result(p2d.size(), 0.0);
    transa = (std::is_same<T, double>::value ? 'T' : 'C');
#ifdef __MPI
    ScalapackConnector::gemm(transa, transb,
        nbands, nbands, nbasis,
        alpha, c, i1, i1, pv.desc_wfc,
        Vc.data(), i1, i1, pv.desc_wfc,
        beta, cVc_result.data(), i1, i1, p2d.desc());
#else
    container::BlasConnector::gemm(transa, transb, nbands, nbands, nbasis, alpha, c, nbasis, Vc.data(), nbasis, beta, cVc_result.data(), nbasis);
#endif
    return cVc_result;
}

template <typename T>
double all_band_energy(const int ik,
                       const std::vector<T>& mat_mo,
                       const Parallel_Blkcyclic& p2d,
                       const ModuleBase::matrix& wg)
{
    double e = 0.0;
    for (int i = 0; i < p2d.nrloc(); ++i)
    {
        for (int j = 0; j < p2d.ncloc(); ++j)
        {
            if (p2d.irglb(i) == p2d.icglb(j))
            {
                e += get_real(mat_mo[j * p2d.nrloc() + i]) * wg(ik, p2d.irglb(i));
            }
        }
    }
    Parallel_Reduce::reduce_all(e);
    return e;
}

template <typename T>
std::vector<double> orbital_energy(const int ik,
                                   const int nbands,
                                   const std::vector<T>& mat_mo,
                                   const Parallel_Blkcyclic& p2d)
{
#ifdef __DEBUG
    assert(nbands >= 0);
#endif
    std::vector<double> e(nbands, 0.0);
    for (int i = 0; i < nbands; ++i)
    {
        if (p2d.in_this_processor(i, i))
        {
            const int index = p2d.icloc(i) * p2d.nrloc() + p2d.irloc(i);
            e[i] = get_real(mat_mo[index]);
        }
    }
    Parallel_Reduce::reduce_all(e.data(), nbands);
    return e;
}

void write_orb_energy(const K_Vectors& kv,
                      const int nspin0,
                      const int nbands,
                      const std::vector<std::vector<double>>& e_orb,
                      const std::string& out_name,
                      const std::string& header,
                      const std::string& global_out_dir,
                      const bool app)
{
    assert(e_orb.size() == kv.get_nks());
    assert(!out_name.empty());
    const int nk = kv.get_nks() / nspin0;
    std::ofstream ofs;
    ofs.open(global_out_dir + out_name,
        app ? std::ios::app : std::ios::out);
    if (!header.empty())
    {
        ofs << header;
    }
    ofs << nk << "\n" << nspin0 << "\n" << nbands << "\n";
    ofs << std::scientific << std::setprecision(16);
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int is = 0; is < nspin0; ++is)
        {
            for (auto e : e_orb[is * nk + ik])
            { // Hartree and eV
                ofs << e / 2. << "\t" << e * ModuleBase::Ry_to_eV << "\n";
            }
        }
    }
}

std::string build_eband_header(const int istep,
                               const double etot_eV,
                               const std::string& term_symbol,
                               const std::string& time_stamp)
{
    const int istep_one_based = istep + 1;
    const std::string step_lines = FmtCore::format(
        "# RELAX STEP %d, Energy: %.8f eV\n"
        "# e_{ikn} = <psi_{ikn} | %s | psi_{ikn}> in the KS orbital basis\n"
        "# units: columns are e/2 in Hartree and e in eV; data rows follow nk, nspin0, nbands\n",
        istep_one_based,
        etot_eV,
        term_symbol.c_str());
    return provenance_header(time_stamp) + step_lines;
}

// Explicit template instantiations
template std::vector<double> cVc<double>(const double*, const double*, const int, const int, const Parallel_Orbitals&, const Parallel_Blkcyclic&);
template std::vector<std::complex<double>> cVc<std::complex<double>>(const std::complex<double>*, const std::complex<double>*, const int, const int, const Parallel_Orbitals&, const Parallel_Blkcyclic&);

template double all_band_energy<double>(const int, const std::vector<double>&, const Parallel_Blkcyclic&, const ModuleBase::matrix&);
template double all_band_energy<std::complex<double>>(const int, const std::vector<std::complex<double>>&, const Parallel_Blkcyclic&, const ModuleBase::matrix&);

template std::vector<double> orbital_energy<double>(const int, const int, const std::vector<double>&, const Parallel_Blkcyclic&);
template std::vector<double> orbital_energy<std::complex<double>>(const int, const int, const std::vector<std::complex<double>>&, const Parallel_Blkcyclic&);

} // namespace ModuleIO
