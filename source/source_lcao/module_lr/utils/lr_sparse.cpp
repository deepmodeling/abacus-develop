#include "lr_util_hcontainer.h"
#include "source_io/module_hs/hs_sparse_io.h"
#include "source_lcao/spar_hsr.h"
#include <ostream>

namespace LR_Util
{
template<typename TR>
void save_HR(const hamilt::HContainer<TR>& hR, const std::string& filename,
             const Parallel_Orbitals& pv, const std::string& out_dir,
             const double sparse_thr, std::ostream& running_log)
{
    const auto all_R = sparse_format::get_R_range(hR);
    ModuleIO::SparseRMatrix<TR> sparse_matrix;
    sparse_format::cal_HContainer(pv, sparse_thr, hR, sparse_matrix);
    ModuleIO::SparseWriteOptions options;
    options.filename = out_dir + filename;
    options.label = "EDM";
    options.threshold = sparse_thr;
    options.istep = 0;
    ModuleIO::save_sparse(sparse_matrix, all_R, pv, options);
    running_log << "Write file: " << options.filename << std::endl;
}

template<typename TK, typename TR>
void save_DMR(const module_dm::DensityMatrix<TK, TR>& dm, const std::string& filename,
              const Parallel_Orbitals& pv, const std::string& out_dir,
              const double sparse_thr, std::ostream& running_log)
{
    const auto& blocks = dm.get_dmr_vec();
    for (std::size_t spin = 0; spin < blocks.size(); ++spin)
    {
        const std::string spin_filename = filename + "_s" + std::to_string(spin);
        save_HR(*blocks[spin], spin_filename, pv, out_dir, sparse_thr, running_log);
    }
}

using Complex = std::complex<double>;
template void save_HR(const hamilt::HContainer<double>&, const std::string&,
                      const Parallel_Orbitals&, const std::string&, double, std::ostream&);
template void save_HR(const hamilt::HContainer<Complex>&, const std::string&,
                      const Parallel_Orbitals&, const std::string&, double, std::ostream&);
template void save_DMR(const module_dm::DensityMatrix<double, double>&, const std::string&,
                       const Parallel_Orbitals&, const std::string&, double, std::ostream&);
template void save_DMR(const module_dm::DensityMatrix<Complex, double>&, const std::string&,
                       const Parallel_Orbitals&, const std::string&, double, std::ostream&);
template void save_DMR(const module_dm::DensityMatrix<Complex, Complex>&, const std::string&,
                       const Parallel_Orbitals&, const std::string&, double, std::ostream&);
}
