#include "write_HS_R.h"

#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/LCAO_HS_arrays.hpp"
#include "source_lcao/spar_dh.h"
#include "source_lcao/spar_hsr.h"
#include "source_lcao/spar_st.h"
#include "write_HS_sparse.h"

namespace {
// Helper: Convert sparse map to HContainer
template <typename T>
hamilt::HContainer<T>* sparse_map_to_hcontainer(
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, T>>>& sparse_map,
    const Parallel_Orbitals& pv,
    const int nbasis)
{
    hamilt::HContainer<T>* hc = new hamilt::HContainer<T>(&pv);
    hc->set_zero();

    for (const auto& r_entry : sparse_map)
    {
        const auto& R = r_entry.first;
        for (const auto& row_entry : r_entry.second)
        {
            const size_t row = row_entry.first;
            for (const auto& col_entry : row_entry.second)
            {
                hc->set_value(R.x, R.y, R.z, row, col_entry.first, col_entry.second);
            }
        }
    }

    return hc;
}
} // anonymous namespace

// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix.
// If the absolute value of the matrix element is less than or equal to the
// 'sparse_thr', it will be ignored.

void ModuleIO::output_dSR(const int& istep,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_dSR");
    ModuleBase::timer::start("ModuleIO", "output_dSR");

    sparse_format::cal_dS(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);

    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary, "s");

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_dSR");
    return;
}

void ModuleIO::output_dHR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_dHR");
    ModuleBase::timer::start("ModuleIO", "output_dHR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                         #Print out dH/dR#                          |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    const int nspin = PARAM.inp.nspin;

    if (nspin == 1 || nspin == 4)
    {
        // mohan add 2024-04-01
        const int cspin = 0;

        sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff);
    }
    else if (nspin == 2)
    {
        for (int cspin = 0; cspin < 2; cspin++)
        {
            sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff);
        }
    }
    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary);

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_dHR");
    return;
}

template <typename TK>
void ModuleIO::output_SR(Parallel_Orbitals& pv,
                         const Grid_Driver& grid,
                         hamilt::Hamilt<TK>* p_ham,
                         const std::string& SR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_SR");
    ModuleBase::timer::start("ModuleIO", "output_SR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                 #Print out overlap matrix S(R)#                    |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::cout << " Overlap matrix file is in " << SR_filename << std::endl;
    GlobalV::ofs_running << " Overlap matrix file is in " << SR_filename << std::endl;

    LCAO_HS_Arrays HS_Arrays;

    sparse_format::cal_SR(pv,
                          HS_Arrays.all_R_coor,
                          HS_Arrays.SR_sparse,
                          HS_Arrays.SR_soc_sparse,
                          grid,
                          sparse_thr,
                          p_ham);

    const int istep = 0;

    if (PARAM.inp.nspin == 4)
    {
        ModuleIO::save_sparse(HS_Arrays.SR_soc_sparse,
                              HS_Arrays.all_R_coor,
                              sparse_thr,
                              binary,
                              SR_filename,
                              pv,
                              "S",
                              istep);
    }
    else
    {
        ModuleIO::save_sparse(HS_Arrays.SR_sparse,
                              HS_Arrays.all_R_coor,
                              sparse_thr,
                              binary,
                              SR_filename,
                              pv,
                              "S",
                              istep);
    }

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_SR");
    return;
}

void ModuleIO::output_TR(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         LCAO_HS_Arrays& HS_Arrays,
                         const Grid_Driver& grid,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const std::string& TR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_TR");
    ModuleBase::timer::start("ModuleIO", "output_TR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |           #Print out kinetic energy term matrix T(R)#              |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::stringstream sst;
    if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
    {
        sst << PARAM.globalv.global_matrix_dir << TR_filename << "g" << istep;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }
    else
    {
        sst << PARAM.globalv.global_out_dir << TR_filename;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }

    sparse_format::cal_TR(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);

    ModuleIO::save_sparse(HS_Arrays.TR_sparse,
                          HS_Arrays.all_R_coor,
                          sparse_thr,
                          binary,
                          sst.str().c_str(),
                          pv,
                          "T",
                          istep);

    sparse_format::destroy_T_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_TR");
    return;
}

template void ModuleIO::output_SR<double>(Parallel_Orbitals& pv,
                                          const Grid_Driver& grid,
                                          hamilt::Hamilt<double>* p_ham,
                                          const std::string& SR_filename,
                                          const bool& binary,
                                          const double& sparse_thr);
template void ModuleIO::output_SR<std::complex<double>>(Parallel_Orbitals& pv,
                                                        const Grid_Driver& grid,
                                                        hamilt::Hamilt<std::complex<double>>* p_ham,
                                                        const std::string& SR_filename,
                                                        const bool& binary,
                                                        const double& sparse_thr);

#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_io/module_output/ucell_io.h"

std::string ModuleIO::hsr_gen_fname(const std::string& prefix,
                                     const int ispin,
                                     const bool append,
                                     const int istep)
{
    if (!append && istep >= 0)
    {
        return prefix + std::to_string(ispin + 1) + "g" + std::to_string(istep + 1) + "_nao.csr";
    }
    else
    {
        return prefix + std::to_string(ispin + 1) + "_nao.csr";
    }
}

std::string ModuleIO::dhr_gen_fname(const std::string& prefix,
                                     const int ispin,
                                     const bool append,
                                     const int istep)
{
    std::string fname = prefix + "rs" + std::to_string(ispin + 1);
    if (!append && istep >= 0)
    {
        fname += "g" + std::to_string(istep + 1);
    }
    fname += "_nao.csr";
    return fname;
}

template <typename TR>
void ModuleIO::write_hcontainer_csr(const std::string& fname,
                                     const UnitCell* ucell,
                                     const int precision,
                                     hamilt::HContainer<TR>* mat_serial,
                                     const int istep,
                                     const int ispin,
                                     const int nspin,
                                     const std::string& label)
{
    std::ofstream ofs;
    if (istep <= 0)
    {
        ofs.open(fname);
    }
    else
    {
        ofs.open(fname, std::ios::app);
    }

    ofs << " --- Ionic Step " << istep + 1 << " ---" << std::endl;
    ofs << " # print " << label << " matrix in real space " << label << "(R)" << std::endl;
    ofs << " " << nspin << " # number of spin directions" << std::endl;
    ofs << " " << ispin + 1 << " # spin index" << std::endl;
    ofs << " " << mat_serial->get_nbasis() << " # number of localized basis" << std::endl;
    ofs << " " << mat_serial->size_R_loop() << " # number of Bravais lattice vector R" << std::endl;
    ofs << std::endl;

    ModuleIO::UcellIO::write_ucell(ofs, ucell);
    ofs << std::endl;

    const double sparse_threshold = 1e-10;
    hamilt::Output_HContainer<TR> out(mat_serial, ofs, sparse_threshold, precision);
    out.write();
    ofs.close();
}

template <typename TR>
void ModuleIO::write_hsr(const std::vector<hamilt::HContainer<TR>*>& hr_vec,
                          const hamilt::HContainer<TR>* sr,
                          const UnitCell* ucell,
                          const int precision,
                          const Parallel_2D& paraV,
                          const bool append,
                          const int* iat2iwt,
                          const int nat,
                          const int istep)
{
    const int nspin = hr_vec.size();
    assert(nspin > 0);

    // Output HR (one file per spin)
    for (int ispin = 0; ispin < nspin; ispin++)
    {
        const int nbasis = hr_vec[ispin]->get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, paraV.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        hamilt::HContainer<TR> hr_serial(&serialV);
        hamilt::gatherParallels(*hr_vec[ispin], &hr_serial, 0);
#else
        hamilt::HContainer<TR> hr_serial(*hr_vec[ispin]);
#endif

        if (GlobalV::MY_RANK == 0)
        {
            std::string fname = PARAM.globalv.global_out_dir
                                + hsr_gen_fname("hrs", ispin, append, istep);
            write_hcontainer_csr(fname, ucell, precision, &hr_serial, istep, ispin, nspin, "H");
        }
    }

    // Output SR (single file)
    {
        const int nbasis = sr->get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, paraV.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        hamilt::HContainer<TR> sr_serial(&serialV);
        hamilt::gatherParallels(*sr, &sr_serial, 0);
#else
        hamilt::HContainer<TR> sr_serial(*sr);
#endif

        if (GlobalV::MY_RANK == 0)
        {
            std::string fname = PARAM.globalv.global_out_dir
                                + hsr_gen_fname("srs", 0, append, istep);
            write_hcontainer_csr(fname, ucell, precision, &sr_serial, istep, 0, 1, "S");
        }
    }
}

// Explicit instantiations
template void ModuleIO::write_hcontainer_csr<double>(
    const std::string&, const UnitCell*, const int,
    hamilt::HContainer<double>*, const int, const int, const int, const std::string&);
template void ModuleIO::write_hcontainer_csr<std::complex<double>>(
    const std::string&, const UnitCell*, const int,
    hamilt::HContainer<std::complex<double>>*, const int, const int, const int, const std::string&);

template void ModuleIO::write_hsr<double>(
    const std::vector<hamilt::HContainer<double>*>&,
    const hamilt::HContainer<double>*,
    const UnitCell*, const int, const Parallel_2D&,
    const bool, const int*, const int, const int);
template void ModuleIO::write_hsr<std::complex<double>>(
    const std::vector<hamilt::HContainer<std::complex<double>>*>&,
    const hamilt::HContainer<std::complex<double>>*,
    const UnitCell*, const int, const Parallel_2D&,
    const bool, const int*, const int, const int);


template <typename TR>
void ModuleIO::write_matrix_r(const std::string& matrix_label,
                               const std::string& description,
                               const std::vector<hamilt::HContainer<TR>*>& matrices,
                               const UnitCell* ucell,
                               const int precision,
                               const Parallel_2D& paraV,
                               const bool append,
                               const int* iat2iwt,
                               const int nat,
                               const int istep)
{
    const int nspin = matrices.size();
    assert(nspin > 0);
    
    for (int ispin = 0; ispin < nspin; ispin++)
    {
        const int nbasis = matrices[ispin]->get_nbasis();
        
        // Generate filename
        std::string fname = dhr_gen_fname(matrix_label, ispin, append, istep);
        if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
        {
            fname = PARAM.globalv.global_matrix_dir + fname;
        }
        else
        {
            fname = PARAM.globalv.global_out_dir + fname;
        }
        
        // Gather parallel matrix to serial
#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, paraV.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        
        hamilt::HContainer<TR> matrix_serial(&serialV);
        hamilt::gatherParallels(*matrices[ispin], &matrix_serial, 0);
        
        if (GlobalV::MY_RANK == 0)
        {
            write_hcontainer_csr(fname, ucell, precision, &matrix_serial, istep, ispin, nspin, description);
        }
#else
        write_hcontainer_csr(fname, ucell, precision, matrices[ispin], istep, ispin, nspin, description);
#endif
    }
}

// Template instantiations
template void ModuleIO::write_matrix_r<double>(
    const std::string&,
    const std::string&,
    const std::vector<hamilt::HContainer<double>*>&,
    const UnitCell*,
    const int,
    const Parallel_2D&,
    const bool,
    const int*,
    const int,
    const int);

template void ModuleIO::write_matrix_r<std::complex<double>>(
    const std::string&,
    const std::string&,
    const std::vector<hamilt::HContainer<std::complex<double>>*>&,
    const UnitCell*,
    const int,
    const Parallel_2D&,
    const bool,
    const int*,
    const int,
    const int);
