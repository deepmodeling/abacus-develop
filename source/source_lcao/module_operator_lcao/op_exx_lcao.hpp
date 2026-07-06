#ifndef OPEXXLCAO_HPP
#define OPEXXLCAO_HPP
#ifdef __EXX

#include "op_exx_lcao.h"
#include "source_base/parallel_reduce.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_restart/restart.h"
#include "source_io/module_restart/restart_exx_csr.h"
#include "source_lcao/module_hcontainer/read_hcontainer.h"
#include "source_lcao/module_ri/Exx_LRI_interface.h"
#include "source_lcao/module_ri/RI_2D_Comm.h"
#include "source_lcao/module_rt/td_info.h"

namespace hamilt
{
using TAC = std::pair<int, std::array<int, 3>>;

// allocate according to the read-in HexxR, used in nscf
template <typename Tdata, typename TR>
void reallocate_hcontainer(const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
                           HContainer<TR>* hR,
                           const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest)
{
    auto* pv = hR->get_paraV();
    bool need_allocate = false;
    for (auto& Htmp1: Hexxs[0])
    {
        const int& iat0 = Htmp1.first;
        for (auto& Htmp2: Htmp1.second)
        {
            const int& iat1 = Htmp2.first.first;
            if (pv->get_nrow_atom(iat0) > 0 && pv->get_ncol_atom(iat1) > 0)
            {
                const Abfs::Vector3_Order<int>& R = RI_Util::array3_to_Vector3(
                    (cell_nearest ? cell_nearest->get_cell_nearest_discrete(iat0, iat1, Htmp2.first.second)
                                  : Htmp2.first.second));
                BaseMatrix<TR>* HlocR = hR->find_matrix(iat0, iat1, R.x, R.y, R.z);
                if (HlocR == nullptr)
                { // add R to HContainer
                    need_allocate = true;
                    AtomPair<TR> tmp(iat0, iat1, R.x, R.y, R.z, pv);
                    hR->insert_pair(tmp);
                }
            }
        }
    }
    if (need_allocate)
    {
        hR->allocate(nullptr, true);
    }
}

/// allocate according to BvK cells, used in scf
template <typename TR>
void reallocate_hcontainer(const int nat,
                           HContainer<TR>* hR,
                           const std::array<int, 3>& Rs_period,
                           const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest)
{
    auto* pv = hR->get_paraV();
    auto Rs = RI_Util::get_Born_von_Karmen_cells(Rs_period);
    bool need_allocate = false;
    for (int iat0 = 0; iat0 < nat; ++iat0)
    {
        for (int iat1 = 0; iat1 < nat; ++iat1)
        {
            // complete the atom pairs that has orbitals in this processor but not in hR due to the adj_list
            // but adj_list is not enought for EXX, which is more nonlocal than Nonlocal
            if (pv->get_nrow_atom(iat0) > 0 && pv->get_ncol_atom(iat1) > 0)
            {
                for (auto& cell: Rs)
                {
                    const Abfs::Vector3_Order<int>& R = RI_Util::array3_to_Vector3(
                        (cell_nearest ? cell_nearest->get_cell_nearest_discrete(iat0, iat1, cell) : cell));
                    BaseMatrix<TR>* HlocR = hR->find_matrix(iat0, iat1, R.x, R.y, R.z);

                    if (HlocR == nullptr)
                    { // add R to HContainer
                        need_allocate = true;
                        AtomPair<TR> tmp(iat0, iat1, R.x, R.y, R.z, pv);
                        hR->insert_pair(tmp);
                    }
                }
            }
        }
    }
    if (need_allocate)
    {
        hR->allocate(nullptr, true);
    }
}

template <typename TK, typename TR>
OperatorEXX<OperatorLCAO<TK, TR>>::OperatorEXX(
    HS_Matrix_K<TK>* hsk_in,
    hamilt::HContainer<TR>* hR_in,
    const UnitCell& ucell,
    const K_Vectors& kv_in,
    std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in,
    std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in,
    Add_Hexx_Type add_hexx_type_in)
    : OperatorLCAO<TK, TR>(hsk_in, kv_in.kvec_d, hR_in), ucell(ucell), kv(kv_in), Hexxd(Hexxd_in), Hexxc(Hexxc_in),
      add_hexx_type(add_hexx_type_in)
{
    this->cal_type = calculation_type::lcao_exx;
    // This one-shot constructor never builds cell_nearest, so cal_dH() must not use it:
    // the (d)Hexxs from LibRI are in native cells, and the dH output mirrors the H-term
    // writer (write_h_exx_impl), which also passes a nullptr cell_nearest.
    this->use_cell_nearest = false;
}

template <typename TK, typename TR>
OperatorEXX<OperatorLCAO<TK, TR>>::OperatorEXX(HS_Matrix_K<TK>* hsk_in,
                                               HContainer<TR>* hR_in,
                                               const UnitCell& ucell_in,
                                               const K_Vectors& kv_in,
                                               Exx_LRI_Interface<TK, double>* exd_in,
                                               Exx_LRI_Interface<TK, std::complex<double>>* exc_in,
                                               Add_Hexx_Type add_hexx_type_in,
                                               const int istep_in,
                                               const bool restart_in)
    : OperatorEXX<OperatorLCAO<TK, TR>>(hsk_in,
                                        hR_in,
                                        ucell_in,
                                        kv_in,
                                        exd_in ? &exd_in->get_Hexxs() : nullptr,
                                        exc_in ? &exc_in->get_Hexxs() : nullptr,
                                        add_hexx_type_in)
{
    this->exd = exd_in;
    this->exc = exc_in;
    const_cast<int&>(this->istep) = istep_in;
    this->restart = restart_in;
    ModuleBase::TITLE("OperatorEXX", "OperatorEXX");
    const Parallel_Orbitals* const pv = hR_in->get_paraV();

    if (PARAM.inp.calculation == "nscf" && GlobalC::exx_info.info_global.cal_exx)
    { // for nscf, calculate HexxR from the read-in DM, or read HexxR in
        auto file_name_list_csr = []() -> std::vector<std::string> {
            std::vector<std::string> file_name_list;
            for (int irank = 0; irank < PARAM.globalv.nproc; ++irank)
            {
                for (int is = 0; is < PARAM.inp.nspin; ++is)
                {
                    file_name_list.push_back(PARAM.globalv.global_readin_dir + "HexxR" + std::to_string(irank) + "_"
                                             + std::to_string(is) + ".csr");
                }
            }
            return file_name_list;
        };
        auto file_name_list_cereal = []() -> std::vector<std::string> {
            std::vector<std::string> file_name_list;
            for (int irank = 0; irank < PARAM.globalv.nproc; ++irank)
            {
                file_name_list.push_back("HexxR_" + std::to_string(irank));
            }
            return file_name_list;
        };
        auto check_exist = [](const std::vector<std::string>& file_name_list) -> bool {
            for (const std::string& file_name: file_name_list)
            {
                std::ifstream ifs(file_name);
                if (!ifs.is_open())
                {
                    return false;
                }
            }
            return true;
        };

        if (PARAM.inp.init_chg == "dm" || PARAM.inp.init_chg == "dm_no_renormalize")
        {
            // 1. cal Cs, Vs
            if (GlobalC::exx_info.info_ri.real_number)
            {
                this->exd->cal_exx_ions(ucell, PARAM.inp.out_ri_cv);
            }
            else
            {
                this->exc->cal_exx_ions(ucell, PARAM.inp.out_ri_cv);
            }

            // 2. read DM
            const int nspin_dm = (PARAM.inp.nspin == 2) ? 2 : 1;
            std::vector<hamilt::HContainer<double>*> dmR_vec(nspin_dm);
            for (int is = 0; is < nspin_dm; ++is)
            {
                const std::string dmfile
                    = PARAM.globalv.global_readin_dir + "/dmrs" + std::to_string(is + 1) + "_nao.csr";
                dmR_vec[is] = new hamilt::HContainer<double>(const_cast<Parallel_Orbitals*>(pv));
                hamilt::Read_HContainer<double> reader_dm(dmR_vec[is], dmfile, PARAM.globalv.nlocal, &ucell);
                reader_dm.read();
            }

            // 3. DM->Ds->Hexx (do not use symmetry for nscf)
            XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
            if (GlobalC::exx_info.info_ri.real_number)
            {
                const auto& Ds = RI_2D_Comm::dm_container_to_Ds<double, double>(dmR_vec, ucell, *pv, PARAM.inp.nspin);
                this->exd->cal_exx_elec(Ds, ucell, *pv);
            }
            else
            {
                const auto& Ds = RI_2D_Comm::dm_container_to_Ds<double, std::complex<double>>(dmR_vec,
                                                                                              ucell,
                                                                                              *pv,
                                                                                              PARAM.inp.nspin);
                this->exc->cal_exx_elec(Ds, ucell, *pv);
            }
        }
        else // need to read HexxR
        {
            std::cout << " Attention: The number of MPI processes must be strictly identical between SCF and NSCF when "
                         "computing exact-exchange."
                      << std::endl;
            if (check_exist(file_name_list_csr()))
            {
                // read HexxR first and reallocate hR according to the read-in HexxR
                const std::string file_name_exx_csr
                    = PARAM.globalv.global_readin_dir + "HexxR" + std::to_string(PARAM.globalv.myrank);
                // Read HexxR in CSR format
                if (GlobalC::exx_info.info_ri.real_number)
                {
                    ModuleIO::read_Hexxs_csr(file_name_exx_csr, ucell, PARAM.inp.nspin, PARAM.globalv.nlocal, *Hexxd);
                }
                else
                {
                    ModuleIO::read_Hexxs_csr(file_name_exx_csr, ucell, PARAM.inp.nspin, PARAM.globalv.nlocal, *Hexxc);
                }
            }
            else if (check_exist(file_name_list_cereal()))
            {
                // Read HexxR in binary format (old version)
                const std::string file_name_exx_cereal
                    = PARAM.globalv.global_readin_dir + "HexxR_" + std::to_string(PARAM.globalv.myrank);
                std::ifstream ifs(file_name_exx_cereal, std::ios::binary);
                if (!ifs)
                {
                    ModuleBase::WARNING_QUIT("OperatorEXX", "Can't open EXX file < " + file_name_exx_cereal + " >.");
                }
                if (GlobalC::exx_info.info_ri.real_number)
                {
                    ModuleIO::read_Hexxs_cereal(file_name_exx_cereal, *Hexxd);
                }
                else
                {
                    ModuleIO::read_Hexxs_cereal(file_name_exx_cereal, *Hexxc);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("OperatorEXX", "Can't open EXX file in " + PARAM.globalv.global_readin_dir);
            }
        }
        // reallocate hR according to Hexx(R)
        if (this->add_hexx_type == Add_Hexx_Type::R)
        {
            if (GlobalC::exx_info.info_ri.real_number)
            {
                reallocate_hcontainer(*this->Hexxd, this->hR);
            }
            else
            {
                reallocate_hcontainer(*this->Hexxc, this->hR);
            }
        }
        this->use_cell_nearest = false;
    }
    else
    { // if scf and Add_Hexx_Type::R, init cell_nearest and reallocate hR according to BvK cells
        if (this->add_hexx_type == Add_Hexx_Type::R)
        {
            // if k points has no shift, use cell_nearest to reduce the memory cost
            this->use_cell_nearest = (ModuleBase::Vector3<double>(std::fmod(this->kv.get_koffset(0), 1.0),
                                                                  std::fmod(this->kv.get_koffset(1), 1.0),
                                                                  std::fmod(this->kv.get_koffset(2), 1.0))
                                          .norm()
                                      < 1e-10);

            const std::array<int, 3> Rs_period = {this->kv.nmp[0], this->kv.nmp[1], this->kv.nmp[2]};
            if (this->use_cell_nearest)
            {
                this->cell_nearest = init_cell_nearest(ucell, Rs_period);
                reallocate_hcontainer(ucell.nat, this->hR, Rs_period, &this->cell_nearest);
            }
            else
            {
                reallocate_hcontainer(ucell.nat, this->hR, Rs_period);
            }
        }

        if (this->restart)
        { ///  Now only Hexx depends on DM, so we can directly read Hexx to reduce the computational cost.
            /// If other operators depends on DM, we can also read DM and then calculate the operators to save the
            /// memory to store operator terms.

            if (this->add_hexx_type == Add_Hexx_Type::k)
            {
                /// read in Hexx(k)
                if (std::is_same<TK, double>::value)
                {
                    this->Hexxd_k_load.resize(this->kv.get_nks());
                    for (int ik = 0; ik < this->kv.get_nks(); ik++)
                    {
                        this->Hexxd_k_load[ik].resize(pv->get_local_size(), 0.0);
                        this->restart = GlobalC::restart.load_disk("Hexx",
                                                                   ik,
                                                                   pv->get_local_size(),
                                                                   this->Hexxd_k_load[ik].data(),
                                                                   false);
                        if (!this->restart)
                        {
                            break;
                        }
                    }
                }
                else
                {
                    this->Hexxc_k_load.resize(this->kv.get_nks());
                    for (int ik = 0; ik < this->kv.get_nks(); ik++)
                    {
                        this->Hexxc_k_load[ik].resize(pv->get_local_size(), 0.0);
                        this->restart = GlobalC::restart.load_disk("Hexx",
                                                                   ik,
                                                                   pv->get_local_size(),
                                                                   this->Hexxc_k_load[ik].data(),
                                                                   false);
                        if (!this->restart)
                        {
                            break;
                        }
                    }
                }
            }
            else if (this->add_hexx_type == Add_Hexx_Type::R)
            {
                // read in Hexx(R)
                const std::string restart_HR_path
                    = GlobalC::restart.folder + "HexxR" + std::to_string(PARAM.globalv.myrank);
                int all_exist = 1;
                for (int is = 0; is < PARAM.inp.nspin; ++is)
                {
                    std::ifstream ifs(restart_HR_path + "_" + std::to_string(is) + ".csr");
                    if (!ifs)
                    {
                        all_exist = 0;
                        break;
                    }
                }
                // Add MPI communication to synchronize all_exist across processes
#ifdef __MPI
                Parallel_Reduce::reduce_min(all_exist);
#endif
                if (all_exist)
                {
                    // Read HexxR in CSR format
                    if (GlobalC::exx_info.info_ri.real_number)
                    {
                        ModuleIO::read_Hexxs_csr(restart_HR_path, ucell, PARAM.inp.nspin, PARAM.globalv.nlocal, *Hexxd);
                    }
                    else
                    {
                        ModuleIO::read_Hexxs_csr(restart_HR_path, ucell, PARAM.inp.nspin, PARAM.globalv.nlocal, *Hexxc);
                    }
                }
                else
                {
                    // Read HexxR in binary format (old version)
                    const std::string restart_HR_path_cereal
                        = GlobalC::restart.folder + "HexxR_" + std::to_string(PARAM.globalv.myrank);
                    std::ifstream ifs(restart_HR_path_cereal, std::ios::binary);
                    int all_exist_cereal = ifs ? 1 : 0;
#ifdef __MPI
                    Parallel_Reduce::reduce_min(all_exist_cereal);
#endif
                    if (!all_exist_cereal)
                    {
                        // no HexxR file in CSR or binary format
                        this->restart = false;
                    }
                    else
                    {
                        if (GlobalC::exx_info.info_ri.real_number)
                        {
                            ModuleIO::read_Hexxs_cereal(restart_HR_path_cereal, *Hexxd);
                        }
                        else
                        {
                            ModuleIO::read_Hexxs_cereal(restart_HR_path_cereal, *Hexxc);
                        }
                    }
                }
            }

            if (!this->restart)
            {
                std::cout << "WARNING: Hexx not found, restart from the non-exx loop." << std::endl
                          << "If the loaded charge density is EXX-solved, this may lead to poor convergence."
                          << std::endl;
            }
            GlobalC::restart.info_load.load_H_finish = this->restart;
        }
    }
}
template <typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("OperatorEXX", "contributeHR");
    // Peize Lin add 2016-12-03

    // 1. For NSCF
    if (PARAM.inp.calculation == "nscf")
    {
        // Do nothing here, allow the code to proceed and calculate EXX.
    }
    // 2. For the first ionic step of SCF, relaxation, or MD:
    else if (this->istep == 0)
    {
        const int two_level_step
            = GlobalC::exx_info.info_ri.real_number ? this->exd->get_two_level_step() : this->exc->get_two_level_step();

        // Check if we are in the pre-convergence stage of the two-level SCF (i.e., the pure GGA loop)
        bool in_gga_pre_loop = (two_level_step == 0);

        // Check if a high-quality initial guess is missing (neither reading wavefunctions from a file nor restarting)
        bool lacks_good_guess = (PARAM.inp.init_wfc != "file" && !this->restart);

        // If in the pre-convergence loop and lacking a good initial guess, skip adding the EXX contribution
        if (in_gga_pre_loop && lacks_good_guess)
        {
            return; // In the non-EXX loop, skip adding EXX contribution
        }
    }
    // 3. For subsequent ionic steps (istep > 0), add EXX normally

    if (this->add_hexx_type == Add_Hexx_Type::k)
    {
        return;
    }

    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        // add H(R) normally
        if (GlobalC::exx_info.info_ri.real_number)
        {
            RI_2D_Comm::add_HexxR(this->current_spin,
                                  GlobalC::exx_info.info_global.hybrid_alpha,
                                  *this->Hexxd,
                                  *this->hR->get_paraV(),
                                  PARAM.globalv.npol,
                                  *this->hR,
                                  this->use_cell_nearest ? &this->cell_nearest : nullptr);
        }
        else
        {
            RI_2D_Comm::add_HexxR(this->current_spin,
                                  GlobalC::exx_info.info_global.hybrid_alpha,
                                  *this->Hexxc,
                                  *this->hR->get_paraV(),
                                  PARAM.globalv.npol,
                                  *this->hR,
                                  this->use_cell_nearest ? &this->cell_nearest : nullptr);
        }
    }
    if (PARAM.inp.nspin == 2)
    {
        this->current_spin = 1 - this->current_spin;
    }
}

template <typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorEXX", "constributeHk");
    const bool has_workflow = GlobalC::exx_info.info_ri.real_number ? (this->exd != nullptr) : (this->exc != nullptr);
    int two_level_step = 0;
    if (has_workflow)
    {
        two_level_step
            = GlobalC::exx_info.info_ri.real_number ? this->exd->get_two_level_step() : this->exc->get_two_level_step();
    }

    // Peize Lin add 2016-12-03

    // Taoni Bao add 2026-05-15
    // In RT-TDDFT, contributeHk is used, but two_level_step is reset to 0 at each ionic step.
    // In order to add EXX correctly in for istep > 0, this->istep == 0 is needed to avoid skipping EXX calculation.
    // 1. For NSCF
    if (PARAM.inp.calculation == "nscf" || !has_workflow)
    {
        // Do nothing here, allow the code to proceed and calculate EXX.
    }
    // 2. For the first ionic step:
    else if (this->istep == 0)
    {
        // If EXX is once turned on (two_level_step > 0), let OperatorEXX remember this
        if (two_level_step > 0)
        {
            this->initial_gga_done = true;
        }

        // Check if we are in the pre-convergence stage of the two-level SCF (i.e., the pure GGA loop)
        bool in_gga_pre_loop = (two_level_step == 0);

        // Check if a high-quality initial guess is missing
        bool lacks_good_guess = (!this->restart);

        // If in the pre-convergence loop and lacking a good initial guess, skip adding the EXX contribution
        // Taoni Bao add 2026-05-18, only skip EXX if initial GGA loop is not done
        // Fix RT-TDDFT EXX missing problem in the evolution
        if (in_gga_pre_loop && lacks_good_guess && !this->initial_gga_done)
        {
            return; // In the non-EXX loop, skip adding EXX contribution
        }
    }
    // 3. For subsequent ionic steps (istep > 0), add EXX normally

    if (this->add_hexx_type == Add_Hexx_Type::R)
    {
        OperatorLCAO<TK, TR>::contributeHk(ik);
    }

    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        if (this->restart)
        {
            if (two_level_step == 0)
            {
                this->add_loaded_Hexx(ik);
                return;
            }
            else // clear loaded Hexx and release memory
            {
                if (this->Hexxd_k_load.size() > 0)
                {
                    this->Hexxd_k_load.clear();
                    this->Hexxd_k_load.shrink_to_fit();
                }
                else if (this->Hexxc_k_load.size() > 0)
                {
                    this->Hexxc_k_load.clear();
                    this->Hexxc_k_load.shrink_to_fit();
                }
            }
        }
        // cal H(k) from H(R) normally
        if (PARAM.inp.esolver_type == "tddft" && PARAM.inp.td_stype == 2)
        {
            RI_2D_Comm::add_Hexx_td(ucell,
                                    this->kv,
                                    ik,
                                    GlobalC::exx_info.info_global.hybrid_alpha,
                                    *this->Hexxc,
                                    *this->hR->get_paraV(),
                                    TD_info::td_vel_op->cart_At,
                                    TD_info::td_vel_op->get_phase_hybrid(),
                                    this->hsk->get_hk());
        }
        else
        {
            if (GlobalC::exx_info.info_ri.real_number)
            {
                RI_2D_Comm::add_Hexx(ucell,
                                     this->kv,
                                     ik,
                                     GlobalC::exx_info.info_global.hybrid_alpha,
                                     *this->Hexxd,
                                     *this->hR->get_paraV(),
                                     this->hsk->get_hk());
            }
            else
            {
                RI_2D_Comm::add_Hexx(ucell,
                                     this->kv,
                                     ik,
                                     GlobalC::exx_info.info_global.hybrid_alpha,
                                     *this->Hexxc,
                                     *this->hR->get_paraV(),
                                     this->hsk->get_hk());
            }
        }
    }
}

template <typename TK, typename TR>
template <typename Tdata>
void OperatorEXX<OperatorLCAO<TK, TR>>::cal_dH(
    const int ispin,
    std::array<std::vector<hamilt::HContainer<double>*>, 3>& dhR,
    const std::array<std::vector<std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>>, 3>& dHexxs)
{
    // dhR is the set of per-atom-I HContainers to fill (not this->hR, which may be a dummy here).
    const Parallel_Orbitals* const paraV = dhR[0][0]->get_paraV();
    const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest
        = this->use_cell_nearest ? &this->cell_nearest : nullptr;
    for (int idir = 0; idir < 3; ++idir)
    {
        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            // add_HexxR only fills existing matrices, so first allocate the atom-pair
            // structure of this per-I container from the exx-form data (same cell mapping).
            reallocate_hcontainer(dHexxs[idir][iat], dhR[idir][iat], cell_nearest);
            RI_2D_Comm::add_HexxR(ispin,
                GlobalC::exx_info.info_global.hybrid_alpha,
                dHexxs[idir][iat],
                *paraV,
                PARAM.globalv.npol,
                *dhR[idir][iat],
                cell_nearest);
        }
    }
}

} // namespace hamilt
#endif // __EXX
#endif // OPEXXLCAO_HPP
