#include "hamilt_bse.h"

#include "hamilt_bse_solver.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"

#include <cassert>

namespace BSE
{
template <typename T>
HamiltBSE<T>::HamiltBSE(const int& nspin,
                const int& naos,
                const std::vector<int>& nocc,
                const std::vector<int>& nvirt,
                const UnitCell& ucell_in,
                const std::vector<double>& orb_cutoff_in,
                const Grid_Driver& gd_in,
                const psi::Psi<T>& psi_in,
                const psi::Psi<T>& psi_glb_in,
                const ModuleBase::matrix& eig_gw_in,
                MolecularLRI<T>& mo_lri_in,
                std::weak_ptr<LR::PotHxcLR> pot_in,
                const K_Vectors& kv_in,
                const std::vector<Parallel_2D>& pX_in,// vector for spin, parallel as {nvirt, nocc}
                const Parallel_2D& pc_in, //parallel as {nbasis, nbands}
                const Parallel_Orbitals& pmat_in, //parallel as {nbasis, nbasis}
                const std::vector<std::string>& spin_types_in, //can be singlet and triplet
                const std::string& tda, // can be: "tda", "full", "both"
                const bool bse_ri_hartree_in,
                const bool bse_mem_save_in,
                const int bse_continue_in,
                const bool out_bse_ab_in,
                const std::string& out_dir_in,
                const std::string& readin_dir_in,
                const std::string& ri_hartree_benchmark_in)
    : nspin(nspin), naos(naos), nocc(nocc), nvirt(nvirt), ucell(ucell_in),
    orb_cutoff(orb_cutoff_in), gd(gd_in), psi_ks(psi_in), psi_ks_glb(psi_glb_in), eig_gw(eig_gw_in),
    mo_lri(mo_lri_in),
    pot(pot_in), kv(kv_in),
    pX(pX_in), pc(pc_in), pmat(pmat_in),
    spin_types(spin_types_in), bse_ri_hartree(bse_ri_hartree_in), bse_mem_save(bse_mem_save_in),
    bse_continue(bse_continue_in), out_bse_ab(out_bse_ab_in), out_dir(out_dir_in), readin_dir(readin_dir_in),
    ri_hartree_benchmark(ri_hartree_benchmark_in)
{
    ModuleBase::TITLE("BSE", "HamiltBSE");
    if (this->pX[0].get_local_size() == 0) {
        std::cerr<< "Warning: Parallel_2D in RANK "+std::to_string(GlobalV::MY_RANK) +" has no local size, please use less mpi." << std::endl;
        std::cerr<< " [File:"<<__FILE__<< ", Function: " << __FUNCTION__ << ", Line: " << __LINE__ << "]" << std::endl;
    }
    assert(naos == pmat.get_global_row_size() && naos == pmat.get_global_col_size());
    this->nk = this->nspin == 2 ? this->kv.get_nks() / 2 : this->kv.get_nks();
    this->ndim = nk * nocc[0] * nvirt[0];

    int nb2d;
    if (this->ndim > 1000)
    {
        nb2d = 64;
    }
    else if (this->ndim > 500)
    {
        nb2d = 32;
    }
    else if (this->ndim > 0)
    {
        nb2d = 1;
    }
    else throw std::runtime_error("ndim in HamiltBSE is zero or negative.");
    LR_Util::setup_2d_division(this->pA, nb2d, ndim, ndim
        #ifdef __MPI
            , this->pX[0].blacs_ctxt
        #endif
        );

    if (!this->bse_ri_hartree && this->ri_hartree_benchmark == "none")
    {
        this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1/*nspin*/, kv_in.kvec_d, nk);
        this->DM_trans->set_DMK_zero();
        LR_Util::initialize_DMR(*this->DM_trans, this->pmat, this->ucell, this->gd, this->orb_cutoff);
    }
    if (this->bse_mem_save) { assert(this->bse_continue == 0 && this->bse_ri_hartree); }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "HamiltBSE is ready to calculate");

    if (this->bse_continue >= 1) {
        BSE_Util::print_mem_estimate("V matrix of A", this->pA.get_local_size(), sizeof(T));
        this->VA_local.resize(this->pA.get_local_size(), 0.0);
        this->read_AB_matrix(this->readin_dir + "A_V_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", this->VA_local.data(), this->ndim, this->ndim);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read_V_for_A");
    }
    if (this->bse_continue >= 2) {
        BSE_Util::print_mem_estimate("W matrix of A", this->pA.get_local_size(), sizeof(T));
        this->WA_local.resize(this->pA.get_local_size(), 0.0);
        this->read_AB_matrix(this->readin_dir + "A_W_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", this->WA_local.data(), this->ndim, this->ndim);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read_W_for_A");
    }
    if (this->bse_continue >= 3) {
        BSE_Util::print_mem_estimate("V matrix of B", this->pA.get_local_size(), sizeof(T));
        this->VB_local.resize(this->pA.get_local_size(), 0.0);
        this->read_AB_matrix(this->readin_dir + "B_V_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", this->VB_local.data(), this->ndim, this->ndim);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read_V_for_B");
    }
    if (this->bse_continue >= 4) {
        BSE_Util::print_mem_estimate("W matrix of B", this->pA.get_local_size(), sizeof(T));
        this->WB_local.resize(this->pA.get_local_size(), 0.0);
        this->read_AB_matrix(this->readin_dir + "B_W_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", this->WB_local.data(), this->ndim, this->ndim);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read_W_for_B");
    }
    
    if (!this->bse_mem_save)
    {
        for (const auto& st : this->spin_types) {
            if (st == "singlet" || st == "triplet") {
                // Hartree term V (exchange electron and hole)
                if (st == "singlet"){
                    this->cal_V_for_A();
                    if (tda == "both" || tda == "full") { this->cal_V_for_B(); }
                }
                else if (st == "triplet") {
                    std::cout << "Hatree term is not needed for triplet." << std::endl;
                }

                // direct term W (electron-electron and hole-hole)
                this->cal_W_for_A();
                if (tda == "both" || tda == "full") { this->cal_W_for_B(); }
            }
            else if(st == "rpa") {
                this->cal_V_for_A();
                if (tda == "both" || tda == "full") { this->cal_V_for_B(); }
            }
            else if(st != "ipa") {
                throw std::runtime_error("Unsupported type in BSE: " + st);
            }
        }
        this->mo_lri.LR_lri.free_Vs();
        this->mo_lri.LR_lri.free_Ws();
        malloc_trim(0);
    }
}

template <typename T>
void HamiltBSE<T>::cal_V_for_A(){
    ModuleBase::TITLE("HamiltBSE", "cal_V_for_A");
    ModuleBase::timer::start("HamiltBSE", "cal_V_for_A");
    std::cout<<"in cal_V_for_A"<<std::endl;
    if (!this->VA_local.empty()) {
        std::cout<< "V for A has been calculated, skip." <<std::endl;
        return;
    }
    BSE_Util::print_mem_estimate("V matrix of A", this->pA.get_local_size(), sizeof(T));
    this->VA_local.resize(this->pA.get_local_size(), 0.0);
    if (this->ri_hartree_benchmark == "aims" || this->ri_hartree_benchmark == "abacus") {
        throw std::runtime_error("this BSE routine only supports aims-librpa/abacus-librpa benchmark");
    }
    else if (this->bse_ri_hartree || this->ri_hartree_benchmark =="aims-librpa" || this->ri_hartree_benchmark == "abacus-librpa") {
        std::cout << "Calculating Hartree term for A with RI approximation" << std::endl;
        this->mo_lri.cal_hartree_for_A(this->VA_local, this->pA);
    }
    else if (this->ri_hartree_benchmark == "none") { // do things like OperatorLRHxc
        std::cout << "Calculating Hartree term for A with grid integration" << std::endl;
        this->cal_V_by_grid(true);
    }
    if (this->out_bse_ab){
        this->write_AB_matrix(this->out_dir+"A_V_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", 6, this->VA_local.data(), this->ndim, this->ndim);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_V_for_A");
    ModuleBase::timer::end("HamiltBSE", "cal_V_for_A");
}

template <typename T>
void HamiltBSE<T>::cal_V_for_B(){
    ModuleBase::TITLE("HamiltBSE", "cal_V_for_B");
    ModuleBase::timer::start("HamiltBSE", "cal_V_for_B");
    std::cout<<"in cal_V_for_B"<<std::endl;
    if (!this->VB_local.empty()) {
        std::cout<< "V for B has been calculated, skip." <<std::endl;
        return;
    }
    BSE_Util::print_mem_estimate("V matrix of B", this->pA.get_local_size(), sizeof(T));
    this->VB_local.resize(this->pA.get_local_size(), 0.0);
    if (this->ri_hartree_benchmark == "aims" || this->ri_hartree_benchmark == "abacus") {
        throw std::runtime_error("this BSE routine only supports aims-librpa/abacus-librpa benchmark");
    }
    else if (this->bse_ri_hartree || this->ri_hartree_benchmark =="aims-librpa" || this->ri_hartree_benchmark == "abacus-librpa") {
        std::cout << "Calculating Hartree term for B with RI approximation" << std::endl;
        this->mo_lri.cal_hartree_for_B(this->VB_local, this->pA);
    }
    else if (this->ri_hartree_benchmark == "none") { // do things like OperatorLRHxc
        std::cout << "Calculating Hartree term for B with grid integration" << std::endl;
        this->cal_V_by_grid(false);
    }
    if (this->out_bse_ab){
        this->write_AB_matrix(this->out_dir+"B_V_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", 6, this->VB_local.data(), this->ndim, this->ndim);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_V_for_B");
    ModuleBase::timer::end("HamiltBSE", "cal_V_for_B");
}

template <typename T>
void HamiltBSE<T>::cal_W_for_A(){
    ModuleBase::TITLE("HamiltBSE", "cal_W_for_A");
    ModuleBase::timer::start("HamiltBSE", "cal_W_for_A");
    std::cout<<"in cal_W_for_A"<<std::endl;
    if (!this->WA_local.empty()) {
        std::cout<< "W for A has been calculated, skip." <<std::endl;
        return;
    }
    BSE_Util::print_mem_estimate("W matrix of A", this->pA.get_local_size(), sizeof(T));
    this->WA_local.resize(this->pA.get_local_size(), 0.0);
    this->mo_lri.cal_W_for_A(this->WA_local, this->pA);    
    if (this->out_bse_ab){
        this->write_AB_matrix(this->out_dir+"A_W_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", 6, this->WA_local.data(), this->ndim, this->ndim);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_W_for_A");
    ModuleBase::timer::end("HamiltBSE", "cal_W_for_A");
}

template <typename T>
void HamiltBSE<T>::cal_W_for_B(){
    ModuleBase::TITLE("HamiltBSE", "cal_W_for_B");
    ModuleBase::timer::start("HamiltBSE", "cal_W_for_B");
    std::cout<<"in cal_W_for_B"<<std::endl;
    if (!this->WB_local.empty()) {
        std::cout<< "W for B has been calculated, skip." <<std::endl;
        return;
    }
    BSE_Util::print_mem_estimate("W matrix of B", this->pA.get_local_size(), sizeof(T));
    this->WB_local.resize(this->pA.get_local_size(), 0.0);
    this->mo_lri.cal_W_for_B(this->WB_local, this->pA);
    
    if (this->out_bse_ab){
        this->write_AB_matrix(this->out_dir+"B_W_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", 6, this->WB_local.data(), this->ndim, this->ndim);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_W_for_B");
    ModuleBase::timer::end("HamiltBSE", "cal_W_for_B");
}


template <typename T>
void HamiltBSE<T>::init_bse_matrix(const bool is_full, const int & st_index){
    ModuleBase::TITLE("HamiltBSE", "init_bse_matrix");
    ModuleBase::timer::start("HamiltBSE", "init_bse_matrix");

    double alpha, beta;
    const std::string& st = this->spin_types[st_index];
    if (st == "singlet")      { alpha = 2.0; beta = -1.0; }
    else if (st == "triplet") { alpha = 0.0; beta = -1.0; }
    else if (st == "rpa")     { alpha = 2.0; beta =  0.0; }
    else if (st == "ipa")     { alpha = 0.0; beta =  0.0; }
    else { throw std::runtime_error("Unsupported type in BSE: " + st); }

    std::string tda_type = is_full ? "full" : "TDA";
    std::cout<<"| init "<< tda_type << " BSE for type: "<< st << std::endl;
    std::cout<<"| A(ai,bj) = (Ea-Ei) δ_ij δ_ab + α (ai|V|jb)  +  β (ji|W|ab)" << std::endl;
    std::cout<<"|   term coefficient: (Exchange) α: "<<std::setw(2)<<alpha<<", (Direct) β: "<<beta<<std::endl;
    if (is_full) {
        std::cout<<"| B(ai,jb) = α (ai|V|bj)  +  β (bi|W|aj)" << std::endl;
    }

    BSE_Util::print_mem_estimate("BSE A matrix", this->pA.get_local_size(), sizeof(T));
    this->BSE_A_local.assign(this->pA.get_local_size(), 0.0);
    if (is_full)
    {
        BSE_Util::print_mem_estimate("BSE B matrix", this->pA.get_local_size(), sizeof(T));
        this->BSE_B_local.assign(this->pA.get_local_size(), 0.0);
    }
    // 1) add diagonal GW energy term
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
    for(int ik = 0;ik < nk;++ik)
    {
        for (int i = 0;i < nocc[0];++i)
        {
            for(int a = 0;a < nvirt[0];++a)
            {
                int index = ik * nocc[0] * nvirt[0] + i * nvirt[0] + a;
                int col_loc = this->pA.global2local_col(index);
                int row_loc = this->pA.global2local_row(index);
                if (col_loc == -1 || row_loc == -1) continue;
                this->BSE_A_local[col_loc * this->pA.get_row_size() + row_loc]
                    = this->eig_gw(ik, nocc[0] + a) - this->eig_gw(ik, i);
            }
        }
    }
    // 2) add V/W contributions
    if (this->bse_mem_save)
    {
        std::cout << "| bse_mem_save is true, V and W matrix will be added to BSE matrix directly." << std::endl;
        if (alpha != 0.0) {
            this->mo_lri.cal_hartree_for_A(this->BSE_A_local, this->pA, alpha);
            if (is_full) { this->mo_lri.cal_hartree_for_B(this->BSE_B_local, this->pA, alpha); }
        }
        if (beta != 0.0) {
            this->mo_lri.cal_W_for_A(this->BSE_A_local, this->pA, beta);
            if (is_full) { this->mo_lri.cal_W_for_B(this->BSE_B_local, this->pA, beta); }
        }
        GlobalV::ofs_running << "| V and W matrix has been added." << std::endl;
        std::cout << "| V and W matrix has been added." << std::endl;
    }
    else
    {
        if (alpha != 0.0) {
            assert(this->VA_local.size() == this->BSE_A_local.size());
            if (is_full) assert(this->VB_local.size() == this->BSE_A_local.size());
        }
        if (beta != 0.0) {
            assert(this->WA_local.size() == this->BSE_A_local.size());
            if (is_full) assert(this->WB_local.size() == this->BSE_A_local.size());
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < this->BSE_A_local.size(); ++i) {
            if (alpha != 0.0) this->BSE_A_local[i] += alpha * this->VA_local[i];
            if (beta != 0.0)  this->BSE_A_local[i] += beta  * this->WA_local[i];
            if (is_full) {
                if (alpha != 0.0) this->BSE_B_local[i] +=(alpha * this->VB_local[i]);
                if (beta != 0.0)  this->BSE_B_local[i] +=(beta * this->WB_local[i]);
            }
        }
    }    
    // 3) check hermiticity/symmetry and (optionally) write to file
    GlobalV::ofs_running << "CHECK hermiticity/symmetry" << std::endl;
    std::cout << "| CHECK hermiticity/symmetry" << std::endl;
    constexpr double threshold = 1.0e-6;
    if (LR_Util::is_hermitian(this->BSE_A_local.data(), this->pA, threshold))
    {
        std::cout << "|  CHECK PASS: Matrix A is hermitian under threshold " << threshold << std::endl;
    }
    else
    {
        std::cout << "|  CHECK WARNING: Matrix A is not hermitian under threshold " << threshold << std::endl;
    }
    if (this->out_bse_ab)
    {
        this->write_AB_matrix("A_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", 6, this->BSE_A_local.data(), this->ndim, this->ndim);
    }

    if (is_full)
    {
        if (LR_Util::is_symmetric(this->BSE_B_local.data(), this->pA, threshold))
        {
            std::cout << "|  CHECK PASS: Matrix B is symmetric under threshold " << threshold << std::endl;
        }
        else
        {
            std::cout << "|  CHECK WARNING: Matrix B is not symmetric under threshold " << threshold << std::endl;
        }
        if (this->out_bse_ab)
        {
            this->write_AB_matrix("B_matrix_"+std::to_string(GlobalV::MY_RANK)+".dat", 6, this->BSE_B_local.data(), this->ndim, this->ndim);
        }
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "init_bse_matrix");
    ModuleBase::timer::end("HamiltBSE", "init_bse_matrix");
}

template <typename T>
void HamiltBSE<T>::tda_solver(const int & st_index, const int& nstates, double* ene_out, T* X_out){
    ModuleBase::TITLE("HamiltBSE", "tda_solver");
    ModuleBase::timer::start("HamiltBSE", "elpa_tda_solver");

    this->init_bse_matrix(false, st_index);

    std::vector<T> X_tda(this->pA.get_local_size(), 0.0);
    std::vector<double> ev(nstates, 0.0);

    BSE::solve_tda(GlobalV::MY_RANK,
                    this->BSE_A_local,
                    this->pA,
                    ev,
                    X_tda);
    // copy to output
    std::copy_n(ev.data(), nstates, ene_out);

    LR_Util::pA2pX(X_out, X_tda.data(), nstates, this->nk,
                    this->nocc, this->nvirt, this->pX, this->pA, 0, 0, false/*openshell*/);
    
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "BSE TDA solver");
    ModuleBase::timer::end("HamiltBSE", "elpa_tda_solver");
}

template <>
void HamiltBSE<double>::full_solver(const int& st_index, const int& nstates,
                                    double* ene_out,
                                    double* X_out,
                                    double* Y_out){
    ModuleBase::TITLE("HamiltBSE", "full_solver(double)");
    ModuleBase::timer::start("HamiltBSE", "elpa_full_solver(double)");

    this->init_bse_matrix(true, st_index);
    Parallel_2D pM;
    LR_Util::setup_2d_division(pM, this->pA.get_block_size(), 2*this->ndim, 2*this->ndim
        #ifdef __MPI
                , this->pA.blacs_ctxt
        #endif
            );
    const auto to_complex = [](const std::vector<double>& vin) {
        return std::vector<std::complex<double>>(vin.begin(), vin.end());
    };
    std::vector<std::complex<double>> BSE_A_local_complex = to_complex(this->BSE_A_local);
    std::vector<std::complex<double>> BSE_B_local_complex = to_complex(this->BSE_B_local);
    std::vector<std::complex<double>> local_v_full(pM.get_local_size(), 0.0);
    std::vector<double> ev(2 * this->ndim, 0.0);
    BSE::solve_full(GlobalV::MY_RANK,
                    BSE_A_local_complex,
                    BSE_B_local_complex,
                    this->pA,
                    pM,
                    ev,
                    local_v_full);

    std::vector<double> local_v_full_real(pM.get_local_size(), 0.0);
    for (size_t i = 0; i < local_v_full.size(); ++i) {
        assert(std::abs(local_v_full[i].imag()) < 1e-10);
        local_v_full_real[i] = local_v_full[i].real();
    }

    // copy positive eigenvalues to output
    std::copy_n(&ev[this->ndim], nstates, ene_out);
    LR_Util::pA2pX(X_out, local_v_full_real.data(), nstates, this->nk,
                    this->nocc, this->nvirt, this->pX, pM, 0, this->ndim, false/*openshell*/);
    LR_Util::pA2pX(Y_out, local_v_full_real.data(), nstates, this->nk,
                    this->nocc, this->nvirt, this->pX, pM, this->ndim, this->ndim, false/*openshell*/);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "BSE Full solver");
    ModuleBase::timer::end("HamiltBSE", "elpa_full_solver(double)");
}

template <>
void HamiltBSE<std::complex<double>>::full_solver(const int& st_index, const int& nstates,
                                                    double* ene_out,
                                                    std::complex<double>* X_out,
                                                    std::complex<double>* Y_out){
    ModuleBase::TITLE("HamiltBSE", "full_solver(complex)");
    ModuleBase::timer::start("HamiltBSE", "elpa_full_solver(complex)");

    this->init_bse_matrix(true, st_index);
    Parallel_2D pM;
    LR_Util::setup_2d_division(pM, this->pA.get_block_size(), 2*this->ndim, 2*this->ndim
        #ifdef __MPI
                , this->pA.blacs_ctxt
        #endif
            );
    
    std::vector<std::complex<double>> local_v_full(pM.get_local_size(), 0.0);
    std::vector<double> ev(2 * this->ndim, 0.0);
    ModuleBase::TITLE("HamiltBSE", "full_solver(complex)2");
    BSE::solve_full(GlobalV::MY_RANK,
                    this->BSE_A_local,
                    this->BSE_B_local,
                    this->pA,
                    pM,
                    ev,
                    local_v_full);
    ModuleBase::TITLE("HamiltBSE", "full_solver(complex)3");
    // copy positive eigenvalues to output
    std::copy_n(&ev[this->ndim], nstates, ene_out);
    LR_Util::pA2pX(X_out, local_v_full.data(), nstates, this->nk,
                    this->nocc, this->nvirt, this->pX, pM, 0, this->ndim, false/*openshell*/);
    LR_Util::pA2pX(Y_out, local_v_full.data(), nstates, this->nk,
                    this->nocc, this->nvirt, this->pX, pM, this->ndim, this->ndim, false/*openshell*/);
    
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "BSE FULL solver");
    ModuleBase::timer::end("HamiltBSE", "elpa_full_solver(complex)");
}

template <typename T>
void HamiltBSE<T>::cal_V_by_grid(bool is_A)
{
    // 1. initialize HContainer VR
    const int is = 0; //spin index, NOTE: only support 1 spin now
    const auto psi_is = LR_Util::get_psi_spin(psi_ks, is, nk);
    std::unique_ptr<hamilt::HContainer<T>> VR = std::unique_ptr<hamilt::HContainer<T>>(new hamilt::HContainer<T>(&this->pmat));
    LR_Util::initialize_HR<T, T>(*VR, this->ucell, this->gd, this->orb_cutoff);
#ifdef __MPI
    Parallel_2D pV_col; //{nvirt*nocc, 1}
    LR_Util::setup_2d_division(pV_col, this->pX[is].get_block_size(), nvirt[is]*nocc[is], 1, this->pX[is].blacs_ctxt);
    std::vector<T> V_col_local(this->nk * pV_col.get_local_size(), 0.0); // V_col(bjk2)
#endif

    int imo1, imo2;
    for (int ik2 = 0; ik2 < nk; ++ik2) {                    
        for (int j = 0; j < nocc[is]; ++j) {
            for (int b = 0; b < nvirt[is]; ++b) {//calculate row {aik1} for each column {bjk2}
                int bjk = ik2 * nocc[is] * nvirt[is] + j * nvirt[is] + b; // column index in BSE matrix
                // 2. calculate transition matrix 
                if (is_A) //jk2→bk2, D(k)=c_b(k)c^†_j(k)
                    { imo1 = j; imo2 = b + nocc[is]; }
                else //bjk2←kb2, D(k)=c_j(k)c^†_b(k)
                    { imo1 = b + nocc[is]; imo2 = j; }
    #ifdef __MPI
                ct::Tensor dm_trans_2d = 
                    BSE_Util::cal_dm_trans_onebase_pblas(psi_is, pc, ik2, naos, imo1, imo2, pmat, (T)1.0 / (T)nk);
    #else
                ct::Tensor dm_trans_2d = 
                    BSE_Util::cal_dm_trans_onebase_blas(psi_is, ik2, naos, imo1, imo2, (T)1.0 / (T)nk);
    #endif
                // LR_Util::print_tensor<T>(dm_trans_2d, "dm_trans_2d", &pmat);
                this->DM_trans->set_DMK_pointer(ik2, dm_trans_2d.data<T>());
                // 3. D(k)→D(R)
                this->DM_trans->cal_DMR(ik2);
                // LR_Util::print_DMR(*DM_trans, ucell.nat, "DMR");

                // 4. D(R)→V(R)
                this->grid_calculation(*VR);

                // 5. V(R)→V(k) 
                std::vector<ct::Tensor> v_k_2d(nk, LR_Util::newTensor<T>({ pmat.get_col_size(), pmat.get_row_size() }));
                for (auto& v : v_k_2d) v.zero();
                int nrow = this->pmat.get_row_size();
                for (int ik1 = 0;ik1 < nk;++ik1) {
                    folding_HR(*VR, v_k_2d[ik1].data<T>(), this->kv.kvec_d[ik1], nrow, 1);
                }
                // for (int ik1 = 0;ik1 < nk;++ik1)
                //     LR_Util::print_tensor<T>(v_k_2d[ik1], "V(k)[ik=" + std::to_string(ik1) + "]", &this->pmat);
    #ifdef __MPI
                LR::ao_to_mo_pblas(v_k_2d, this->pmat, psi_is, this->pc, this->naos,
                                nocc[is], nvirt[is], pV_col, V_col_local.data(), false, LR_Util::MO_TYPE::VO);

                for (int ik1 = 0; ik1 < this->nk; ++ik1) {
                    Cpxgemr2d(nvirt[is]*nocc[is], 1,
                            V_col_local.data() + ik1 * pV_col.get_local_size(), 1, 1, pV_col.desc,
                            this->VA_local.data(), 
                            ik1 * nocc[is] * nvirt[is] + 1 , bjk + 1, this->pA.desc,
                            this->pA.blacs_ctxt);
                }
    #else
                LR::ao_to_mo_blas(v_k_2d, psi_is, nocc[is], nvirt[is], this->VA_local.data()+bjk * this->ndim, false, LR_Util::MO_TYPE::VO);
    #endif
            }
        }
    }
}

template<>
void HamiltBSE<double>::grid_calculation(hamilt::HContainer<double>& VR) const
{
    ModuleBase::TITLE("HamiltBSE", "grid_calculation(double)");
    ModuleBase::timer::start("HamiltBSE", "grid_calculation(double)");

    // 4.1. transition density rho on grid
    double** rho_trans;
    const int& nrxx = this->pot.lock()->nrxx;

    LR_Util::_allocate_2order_nested_ptr(rho_trans, 1, nrxx); // nspin=1 for transition density
    ModuleBase::GlobalFunc::ZEROS(rho_trans[0], nrxx);
    ModuleGint::cal_gint_rho(this->DM_trans->get_DMR_vector(), 1, rho_trans, false);

    // 4.2. v_hxc = f_hxc * rho_trans
    ModuleBase::matrix vr_hxc(1, nrxx);   //grid
    std::vector<int> ispin_ks = { 0 }; //for close-shell dft-xc kerenl, actually placeholder for bse 
    this->pot.lock()->cal_v_eff(rho_trans, ucell, vr_hxc, ispin_ks);// in this function, unit changes from Ha to Ry
    LR_Util::_deallocate_2order_nested_ptr(rho_trans, 1);

    // 4.3 V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\nu(r)
    VR.set_zero();
    ModuleGint::cal_gint_vl(vr_hxc.c, &VR);
    // LR_Util::print_HR(VR, this->ucell.nat, "VR(real, 2d)");

    ModuleBase::timer::end("HamiltBSE", "grid_calculation(double)");
}

template<>
void HamiltBSE<std::complex<double>>::grid_calculation(hamilt::HContainer<std::complex<double>>& VR) const
{
    ModuleBase::TITLE("HamiltBSE", "grid_calculation(complex)");
    ModuleBase::timer::start("HamiltBSE", "grid_calculation(complex)");

    elecstate::DensityMatrix<std::complex<double>, double> DM_trans_real_imag(&this->pmat, 1, this->kv.kvec_d, this->nk);
    DM_trans_real_imag.init_DMR(VR);
    hamilt::HContainer<double> HR_real_imag(ucell, &this->pmat);
    LR_Util::initialize_HR<std::complex<double>, double>(HR_real_imag, ucell, gd, orb_cutoff);

    auto dmR_to_hR = [&, this](const char& type) -> void
        {
            LR_Util::get_DMR_real_imag_part(*this->DM_trans, DM_trans_real_imag, ucell.nat, type);
            // if (this->first_print)LR_Util::print_DMR(DM_trans_real_imag, ucell.nat, "DMR(2d, real)");

            // 4.1. transition density rho on grid
            double** rho_trans;
            const int& nrxx = this->pot.lock()->nrxx;

            LR_Util::_allocate_2order_nested_ptr(rho_trans, 1, nrxx); // nspin=1 for transition density
            ModuleBase::GlobalFunc::ZEROS(rho_trans[0], nrxx);
            ModuleGint::cal_gint_rho(DM_trans_real_imag.get_DMR_vector(), 1, rho_trans, false);

            // 4.2. v_hxc = f_hxc * rho_trans
            ModuleBase::matrix vr_hxc(1, nrxx);   //grid
            std::vector<int> ispin_ks = { 0 }; //for close-shell dft-xc kerenl, actually placeholder for bse 
            this->pot.lock()->cal_v_eff(rho_trans, ucell, vr_hxc, ispin_ks);// in this function, unit changes from Ha to Ry
            LR_Util::_deallocate_2order_nested_ptr(rho_trans, 1);

            // 4.3 V^{Hxc}_{\mu,\nu}=\int{dr} \phi_\mu(r) v_{Hxc}(r) \phi_\nu(r)
            HR_real_imag.set_zero();
            ModuleGint::cal_gint_vl(vr_hxc.c, &HR_real_imag);
            // LR_Util::print_HR(HR_real_imag, this->ucell.nat, "VR(real, 2d)");
            LR_Util::set_HR_real_imag_part(HR_real_imag, VR, ucell.nat, type);
        };
    VR.set_zero();
    dmR_to_hR('R');   //real
    if (this->nk > 1) { dmR_to_hR('I'); }   //imag for multi-k
    ModuleBase::timer::end("HamiltBSE", "grid_calculation(complex)");
}

template class HamiltBSE<double>;
template class HamiltBSE<std::complex<double>>;
}//namespace BSE