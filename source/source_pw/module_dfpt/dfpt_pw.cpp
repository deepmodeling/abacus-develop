// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_pw.h"
#include "dfpt_pw_data.h"
#include "dfpt_irrep_data.h"
#include "dfpt_pert.h"
#include "dfpt_stern.h"
#include "dfpt_rho.h"
#include "dfpt_phon.h"
#include "dfpt_q0.h"
#include "dfpt_metal.h"
#include "source_cell/qlist.h"

namespace ModuleDFPT {

class DFPT_PW::Impl {
public:
    Impl() {}
    ~Impl() {}
    
    DFPT_PW_Data data_;
    DFPT_Pert pert_;
    DFPT_Stern stern_;
    DFPT_Rho rho_;
    DFPT_Phon phon_;
    DFPT_Q0 q0_;
    DFPT_Metal metal_;
    ModuleCell::QList qlist_;
    
    psi::Psi<std::complex<double>> gs_psi_;
    UnitCell* ucell_ = nullptr;
    double nelec_ = 0.0;
    double ecutwfc_ = 0.0;
    const Plus_U* dftu_ = nullptr;
    
    int nqx_ = 1, nqy_ = 1, nqz_ = 1;
    double conv_thr_ = 1e-8;
    int max_iter_ = 100;
};

DFPT_PW::DFPT_PW() : pimpl_(new Impl()) {}

DFPT_PW::~DFPT_PW() {
    delete pimpl_;
}

void DFPT_PW::init(UnitCell& ucell, const psi::Psi<std::complex<double>>& psi,
                   double nelec, double ecutwfc, const Plus_U* dftu) {
    pimpl_->ucell_ = &ucell;
    pimpl_->gs_psi_ = psi;
    pimpl_->nelec_ = nelec;
    pimpl_->ecutwfc_ = ecutwfc;
    pimpl_->dftu_ = dftu;
    
    std::vector<int> mp_grid = {pimpl_->nqx_, pimpl_->nqy_, pimpl_->nqz_};
    pimpl_->qlist_.generate_mesh(ucell, ucell.symm, mp_grid, true);
    
    int nq = pimpl_->qlist_.get_nq();
    int nk = psi.get_nk();
    int nbands = psi.get_nbands();
    int npw_max = psi.get_current_ngk();
    int nrxx = 0;
    int nspin = 1;
    int nat = ucell.nat;
    
    pimpl_->phon_.init(ucell);
    pimpl_->data_.init(&pimpl_->qlist_, nk, nbands, npw_max, nrxx, nspin, nat, dftu);
}

bool DFPT_PW::get_with_u() const {
    return pimpl_->data_.with_u();
}

bool DFPT_PW::get_u_active() const {
    return pimpl_->data_.u_active();
}

void DFPT_PW::run() {
    int nq = pimpl_->qlist_.get_nq();
    DFPT_IrrepData irrep_data(pimpl_->data_);
    for (int q_idx = 0; q_idx < nq; ++q_idx) {
        // Special handling for q=0 (uniform electric field responses):
        // The standard position operator r is ill-defined in periodic systems.
        // Developers should NOT pass a conventional position matrix. Instead,
        // matrix elements should be computed using the well-defined periodic
        // commutator [Ĥ_SCF, r̂]. This is implemented in DFPT_Q0 module.
        if (q_idx == 0) {
            pimpl_->q0_.compute_q0_response(pimpl_->data_);
        }

        // Per-irrep self-consistent loop: solve only the representative modes
        // of each little-group irrep. The irrep decomposition is exposed
        // through DFPT_IrrepData (option-2 signatures simulated by the interim
        // wrapper); the screening-potential perturbation, the Sternheimer
        // solver and the first-order density update are wired in subsequent
        // iterations.
        const int nirr = irrep_data.get_nirr(q_idx);
        for (int irrep = 0; irrep < nirr; ++irrep) {
            irrep_data.set_converged(q_idx, irrep, false);
            irrep_data.set_current_iter(q_idx, irrep, 0);
            while (!irrep_data.get_converged(q_idx, irrep)
                   && irrep_data.get_current_iter(q_idx, irrep) < pimpl_->max_iter_) {
                // 1. Compute the perturbation of the screening potential
                //    pimpl_->pert_.build_dv(q_idx, irrep, pimpl_->data_)
                // 2. Solve the Sternheimer equation for the representative
                //    modes of this irrep
                //    pimpl_->stern_.solve(q_idx, irrep, pimpl_->data_)
                // 3. Calculate the first-order density
                //    pimpl_->rho_.compute_drho(q_idx, irrep, pimpl_->data_)
                // 4. Check convergence and iterate until self-consistency
                irrep_data.add_residual(q_idx, irrep, 0.0);
                irrep_data.set_converged(q_idx, irrep, true);
            }
        }

        pimpl_->phon_.assemble(q_idx, pimpl_->data_);
        pimpl_->phon_.diagonalize(q_idx, pimpl_->data_);
    }
}

std::vector<double> DFPT_PW::get_phonon_freq(int q_idx) const {
    return pimpl_->data_.get_phon_freq(q_idx);
}

ModuleBase::matrix DFPT_PW::get_dielectric_tensor() const {
    return pimpl_->data_.get_dielectric();
}

ModuleBase::matrix DFPT_PW::get_born_charges(int atom_idx) const {
    return pimpl_->data_.get_born(atom_idx);
}

void DFPT_PW::set_parameters(const std::string& param_file) {
    (void)param_file;
}

void DFPT_PW::set_qmesh(int nqx, int nqy, int nqz) {
    pimpl_->nqx_ = nqx;
    pimpl_->nqy_ = nqy;
    pimpl_->nqz_ = nqz;
}

void DFPT_PW::set_conv_thr(double thr) {
    pimpl_->conv_thr_ = thr;
    pimpl_->data_.set_conv_thr(thr);
}

void DFPT_PW::set_max_iter(int max_iter) {
    pimpl_->max_iter_ = max_iter;
    pimpl_->data_.set_max_iter(max_iter);
}

} // namespace ModuleDFPT