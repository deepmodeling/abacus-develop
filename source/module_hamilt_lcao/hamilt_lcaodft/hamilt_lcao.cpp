#include "hamilt_lcao.h"

#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include <vector>

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "operator_lcao/deepks_lcao.h"
#endif

#ifdef __EXX
#include "module_ri/Exx_LRI_interface.h"
#include "operator_lcao/op_exx_lcao.h"
#endif

#ifdef __ELPA
#include "module_hsolver/diago_elpa.h"
#endif

#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hsolver/hsolver_lcao.h"
#include "operator_lcao/dftu_lcao.h"
#include "operator_lcao/ekinetic_new.h"
#include "operator_lcao/meta_lcao.h"
#include "operator_lcao/nonlocal_new.h"
#include "operator_lcao/op_dftu_lcao.h"
#include "operator_lcao/op_exx_lcao.h"
#include "operator_lcao/overlap_new.h"
#include "operator_lcao/dspin_lcao.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "operator_lcao/td_ekinetic_lcao.h"
#include "operator_lcao/td_nonlocal_lcao.h"
#include "operator_lcao/veff_lcao.h"

namespace hamilt
{

template <typename TK, typename TR>
HamiltLCAO<TK, TR>::HamiltLCAO(const Parallel_Orbitals* paraV, const K_Vectors& kv_in, const TwoCenterIntegrator& intor_overlap_orb, const std::vector<double>& orb_cutoff)
{
    this->classname = "HamiltLCAO";

    this->kv = &kv_in;

    // Real space Hamiltonian is inited with template TR
    this->hR = new HContainer<TR>(paraV);
    this->sR = new HContainer<TR>(paraV);
    this->hsk = new HS_Matrix_K<TK>(paraV);

    this->getOperator() = new OverlapNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                               this->kv->kvec_d,
                                                               this->hR,
                                                               this->sR,
                                                               &GlobalC::ucell,
                                                               orb_cutoff,
                                                               &GlobalC::GridD,
                                                               &intor_overlap_orb);
}

template <typename TK, typename TR>
HamiltLCAO<TK, TR>::HamiltLCAO(Gint_Gamma* GG_in,
    Gint_k* GK_in,
    const Parallel_Orbitals* paraV,
    elecstate::Potential* pot_in,
    const K_Vectors& kv_in,
    const TwoCenterBundle& two_center_bundle,
    const LCAO_Orbitals& orb,
    elecstate::DensityMatrix<TK, double>* DM_in
#ifdef __EXX
    , const int istep
    , int* exx_two_level_step
    , std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd
    , std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
) {
    this->kv = &kv_in;
    this->classname = "HamiltLCAO";

    // Real space Hamiltonian is inited with template TR
    this->hR = new HContainer<TR>(paraV);
    this->sR = new HContainer<TR>(paraV);
    this->hsk = new HS_Matrix_K<TK>(paraV);

    // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>) is registered without template
    std::vector<std::string> pot_register_in;
    if (PARAM.inp.vl_in_h)
    {
        if (PARAM.inp.vion_in_h)
        {
            pot_register_in.push_back("local");
        }
        if (PARAM.inp.vh_in_h)
        {
            pot_register_in.push_back("hartree");
        }
        pot_register_in.push_back("xc");
        if (PARAM.inp.imp_sol)
        {
            pot_register_in.push_back("surchem");
        }
        if (PARAM.inp.efield_flag)
        {
            pot_register_in.push_back("efield");
        }
        if (PARAM.inp.gate_flag)
        {
            pot_register_in.push_back("gatefield");
        }
        if (PARAM.inp.esolver_type == "tddft")
        {
            pot_register_in.push_back("tddft");
        }
    }

    // Gamma_only case to initialize HamiltLCAO
    //
    // code block to construct Operator Chains
    if (std::is_same<TK, double>::value)
    {
        // fix HR to gamma case, where SR will be fixed in Overlap Operator
        this->hR->fix_gamma();
        // initial operator for Gamma_only case
        // overlap term (<psi|psi>) is indispensable
        // in Gamma_only case, target SK is this->hsk->get_sk(), the target SR is this->sR
        this->getOperator() = new OverlapNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                                   this->kv->kvec_d,
                                                                   this->hR,
                                                                   this->sR,
                                                                   &GlobalC::ucell,
                                                                   orb.cutoffs(),
                                                                   &GlobalC::GridD,
                                                                   two_center_bundle.overlap_orb.get());

        // kinetic term (<psi|T|psi>)
        if (PARAM.inp.t_in_h)
        {
            Operator<TK>* ekinetic = new EkineticNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                                           this->kv->kvec_d,
                                                                           this->hR,
                                                                           &GlobalC::ucell,
                                                                           orb.cutoffs(),
                                                                           &GlobalC::GridD,
                                                                           two_center_bundle.kinetic_orb.get());
            this->getOperator()->add(ekinetic);
        }

        // nonlocal term (<psi|beta>D<beta|psi>)
        // in general case, target HR is this->hR, while target HK is this->hsk->get_hk()
        if (PARAM.inp.vnl_in_h)
        {
            Operator<TK>* nonlocal = new NonlocalNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                                           this->kv->kvec_d,
                                                                           this->hR,
                                                                           &GlobalC::ucell,
                                                                           orb.cutoffs(),
                                                                           &GlobalC::GridD,
                                                                           two_center_bundle.overlap_orb_beta.get());
            this->getOperator()->add(nonlocal);
        }

        // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
        // in general case, target HR is Gint::hRGint, while target HK is this->hsk->get_hk()
        if (PARAM.inp.vl_in_h)
        {
            // only Potential is not empty, Veff and Meta are available
            if (pot_register_in.size() > 0)
            {
                // register Potential by gathered operator
                pot_in->pot_register(pot_register_in);
                // effective potential term
                Operator<TK>* veff = new Veff<OperatorLCAO<TK, TR>>(GG_in,
                    this->hsk,
                    this->kv->kvec_d,
                    pot_in,
                    this->hR, // no explicit call yet
                    &GlobalC::ucell,
                    orb.cutoffs(),
                    &GlobalC::GridD,
                    PARAM.inp.nspin
                );
                this->getOperator()->add(veff);
            }
        }

#ifdef __DEEPKS
        if (PARAM.inp.deepks_scf)
        {
            Operator<TK>* deepks = new DeePKS<OperatorLCAO<TK, TR>>(this->hsk,
                                                                    this->kv->kvec_d,
                                                                    this->hR, // no explicit call yet
                                                                    &GlobalC::ucell,
                                                                    &GlobalC::GridD,
                                                                    two_center_bundle.overlap_orb_alpha.get(),
                                                                    &orb,
                                                                    this->kv->get_nks(),
                                                                    DM_in);
            this->getOperator()->add(deepks);
        }
#endif

        // end node should be OperatorDFTU
        if (PARAM.inp.dft_plus_u)
        {
            Operator<TK>* dftu = nullptr;
            if (PARAM.inp.dft_plus_u == 2)
            {
                dftu = new OperatorDFTU<OperatorLCAO<TK, TR>>(this->hsk,
                                                              kv->kvec_d,
                                                              this->hR, // no explicit call yet
                                                              this->kv->isk);
            }
            else
            {
                dftu = new DFTU<OperatorLCAO<TK, TR>>(this->hsk,
                                                      this->kv->kvec_d,
                                                      this->hR,
                                                      GlobalC::ucell,
                                                      &GlobalC::GridD,
                                                      two_center_bundle.overlap_orb_onsite.get(),
                                                      orb.cutoffs(),
                                                      &GlobalC::dftu);
            }
            this->getOperator()->add(dftu);
        }
    }
    // multi-k-points case to initialize HamiltLCAO, ops will be used
    else if (std::is_same<TK, std::complex<double>>::value)
    {
        // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
        // Meta potential term (\sum_r <psi(r)|tau(r)|psi(r)>)
        // in general case, target HR is Gint::pvpR_reduced, while target HK is this->hsk->get_hk()
        if (PARAM.inp.vl_in_h)
        {
            // only Potential is not empty, Veff and Meta are available
            if (pot_register_in.size() > 0)
            {
                // register Potential by gathered operator
                pot_in->pot_register(pot_register_in);
                // Veff term
                this->getOperator() = new Veff<OperatorLCAO<TK, TR>>(GK_in,
                    this->hsk,
                    kv->kvec_d,
                    pot_in,
                    this->hR,
                    &GlobalC::ucell,
                    orb.cutoffs(),
                    &GlobalC::GridD,
                    PARAM.inp.nspin);
                // reset spin index and real space Hamiltonian matrix
                int start_spin = -1;
                GK_in->reset_spin(start_spin);
            }
        }

        // initial operator for multi-k case
        // overlap term is indispensable
        Operator<TK>* overlap = new OverlapNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                                     this->kv->kvec_d,
                                                                     this->hR,
                                                                     this->sR,
                                                                     &GlobalC::ucell,
                                                                     orb.cutoffs(),
                                                                     &GlobalC::GridD,
                                                                     two_center_bundle.overlap_orb.get());
        if (this->getOperator() == nullptr)
        {
            this->getOperator() = overlap;
        }
        else
        {
            this->getOperator()->add(overlap);
        }

        // kinetic term (<psi|T|psi>),
        // in general case, target HR is this->hR, while target HK is this->hsk->get_hk()
        if (PARAM.inp.t_in_h)
        {
            Operator<TK>* ekinetic = new EkineticNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                                           this->kv->kvec_d,
                                                                           this->hR,
                                                                           &GlobalC::ucell,
                                                                           orb.cutoffs(),
                                                                           &GlobalC::GridD,
                                                                           two_center_bundle.kinetic_orb.get());
            this->getOperator()->add(ekinetic);
        }

        // nonlocal term (<psi|beta>D<beta|psi>)
        // in general case, target HR is this->hR, while target HK is this->hsk->get_hk()
        if (PARAM.inp.vnl_in_h)
        {
            Operator<TK>* nonlocal = new NonlocalNew<OperatorLCAO<TK, TR>>(this->hsk,
                                                                           this->kv->kvec_d,
                                                                           this->hR,
                                                                           &GlobalC::ucell,
                                                                           orb.cutoffs(),
                                                                           &GlobalC::GridD,
                                                                           two_center_bundle.overlap_orb_beta.get());
            // TDDFT velocity gague will calculate full non-local potential including the original one and the
            // correction on its own. So the original non-local potential term should be skipped
            if (PARAM.inp.esolver_type != "tddft" || elecstate::H_TDDFT_pw::stype != 1)
            {
                this->getOperator()->add(nonlocal);
            }
            else
            {
                delete nonlocal;
            }
        }

#ifdef __DEEPKS
        if (PARAM.inp.deepks_scf)
        {
            Operator<TK>* deepks = new DeePKS<OperatorLCAO<TK, TR>>(this->hsk,
                                                                    this->kv->kvec_d,
                                                                    hR,
                                                                    &GlobalC::ucell,
                                                                    &GlobalC::GridD,
                                                                    two_center_bundle.overlap_orb_alpha.get(),
                                                                    &orb,
                                                                    this->kv->get_nks(),
                                                                    DM_in);
            this->getOperator()->add(deepks);
        }
#endif
        // TDDFT_velocity_gague
        if (TD_Velocity::tddft_velocity)
        {
            if(!TD_Velocity::init_vecpot_file) { elecstate::H_TDDFT_pw::update_At();
}
            Operator<TK>* td_ekinetic = new TDEkinetic<OperatorLCAO<TK, TR>>(this->hsk,
                                                                             this->hR,
                                                                             kv,
                                                                             &GlobalC::ucell,
                                                                             orb.cutoffs(),
                                                                             &GlobalC::GridD,
                                                                             two_center_bundle.overlap_orb.get());
            this->getOperator()->add(td_ekinetic);

            Operator<TK>* td_nonlocal = new TDNonlocal<OperatorLCAO<TK, TR>>(this->hsk,
                                                                             this->kv->kvec_d,
                                                                             this->hR,
                                                                             &GlobalC::ucell,
                                                                             orb,
                                                                             &GlobalC::GridD);
            this->getOperator()->add(td_nonlocal);
        }
        if (PARAM.inp.dft_plus_u)
        {
            Operator<TK>* dftu = nullptr;
            if (PARAM.inp.dft_plus_u == 2)
            {
                dftu = new OperatorDFTU<OperatorLCAO<TK, TR>>(this->hsk,
                                                              kv->kvec_d,
                                                              this->hR, // no explicit call yet
                                                              this->kv->isk);
            }
            else
            {
                dftu = new DFTU<OperatorLCAO<TK, TR>>(this->hsk,
                                                      this->kv->kvec_d,
                                                      this->hR,
                                                      GlobalC::ucell,
                                                      &GlobalC::GridD,
                                                      two_center_bundle.overlap_orb_onsite.get(),
                                                      orb.cutoffs(),
                                                      &GlobalC::dftu);
            }
            this->getOperator()->add(dftu);
        }
        if (PARAM.inp.sc_mag_switch)
        {
            Operator<TK>* sc_lambda = 
                new DeltaSpin<OperatorLCAO<TK, TR>>(
                    this->hsk,
                    kv->kvec_d,
                    this->hR,
                    GlobalC::ucell,
                    &GlobalC::GridD,
                    two_center_bundle.overlap_orb_onsite.get(),
                    orb.cutoffs()
            );
            this->getOperator()->add(sc_lambda);
            spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
            sc.set_operator(sc_lambda);
        }
    }

#ifdef __EXX
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        // Peize Lin add 2016-12-03
        // set xc type before the first cal of xc in pelec->init_scf
        // and calculate Cs, Vs
        Operator<TK>* exx = new OperatorEXX<OperatorLCAO<TK, TR>>(this->hsk,
            this->hR,
            *this->kv,
            Hexxd,
            Hexxc,
            Add_Hexx_Type::R,
            istep,
            exx_two_level_step,
            !GlobalC::restart.info_load.restart_exx
            && GlobalC::restart.info_load.load_H);
        this->getOperator()->add(exx);
    }
#endif
    // if NSPIN==2, HR should be separated into two parts, save HR into this->hRS2
    int memory_fold = 1;
    if (PARAM.inp.nspin == 2)
    {
        this->hRS2.resize(this->hR->get_nnr() * 2);
        this->hR->allocate(this->hRS2.data(), 0);
        memory_fold = 2;
    }

    ModuleBase::Memory::record("HamiltLCAO::hR", this->hR->get_memory_size() * memory_fold);
    ModuleBase::Memory::record("HamiltLCAO::sR", this->sR->get_memory_size());

    return;
}

// case for multi-k-points
template <typename TK, typename TR>
void HamiltLCAO<TK, TR>::matrix(MatrixBlock<TK>& hk_in, MatrixBlock<TK>& sk_in)
{
    auto op = dynamic_cast<OperatorLCAO<TK, TR>*>(this->getOperator());
    assert(op != nullptr);
    op->matrixHk(hk_in, sk_in);
}

template <typename TK, typename TR>
void HamiltLCAO<TK, TR>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltLCAO", "updateHk");
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
    // update global spin index
    if (PARAM.inp.nspin == 2)
    {
        // if Veff is added and current_spin is changed, refresh HR
        if (PARAM.inp.vl_in_h && this->kv->isk[ik] != this->current_spin)
        {
            // change data pointer of HR
            this->hR->allocate(this->hRS2.data() + this->hRS2.size() / 2 * this->kv->isk[ik], 0);
            if (this->refresh_times > 0)
            {
                this->refresh_times--;
                dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(this->ops)->set_hr_done(false);
            }
        }
        this->current_spin = this->kv->isk[ik];
    }
    this->getOperator()->init(ik);
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
}

template <typename TK, typename TR>
void HamiltLCAO<TK, TR>::refresh()
{
    ModuleBase::TITLE("HamiltLCAO", "refresh");
    dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(this->ops)->set_hr_done(false);
    if (PARAM.inp.nspin == 2)
    {
        this->refresh_times = 1;
        this->current_spin = 0;
        if (this->hR->get_nnr() != this->hRS2.size() / 2)
        {
            // operator has changed, resize hRS2
            this->hRS2.resize(this->hR->get_nnr() * 2);
        }
        this->hR->allocate(this->hRS2.data(), 0);
    }
}

// get Operator base class pointer
template <typename TK, typename TR>
Operator<TK>*& HamiltLCAO<TK, TR>::getOperator()
{
    return this->ops;
}

template <typename TK, typename TR>
void HamiltLCAO<TK, TR>::updateSk(const int ik, const int hk_type)
{
    ModuleBase::TITLE("HamiltLCAO", "updateSk");
    ModuleBase::timer::tick("HamiltLCAO", "updateSk");
    ModuleBase::GlobalFunc::ZEROS(this->getSk(), this->get_size_hsk());
    if (hk_type == 1) // collumn-major matrix for SK
    {
        const int nrow = this->hsk->get_pv()->get_row_size();
        hamilt::folding_HR(*this->sR, this->getSk(), this->kv->kvec_d[ik], nrow, 1);
    }
    else if (hk_type == 0) // row-major matrix for SK
    {
        const int ncol = this->hsk->get_pv()->get_col_size();
        hamilt::folding_HR(*this->sR, this->getSk(), this->kv->kvec_d[ik], ncol, 0);
    }
    ModuleBase::timer::tick("HamiltLCAO", "updateSk");
}

// case for nspin<4, gamma-only k-point
template class HamiltLCAO<double, double>;
// case for nspin<4, multi-k-points
template class HamiltLCAO<std::complex<double>, double>;
// case for nspin == 4, non-collinear spin case
template class HamiltLCAO<std::complex<double>, std::complex<double>>;
} // namespace hamilt
