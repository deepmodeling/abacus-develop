#include "hamilt_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "operator_lcao/deepks_lcao.h"
#endif
#ifdef __ELPA
#include "module_hsolver/diago_elpa.h"
#endif
#include "operator_lcao/op_dftu_lcao.h"
#include "operator_lcao/meta_lcao.h"
#include "operator_lcao/op_exx_lcao.h"
#include "operator_lcao/overlap_new.h"
#include "operator_lcao/ekinetic_new.h"
#include "operator_lcao/nonlocal_new.h"
#include "operator_lcao/veff_lcao.h"
#include "operator_lcao/sc_lambda_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

namespace hamilt
{

template<typename TK, typename TR>
HamiltLCAO<TK, TR>::HamiltLCAO(LCAO_Matrix* LM_in, const K_Vectors& kv_in)
{
    this->classname = "HamiltLCAO";

    this->kv = &kv_in;

    // Real space Hamiltonian is inited with template TR
    this->hR = new HContainer<TR>(LM_in->ParaV);
    this->sR = new HContainer<TR>(LM_in->ParaV);

    this->getOperator() = new OverlapNew<OperatorLCAO<TK, TR>>(
        LM_in,
        this->kv->kvec_d,
        this->hR,
        &(this->getHk(LM_in)),
        this->sR,
        &(this->getSk(LM_in)),
        &GlobalC::ucell,
        &GlobalC::GridD,
        LM_in->ParaV
    );
}

template<typename TK, typename TR>
HamiltLCAO<TK, TR>::HamiltLCAO(
    Gint_Gamma* GG_in,
    Gint_k* GK_in,
    LCAO_gen_fixedH* genH_in,
    LCAO_Matrix* LM_in,
    Local_Orbital_Charge* loc_in,
    elecstate::Potential* pot_in,
    const K_Vectors& kv_in,
    elecstate::DensityMatrix<TK,double>* DM_in)
{
    this->kv = &kv_in;
    this->classname = "HamiltLCAO";

    // Real space Hamiltonian is inited with template TR
    this->hR = new HContainer<TR>(LM_in->ParaV);
    this->sR = new HContainer<TR>(LM_in->ParaV);

    this->getSk(LM_in).resize(LM_in->ParaV->get_row_size() * LM_in->ParaV->get_col_size());
    this->getHk(LM_in).resize(LM_in->ParaV->get_row_size() * LM_in->ParaV->get_col_size());

    // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>) is registered without template
    std::vector<std::string> pot_register_in;
    if(GlobalV::VL_IN_H)
    {
        if (GlobalV::VION_IN_H)
        {
            pot_register_in.push_back("local");
        }
        if (GlobalV::VH_IN_H)
        {
            pot_register_in.push_back("hartree");
        }
        pot_register_in.push_back("xc");
        if (GlobalV::imp_sol)
        {
            pot_register_in.push_back("surchem");
        }
        if (GlobalV::EFIELD_FLAG)
        {
            pot_register_in.push_back("efield");
        }
        if (GlobalV::GATE_FLAG)
        {
            pot_register_in.push_back("gatefield");
        }
        if (GlobalV::ESOLVER_TYPE == "tddft")
        {
            pot_register_in.push_back("tddft");
        }
    }

    // Gamma_only case to initialize HamiltLCAO
    //
    // code block to construct Operator Chains
    if(std::is_same<TK, double>::value)
    {
        // fix HR to gamma case, where SR will be fixed in Overlap Operator
        this->hR->fix_gamma();
        // initial operator for Gamma_only case
        // overlap term (<psi|psi>) is indispensable
        // in Gamma_only case, target SR is LCAO_Matrix::Sloc, which is same as SK
        this->getOperator() = new OverlapNew<OperatorLCAO<TK, TR>>(
            LM_in,
            this->kv->kvec_d,
            this->hR,
            &(this->getHk(LM_in)),
            this->sR,
            &(this->getSk(LM_in)),
            &GlobalC::ucell,
            &GlobalC::GridD,
            LM_in->ParaV
        );

        // kinetic term (<psi|T|psi>),
        // in Gamma_only case, target HR is LCAO_Matrix::Hloc_fixed, while target HK is LCAO_Matrix::Hloc
        // LCAO_Matrix::Hloc_fixed2 is used for storing
        if(GlobalV::T_IN_H)
        {
            Operator<TK>* ekinetic = new EkineticNew<OperatorLCAO<TK, TR>>(
                LM_in, 
                this->kv->kvec_d, 
                this->hR, 
                &(this->getHk(LM_in)),
                &GlobalC::ucell, 
                &GlobalC::GridD,
                LM_in->ParaV
            );
            this->getOperator()->add(ekinetic);
        }

        // nonlocal term (<psi|beta>D<beta|psi>)
        // in general case, target HR is this->hR, while target HK is LCAO_Matrix::Hloc
        if(GlobalV::VNL_IN_H)
        {
            Operator<TK>* nonlocal = new NonlocalNew<OperatorLCAO<TK, TR>>(
                LM_in, 
                this->kv->kvec_d, 
                this->hR, 
                &(this->getHk(LM_in)),
                &GlobalC::ucell, 
                &GlobalC::GridD,
                LM_in->ParaV
            );
            this->getOperator()->add(nonlocal);
        }

        // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
        // in general case, target HR is Gint::pvpR_grid, while target HK is LCAO_Matrix::Hloc
        if(GlobalV::VL_IN_H)
        {
            //only Potential is not empty, Veff and Meta are available
            if(pot_register_in.size()>0)
            {
                //register Potential by gathered operator
                pot_in->pot_register(pot_register_in);
                //effective potential term
                Operator<TK>* veff = new Veff<OperatorLCAO<TK, TR>>(
                    GG_in,
                    loc_in,
                    LM_in,
                    this->kv->kvec_d,
                    pot_in,
                    this->hR, // no explicit call yet
                    &(this->getHk(LM_in)),
                    &GlobalC::ucell,
                    &GlobalC::GridD,
                    LM_in->ParaV // no explicit call yet
                );
                this->getOperator()->add(veff);
            }
        }

    #ifdef __DEEPKS
        if (GlobalV::deepks_scf)
        {
            Operator<TK>* deepks = new DeePKS<OperatorLCAO<TK, TR>>(loc_in,
                                                                        LM_in,
                                                                        this->kv->kvec_d,
                                                                        this->hR, // no explicit call yet
                                                                        &(this->getHk(LM_in)),
                                                                        &GlobalC::ucell,
                                                                        &GlobalC::GridD,
                                                                        this->kv->nks,
                                                                        DM_in);
            this->getOperator()->add(deepks);
        }
    #endif

        //end node should be OperatorDFTU
        if (GlobalV::dft_plus_u)
        {
            Operator<TK>* dftu = new OperatorDFTU<OperatorLCAO<TK, TR>>(
                LM_in,
                kv->kvec_d,
                this->hR,// no explicit call yet
                &(this->getHk(LM_in)),
                this->kv->isk
            );
            this->getOperator()->add(dftu);
        }
    }
    // multi-k-points case to initialize HamiltLCAO, ops will be used
    else if(std::is_same<TK, std::complex<double>>::value)
    {
        // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
        // Meta potential term (\sum_r <psi(r)|tau(r)|psi(r)>)
        // in general case, target HR is Gint::pvpR_reduced, while target HK is LCAO_Matrix::Hloc2
        if(GlobalV::VL_IN_H)
        {
            //only Potential is not empty, Veff and Meta are available
            if(pot_register_in.size()>0)
            {
                //register Potential by gathered operator
                pot_in->pot_register(pot_register_in);
                //Veff term
                this->getOperator() = new Veff<OperatorLCAO<TK, TR>>(
                    GK_in,
                    loc_in,
                    LM_in,
                    kv->kvec_d,
                    pot_in,
                    this->hR, 
                    &(this->getHk(LM_in)),
                    &GlobalC::ucell,
                    &GlobalC::GridD,
                    LM_in->ParaV
                );
                //reset spin index and real space Hamiltonian matrix
                int start_spin = -1;
                GK_in->reset_spin(start_spin);
                GK_in->destroy_pvpR();
                GK_in->allocate_pvpR();
            }
        }

        // initial operator for multi-k case
        // overlap term is indispensable
        Operator<TK>* overlap = new OverlapNew<OperatorLCAO<TK, TR>>(
            LM_in,
            this->kv->kvec_d,
            this->hR,
            &(this->getHk(LM_in)),
            this->sR,
            &(this->getSk(LM_in)),
            &GlobalC::ucell,
            &GlobalC::GridD,
            LM_in->ParaV
        );
        if(this->getOperator() == nullptr)
        {
            this->getOperator() = overlap;
        }
        else
        {
            this->getOperator()->add(overlap);
        }

        // kinetic term (<psi|T|psi>),
        // in general case, target HR is this->hR, while target HK is LCAO_Matrix::Hloc2
        if(GlobalV::T_IN_H)
        {
            Operator<TK>* ekinetic = new EkineticNew<OperatorLCAO<TK, TR>>(
                LM_in,
                this->kv->kvec_d,
                this->hR,
                &(this->getHk(LM_in)),
                &GlobalC::ucell,
                &GlobalC::GridD,
                LM_in->ParaV
            );
            this->getOperator()->add(ekinetic);
        }

        // nonlocal term (<psi|beta>D<beta|psi>)
        // in general case, target HR is this->hR, while target HK is LCAO_Matrix::Hloc2
        if(GlobalV::VNL_IN_H)
        {
            Operator<TK>* nonlocal = new NonlocalNew<OperatorLCAO<TK, TR>>(
                LM_in,
                this->kv->kvec_d,
                this->hR,
                &(this->getHk(LM_in)),
                &GlobalC::ucell,
                &GlobalC::GridD,
                LM_in->ParaV
            );
            this->getOperator()->add(nonlocal);
        }

    #ifdef __DEEPKS
        if (GlobalV::deepks_scf)
        {
            Operator<TK>* deepks
                = new DeePKS<OperatorLCAO<TK, TR>>(loc_in,
                                                    LM_in,
                                                    this->kv->kvec_d,
                                                    hR,
                                                    &(this->getHk(LM_in)),
                                                    &GlobalC::ucell,
                                                    &GlobalC::GridD,
                                                    this->kv->nks,
                                                    DM_in);
            this->getOperator()->add(deepks);
        }
    #endif
        if (GlobalV::dft_plus_u)
        {
            Operator<TK>* dftu = new OperatorDFTU<OperatorLCAO<TK, TR>>(
                LM_in,
                kv->kvec_d,
                this->hR,// no explicit call yet
                &(this->getHk(LM_in)),
                this->kv->isk
            );
            this->getOperator()->add(dftu);
        }
        if (GlobalV::sc_mag_switch)
        {
            Operator<TK>* sc_lambda = new OperatorScLambda<OperatorLCAO<TK, TR>>(
                LM_in,
                kv->kvec_d,
                this->hR,// no explicit call yet
                &(this->getHk(LM_in)),
                this->kv->isk);
            this->getOperator()->add(sc_lambda);
        }
    }

    ModuleBase::Memory::record("HamiltLCAO::hR", this->hR->get_memory_size());
    ModuleBase::Memory::record("HamiltLCAO::sR", this->sR->get_memory_size());
    
    return;
}

// case for multi-k-points
template <typename TK, typename TR>
void HamiltLCAO<TK, TR>::matrix(MatrixBlock<TK> &hk_in,
                                              MatrixBlock<TK> &sk_in)
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
    //update global spin index
    if (GlobalV::NSPIN == 2)
    {
        // if Veff is added and current_spin is changed, refresh HR
        if(GlobalV::VL_IN_H && this->kv->isk[ik] != GlobalV::CURRENT_SPIN)
        {
            this->refresh();
        }
        GlobalV::CURRENT_SPIN = this->kv->isk[ik];
    }
    this->getOperator()->init(ik);
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
}

template <typename TK, typename TR>
void HamiltLCAO<TK, TR>::refresh()
{
    dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(this->ops)->set_hr_done(false);
}

// get Operator base class pointer
template <typename TK, typename TR>
Operator<TK>*& HamiltLCAO<TK, TR>::getOperator()
{
    return this->ops;
}
// getHk
template <>
std::vector<double>& HamiltLCAO<double, double>::getHk(LCAO_Matrix* LM)
{
    return LM->Hloc;
}

template <>
std::vector<std::complex<double>>& HamiltLCAO<std::complex<double>, double>::getHk(LCAO_Matrix* LM)
{
    return LM->Hloc2;
}
template <>
std::vector<std::complex<double>>& HamiltLCAO<std::complex<double>, std::complex<double>>::getHk(LCAO_Matrix* LM)
{
    return LM->Hloc2;
}

// getSk
template <>
std::vector<double>& HamiltLCAO<double, double>::getSk(LCAO_Matrix* LM)
{
    return LM->Sloc;
}

template <>
std::vector<std::complex<double>>& HamiltLCAO<std::complex<double>, double>::getSk(LCAO_Matrix* LM)
{
    return LM->Sloc2;
}
template <>
std::vector<std::complex<double>>& HamiltLCAO<std::complex<double>, std::complex<double>>::getSk(LCAO_Matrix* LM)
{
    return LM->Sloc2;
}

// getHR
template <typename TK, typename TR>
HContainer<TR>*& HamiltLCAO<TK, TR>::getHR()
{
    return this->hR;
}

// getSR
template <typename TK, typename TR>
HContainer<TR>*& HamiltLCAO<TK, TR>::getSR()
{
    return this->sR;
}

template<typename TK, typename TR>
void HamiltLCAO<TK, TR>::updateSk(const int ik, LCAO_Matrix* LM_in, const int hk_type)
{
    ModuleBase::TITLE("HamiltLCAO", "updateSk");
    ModuleBase::timer::tick("HamiltLCAO", "updateSk");
    ModuleBase::GlobalFunc::ZEROS(this->getSk(LM_in).data(), this->getSk(LM_in).size());
    if(hk_type == 1)// collumn-major matrix for SK
    {
        const int nrow = LM_in->ParaV->get_row_size();
        hamilt::folding_HR(*this->sR, this->getSk(LM_in).data(), this->kv->kvec_d[ik], nrow, 1);
    }
    else if(hk_type == 0) // row-major matrix for SK
    {
        const int ncol = LM_in->ParaV->get_col_size();
        hamilt::folding_HR(*this->sR, this->getSk(LM_in).data(), this->kv->kvec_d[ik], ncol, 0);
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
