#include "hamilt_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_hsolver/diago_elpa.h"
#include "src_lcao/dftu.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"
#ifdef __DEEPKS
#include "module_deepks/LCAO_deepks.h"
#include "ks_lcao/deepks_lcao.h"
#endif
#include "ks_lcao/op_dftu_lcao.h"
#include "ks_lcao/ekinetic_lcao.h"
#include "ks_lcao/meta_lcao.h"
#include "ks_lcao/nonlocal_lcao.h"
#include "ks_lcao/op_exx_lcao.h"
#include "ks_lcao/overlap_lcao.h"
#include "ks_lcao/veff_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_xc/xc_functional.h"

namespace hamilt
{
// case for nspin<4, gamma-only k-point
template class HamiltLCAO<double>;
// case for nspin<4, multi-k-points
// case for nspin == 4, non-collinear spin case
template class HamiltLCAO<std::complex<double>>;

template<typename T>
HamiltLCAO<T>::HamiltLCAO(
    Gint_Gamma* GG_in,
    LCAO_gen_fixedH* genH_in,
    LCAO_Matrix* LM_in,
    Local_Orbital_Charge* loc_in)
{
    this->classname = "HamiltLCAO";
    //reset fixed Hamiltonian matrix in real space
    LM_in->zeros_HSgamma('T');
    //reset Overlap matrix in real space
    LM_in->zeros_HSgamma('S');

    // initial operator for Gamma_only case
    // overlap term (<psi|psi>) is indispensable
    // in Gamma_only case, target SR is LCAO_Matrix::Sloc, which is same as SK
    this->opsd = new Overlap<OperatorLCAO<double>>(
        genH_in,
        LM_in,
        &(LM_in->Sloc),
        &(LM_in->Sloc)
    );

    // kinetic term (<psi|T|psi>), 
    // in Gamma_only case, target HR is LCAO_Matrix::Hloc_fixed, while target HK is LCAO_Matrix::Hloc
    // LCAO_Matrix::Hloc_fixed2 is used for storing 
    if(GlobalV::T_IN_H)
    {
        Operator<double>* ekinetic = new Ekinetic<OperatorLCAO<double>>(
            genH_in,
            LM_in,
            &(LM_in->Hloc_fixed),
            &(LM_in->Hloc)
        );
        this->opsd->add(ekinetic);
    }

    // nonlocal term (<psi|beta>D<beta|psi>)
    // in general case, target HR is LCAO_Matrix::Hloc_fixedR, while target HK is LCAO_Matrix::Hloc
    if(GlobalV::VNL_IN_H)
    {
        Operator<double>* nonlocal = new Nonlocal<OperatorLCAO<double>>(
            genH_in,
            LM_in,
            &(LM_in->Hloc_fixed),
            &(LM_in->Hloc)
        );
        this->opsd->add(nonlocal);
    }

    // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
    // in general case, target HR is Gint::pvpR_grid, while target HK is LCAO_Matrix::Hloc
    if(GlobalV::VL_IN_H)
    {
        //effective potential term
        Operator<double>* veff = new Veff<OperatorLCAO<double>>(
            GG_in,
            loc_in,
            LM_in,
            nullptr, // no explicit call yet
            &(LM_in->Hloc) // no explicit call yet
        );
        this->opsd->add(veff);

        // Meta term
        Operator<double>* meta = new Meta<OperatorLCAO<double>>(
            GG_in,
            loc_in,
            LM_in,
            nullptr, // no explicit call yet
            &(LM_in->Hloc) // no explicit call yet
        );
        this->opsd->add(meta);

        //exact exchange term
        Operator<double>* exx = new OperatorEXX<OperatorLCAO<double>>(
            LM_in,
            nullptr, //no explicit call yet
            &(LM_in->Hloc)
        );
        this->opsd->add(exx);
    }

    if (GlobalV::dft_plus_u)
    {
        Operator<double>* dftu = new OperatorDFTU<OperatorLCAO<double>>(
            LM_in,
            nullptr,// no explicit call yet
            &(LM_in->Hloc)
        );
        this->opsd->add(dftu);
    }

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        Operator<double>* deepks = new DeePKS<OperatorLCAO<double>>(
            loc_in,
            LM_in,
            nullptr,// no explicit call yet
            &(LM_in->Hloc)
        );
        this->opsd->add(deepks);
    }
#endif
}

template<typename T>
HamiltLCAO<T>::HamiltLCAO(
    Gint_k* GK_in,
    LCAO_gen_fixedH* genH_in,
    LCAO_Matrix* LM_in,
    Local_Orbital_Charge* loc_in)
{
    this->classname = "HamiltLCAO";
    
    //reset fixed Hamiltonian matrix in real space
    LM_in->zeros_HSR('T');
    //reset Overlap matrix in real space
    LM_in->zeros_HSR('S');

    // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
    // Meta potential term (\sum_r <psi(r)|tau(r)|psi(r)>)
    // in general case, target HR is Gint::pvpR_reduced, while target HK is LCAO_Matrix::Hloc2
    if(GlobalV::VL_IN_H)
    {
        //Veff term
        this->ops = new Veff<OperatorLCAO<std::complex<double>>>(
            GK_in,
            loc_in,
            LM_in,
            nullptr, // no explicit call yet
            &(LM_in->Hloc2) // no explicit call yet
        );
        //reset spin index and real space Hamiltonian matrix
        int start_spin = -1;
        GK_in->reset_spin(start_spin);
        GK_in->destroy_pvpR();
        GK_in->allocate_pvpR();

        // Meta term
        Operator<std::complex<double>>* meta = new Meta<OperatorLCAO<std::complex<double>>>(
            GK_in,
            loc_in,
            LM_in,
            nullptr, // no explicit call yet
            &(LM_in->Hloc2) // no explicit call yet
        );
        this->ops->add(meta);

        //exact exchange term
        Operator<std::complex<double>>* exx = new OperatorEXX<OperatorLCAO<std::complex<double>>>(
            LM_in,
            nullptr, //no explicit call yet
            &(LM_in->Hloc2)
        );
        this->ops->add(exx);
    }

    // initial operator for multi-k case
    // overlap term is indispensable
    Operator<std::complex<double>>* overlap = new Overlap<OperatorLCAO<std::complex<double>>>(
        genH_in,
        LM_in,
        &(LM_in->SlocR),
        &(LM_in->Sloc2)
    );
    this->ops->add(overlap);

    // kinetic term (<psi|T|psi>), 
    // in general case, target HR is LCAO_Matrix::Hloc_fixedR, while target HK is LCAO_Matrix::Hloc2
    if(GlobalV::T_IN_H)
    {
        Operator<std::complex<double>>* ekinetic = new Ekinetic<OperatorLCAO<std::complex<double>>>(
            genH_in,
            LM_in,
            &(LM_in->Hloc_fixedR),
            &(LM_in->Hloc2)
        );
        this->ops->add(ekinetic);
    }

    // nonlocal term (<psi|beta>D<beta|psi>)
    // in general case, target HR is LCAO_Matrix::Hloc_fixedR, while target HK is LCAO_Matrix::Hloc2
    if(GlobalV::VNL_IN_H)
    {
        Operator<std::complex<double>>* nonlocal = new Nonlocal<OperatorLCAO<std::complex<double>>>(
            genH_in,
            LM_in,
            &(LM_in->Hloc_fixedR),
            &(LM_in->Hloc2)
        );
        this->ops->add(nonlocal);
    }

    if (GlobalV::dft_plus_u)
    {
        Operator<std::complex<double>>* dftu = new OperatorDFTU<OperatorLCAO<std::complex<double>>>(
            LM_in,
            nullptr,// no explicit call yet
            &(LM_in->Hloc2)
        );
        this->ops->add(dftu);
    }

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        Operator<std::complex<double>>* deepks = new DeePKS<OperatorLCAO<std::complex<double>>>(
            loc_in,
            LM_in,
            nullptr,// no explicit call yet
            &(LM_in->Hloc2)
        );
        this->ops->add(deepks);
    }
#endif

}

// case for multi-k-points
template <>
void HamiltLCAO<std::complex<double>>::matrix(MatrixBlock<std::complex<double>> &hk_in,
                                              MatrixBlock<std::complex<double>> &sk_in)
{
    auto op = dynamic_cast<OperatorLCAO<std::complex<double>>*>(this->ops);
    assert(op != NULL);
    op->matrixHk(hk_in, sk_in);
}

// case for nspin<4, gamma_only
template <> void HamiltLCAO<double>::matrix(MatrixBlock<double> &hk_in, MatrixBlock<double> &sk_in)
{
    auto op = dynamic_cast<OperatorLCAO<double>*>(this->opsd);
    assert(op != NULL);
    op->matrixHk(hk_in, sk_in);
}

template <> void HamiltLCAO<double>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltLCAO", "updateHk");
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
    //update global spin index
    if (GlobalV::NSPIN == 2)
    {
        GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
    }
    this->opsd->init(ik);
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
}

template <> void HamiltLCAO<std::complex<double>>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltLCAO", "updateHk");
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
    //update global spin index
    if (GlobalV::NSPIN == 2)
    {
        GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
    }
    this->ops->init(ik);
    ModuleBase::timer::tick("HamiltLCAO", "updateHk");
}

} // namespace hamilt
