#include "source_lcao/hamilt_lcao_factory.h"

#include "source_base/global_function.h"
#include "source_estate/module_pot/h_tddft_pw.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"

#ifdef __MLALGO
#include "module_operator_lcao/deepks_lcao.h"
#endif

// operator nodes, in construction-chain order:
// overlap -> kinetic -> nonlocal -> veff -> dftu -> tddft
#include "module_operator_lcao/overlap.h"
#include "module_operator_lcao/ekinetic.h"
#include "module_operator_lcao/nonlocal.h"
#include "module_operator_lcao/veff_lcao.h"
#include "module_dftu/dftu_nao_op.h"
#include "module_dftu/dftu_nao_op_legacy.h"
#include "module_operator_lcao/dspin_lcao.h"
#include "module_operator_lcao/td_ekinetic_lcao.h"
#include "module_operator_lcao/td_nonlocal_lcao.h"
#include "module_operator_lcao/td_pot_hybrid.h"

#include <complex>
#include <type_traits>

namespace hamilt
{

namespace
{

/**
 * @brief append the DFT+U operator node shared by gamma and multi-k chains.
 *
 * @param ops current chain head, must be non-null so the node can be added
 */
template <typename TK, typename TR>
void add_dftu_op(Operator<TK>*& ops,
                 const UnitCell& ucell,
                 const Grid_Driver& grid_d,
                 const TwoCenterBundle& two_center_bundle,
                 const LCAO_Orbitals& orb,
                 elecstate::DensityMatrix<TK, double>* DM_in,
                 Plus_U_Base* p_dftu,
                 const Input_para& inp,
                 const K_Vectors* kv,
                 HS_Matrix_K<TK>* hsk,
                 HContainer<TR>* hR)
{
    Operator<TK>* plus_u = nullptr;
    if (inp.dft_plus_u == 1)
    {
        // radius-adjustable localized projections (with onsite_radius)
        plus_u = new DFTU_onsite<OperatorLCAO<TK, TR>>(hsk,
                                                kv->kvec_d, hR,
                                                ucell, &grid_d,
                                                two_center_bundle.overlap_orb_onsite.get(),
                                                orb.cutoffs(), p_dftu,
                                                inp.nspin, inp.onsite_radius, DM_in);
    }
    else if (inp.dft_plus_u == 2)
    {
        // first-zeta NAO projections (old method, kept for testing)
        plus_u = new DFTU_firstzeta<OperatorLCAO<TK, TR>>(hsk,
                                                        kv->kvec_d, hR,
                                                        ucell, p_dftu,
                                                        kv->isk);
    }
    else
    {
        // add_dftu_op is only reachable when dft_plus_u != 0; any other
        // value is an unsupported DFT+U method and must abort here
        ModuleBase::WARNING_QUIT("add_dftu_op",
                                 "unsupported dft_plus_u value: only 1 (radius-adjustable) or 2 (first-zeta NAO) are allowed");
    }
    ops->add(plus_u);
}

#ifdef __MLALGO
/**
 * @brief append the DeePKS operator node shared by gamma and multi-k chains.
 *
 * @param ops current chain head, must be non-null so the node can be added
 * @return the DeePKS V_delta(R) container exposed by the appended node
 */
template <typename TK, typename TR>
HContainer<TR>* add_deepks_op(Operator<TK>*& ops,
                              const UnitCell& ucell,
                              const Grid_Driver& grid_d,
                              const TwoCenterBundle& two_center_bundle,
                              const LCAO_Orbitals& orb,
                              elecstate::DensityMatrix<TK, double>* DM_in,
                              Setup_DeePKS<TK>& deepks,
                              const K_Vectors* kv,
                              HS_Matrix_K<TK>* hsk,
                              HContainer<TR>* hR)
{
    Operator<TK>* deepks_op = new DeePKS<OperatorLCAO<TK, TR>>(hsk,
                                                               kv->kvec_d, hR,
                                                               &ucell, &grid_d,
                                                               two_center_bundle.overlap_orb_alpha.get(),
                                                               &orb, kv->get_nks(),
                                                               DM_in, &deepks.ld);
    ops->add(deepks_op);
    return dynamic_cast<DeePKS<OperatorLCAO<TK, TR>>*>(deepks_op)->get_V_delta_R();
}
#endif

} // anonymous namespace

// build the operator chain for the gamma-only case (TK == double)
template <typename TK, typename TR>
LcaoOpsBundle<TK, TR> build_gamma_ops(const UnitCell& ucell,
                                      const Grid_Driver& grid_d,
                                      const Parallel_Orbitals* paraV,
                                      elecstate::Potential* pot_in,
                                      const TwoCenterBundle& two_center_bundle,
                                      const LCAO_Orbitals& orb,
                                      elecstate::DensityMatrix<TK, double>* DM_in,
                                      Plus_U_Base* p_dftu,
                                      Setup_DeePKS<TK>& deepks,
                                      const Input_para& inp,
                                      const std::vector<std::string>& pot_register_in,
                                      const K_Vectors* kv,
                                      HS_Matrix_K<TK>* hsk,
                                      HContainer<TR>* hR,
                                      HContainer<TR>* sR)
{
    LcaoOpsBundle<TK, TR> bundle;

    // fix HR to gamma case, where SR will be fixed in Overlap Operator
    hR->fix_gamma();
    // initial operator for Gamma_only case
    // overlap term (<psi|psi>) is indispensable
    // in Gamma_only case, target SK is hsk->get_sk(), the target SR is sR
    Operator<TK>* ops = new Overlap<OperatorLCAO<TK, TR>>(hsk,
                                                          kv->kvec_d, hR, sR,
                                                          &ucell, orb.cutoffs(), &grid_d,
                                                          two_center_bundle.overlap_orb.get());

    // kinetic term (<psi|T|psi>)
    if (inp.t_in_h)
    {
        Operator<TK>* ekinetic = new EKinetic<OperatorLCAO<TK, TR>>(hsk,
                                                                    kv->kvec_d, hR,
                                                                    &ucell, orb.cutoffs(), &grid_d,
                                                                    two_center_bundle.kinetic_orb.get());
        ops->add(ekinetic);
    }

    // nonlocal term (<psi|beta>D<beta|psi>)
    // in general case, target HR is hR, while target HK is hsk->get_hk()
    if (inp.vnl_in_h)
    {
        Operator<TK>* nonlocal = new Nonlocal<OperatorLCAO<TK, TR>>(hsk,
                                                                    kv->kvec_d, hR,
                                                                    &ucell, orb.cutoffs(), &grid_d,
                                                                    two_center_bundle.overlap_orb_beta.get());
        ops->add(nonlocal);
    }

    // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
    // in general case, target HR is Gint::hRGint, while target HK is hsk->get_hk()
    if (inp.vl_in_h)
    {
        // only Potential is not empty, Veff and Meta are available
        if (pot_register_in.size() > 0)
        {
            // register Potential by gathered operator
            pot_in->pot_register(pot_register_in);
            // effective potential term
            Operator<TK>* veff = new Veff<OperatorLCAO<TK, TR>>(hsk,
                                                                kv->kvec_d, pot_in,
                                                                hR, // no explicit call yet
                                                                &ucell, orb.cutoffs(), &grid_d,
                                                                inp.nspin);
            ops->add(veff);
        }
    }

#ifdef __MLALGO
    if (inp.deepks_scf)
    {
        bundle.v_delta_R = add_deepks_op<TK, TR>(ops, ucell, grid_d, two_center_bundle,
                                                 orb, DM_in, deepks, kv, hsk, hR);
    }
#endif

    // end node should be DFTU_firstzeta
    if (inp.dft_plus_u)
    {
        add_dftu_op<TK, TR>(ops, ucell, grid_d, two_center_bundle, orb, DM_in,
                            p_dftu, inp, kv, hsk, hR);
    }

    bundle.ops = dynamic_cast<OperatorLCAO<TK, TR>*>(ops);
    return bundle;
}

// build the operator chain for the multi-k case (TK == complex<double>)
template <typename TK, typename TR>
LcaoOpsBundle<TK, TR> build_multik_ops(const UnitCell& ucell,
                                       const Grid_Driver& grid_d,
                                       const Parallel_Orbitals* paraV,
                                       elecstate::Potential* pot_in,
                                       const TwoCenterBundle& two_center_bundle,
                                       const LCAO_Orbitals& orb,
                                       elecstate::DensityMatrix<TK, double>* DM_in,
                                       Plus_U_Base* p_dftu,
                                       Setup_DeePKS<TK>& deepks,
                                       const Input_para& inp,
                                       const std::vector<std::string>& pot_register_in,
                                       const K_Vectors* kv,
                                       HS_Matrix_K<TK>* hsk,
                                       HContainer<TR>* hR,
                                       HContainer<TR>* sR)
{
    LcaoOpsBundle<TK, TR> bundle;

    Operator<TK>* ops = nullptr;

    // Effective potential term (\sum_r <psi(r)|Veff(r)|psi(r)>)
    // Meta potential term (\sum_r <psi(r)|tau(r)|psi(r)>)
    // in general case, target HR is Gint::pvpR_reduced, while target HK is hsk->get_hk()
    if (inp.vl_in_h)
    {
        // only Potential is not empty, Veff and Meta are available
        if (pot_register_in.size() > 0)
        {
            // register Potential by gathered operator
            pot_in->pot_register(pot_register_in);
            // Veff term
            ops = new Veff<OperatorLCAO<TK, TR>>(hsk,
                                                 kv->kvec_d, pot_in,
                                                 hR,
                                                 &ucell, orb.cutoffs(), &grid_d,
                                                 inp.nspin);
        }
    }

    // initial operator for multi-k case
    // overlap term is indispensable
    Operator<TK>* overlap = new Overlap<OperatorLCAO<TK, TR>>(hsk,
                                                              kv->kvec_d, hR, sR,
                                                              &ucell, orb.cutoffs(), &grid_d,
                                                              two_center_bundle.overlap_orb.get());
    if (ops == nullptr)
    {
        ops = overlap;
    }
    else
    {
        ops->add(overlap);
    }

    // kinetic term (<psi|T|psi>),
    // in general case, target HR is hR, while target HK is hsk->get_hk()
    if (inp.t_in_h)
    {
        Operator<TK>* ekinetic = new EKinetic<OperatorLCAO<TK, TR>>(hsk,
                                                                    kv->kvec_d, hR,
                                                                    &ucell, orb.cutoffs(), &grid_d,
                                                                    two_center_bundle.kinetic_orb.get());
        ops->add(ekinetic);
    }

    // nonlocal term (<psi|beta>D<beta|psi>)
    // in general case, target HR is hR, while target HK is hsk->get_hk()
    // TDDFT velocity gauge will calculate full non-local potential including the original one and the
    // correction on its own, so the original non-local potential term should be skipped then
    if (inp.vnl_in_h)
    {
        Operator<TK>* nonlocal = new Nonlocal<OperatorLCAO<TK, TR>>(hsk,
                                                                    kv->kvec_d, hR,
                                                                    &ucell, orb.cutoffs(), &grid_d,
                                                                    two_center_bundle.overlap_orb_beta.get());
        // Nonlocal::initialize_HR inserts atom pairs into hR (using a cutoff
        // that includes the nonlocal pseudopotential radius). TDEkinetic and
        // TDNonlocal both build hR_tmp by iterating over hR's pairs, so the
        // Nonlocal constructor must run even when the operator itself is not
        // added to the chain (TDDFT velocity gauge case).
        if (inp.esolver_type != "tddft" || elecstate::H_TDDFT_pw::stype != 1)
        {
            ops->add(nonlocal);
        }
        else
        {
            delete nonlocal;
        }
    }

#ifdef __MLALGO
    if (inp.deepks_scf)
    {
        bundle.v_delta_R = add_deepks_op<TK, TR>(ops, ucell, grid_d, two_center_bundle,
                                                 orb, DM_in, deepks, kv, hsk, hR);
    }
#endif
    // TDDFT_velocity_gauge
    // These operators are complex-only (no double instantiation of
    // TDEkinetic/TDNonlocal). The std::is_same guard lets the compiler
    // dead-branch-eliminate this block in the double instantiation, avoiding
    // references to missing double symbols.
    if (std::is_same<TK, std::complex<double>>::value && inp.esolver_type == "tddft" && inp.td_stype == 1)
    {
        Operator<TK>* td_ekinetic = new TDEkinetic<OperatorLCAO<TK, TR>>(hsk,
                                                                         hR, kv,
                                                                         &ucell, orb.cutoffs(), &grid_d,
                                                                         two_center_bundle.overlap_orb.get());
        ops->add(td_ekinetic);

        Operator<TK>* td_nonlocal = new TDNonlocal<OperatorLCAO<TK, TR>>(hsk,
                                                                         kv->kvec_d, hR,
                                                                         &ucell, orb, &grid_d);
        ops->add(td_nonlocal);
    }
    if (inp.esolver_type == "tddft" && inp.td_stype == 2)
    {
        Operator<TK>* td_pot_hybrid = new TD_pot_hybrid<OperatorLCAO<TK, TR>>(hsk,
                                                                              kv, hR, sR,
                                                                              orb, &ucell, orb.cutoffs(), &grid_d,
                                                                              two_center_bundle.kinetic_orb.get());
        ops->add(td_pot_hybrid);
    }
    if (inp.dft_plus_u)
    {
        add_dftu_op<TK, TR>(ops, ucell, grid_d, two_center_bundle, orb, DM_in,
                            p_dftu, inp, kv, hsk, hR);
    }
    if (inp.sc_mag_switch)
    {
        Operator<TK>* sc_lambda = new DeltaSpin<OperatorLCAO<TK, TR>>(hsk,
                                                                      kv->kvec_d, hR,
                                                                      ucell, &grid_d,
                                                                      two_center_bundle.overlap_orb_onsite.get(),
                                                                      orb.cutoffs());
        ops->add(sc_lambda);
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.set_operator(sc_lambda);
    }

    bundle.ops = dynamic_cast<OperatorLCAO<TK, TR>*>(ops);
    return bundle;
}

// explicit instantiation: gamma-only, multi-k, and non-collinear spin cases
template struct LcaoOpsBundle<double, double>;
template struct LcaoOpsBundle<std::complex<double>, double>;
template struct LcaoOpsBundle<std::complex<double>, std::complex<double>>;

template LcaoOpsBundle<double, double> build_gamma_ops<double, double>(
    const UnitCell&, const Grid_Driver&, const Parallel_Orbitals*,
    elecstate::Potential*, const TwoCenterBundle&, const LCAO_Orbitals&,
    elecstate::DensityMatrix<double, double>*, Plus_U_Base*, Setup_DeePKS<double>&,
    const Input_para&, const std::vector<std::string>&, const K_Vectors*,
    HS_Matrix_K<double>*, HContainer<double>*, HContainer<double>*);

template LcaoOpsBundle<double, double> build_multik_ops<double, double>(
    const UnitCell&, const Grid_Driver&, const Parallel_Orbitals*,
    elecstate::Potential*, const TwoCenterBundle&, const LCAO_Orbitals&,
    elecstate::DensityMatrix<double, double>*, Plus_U_Base*, Setup_DeePKS<double>&,
    const Input_para&, const std::vector<std::string>&, const K_Vectors*,
    HS_Matrix_K<double>*, HContainer<double>*, HContainer<double>*);

template LcaoOpsBundle<std::complex<double>, double> build_gamma_ops<std::complex<double>, double>(
    const UnitCell&, const Grid_Driver&, const Parallel_Orbitals*,
    elecstate::Potential*, const TwoCenterBundle&, const LCAO_Orbitals&,
    elecstate::DensityMatrix<std::complex<double>, double>*, Plus_U_Base*,
    Setup_DeePKS<std::complex<double>>&, const Input_para&,
    const std::vector<std::string>&, const K_Vectors*,
    HS_Matrix_K<std::complex<double>>*, HContainer<double>*, HContainer<double>*);

template LcaoOpsBundle<std::complex<double>, double> build_multik_ops<std::complex<double>, double>(
    const UnitCell&, const Grid_Driver&, const Parallel_Orbitals*,
    elecstate::Potential*, const TwoCenterBundle&, const LCAO_Orbitals&,
    elecstate::DensityMatrix<std::complex<double>, double>*, Plus_U_Base*,
    Setup_DeePKS<std::complex<double>>&, const Input_para&,
    const std::vector<std::string>&, const K_Vectors*,
    HS_Matrix_K<std::complex<double>>*, HContainer<double>*, HContainer<double>*);

template LcaoOpsBundle<std::complex<double>, std::complex<double>>
build_gamma_ops<std::complex<double>, std::complex<double>>(
    const UnitCell&, const Grid_Driver&, const Parallel_Orbitals*,
    elecstate::Potential*, const TwoCenterBundle&, const LCAO_Orbitals&,
    elecstate::DensityMatrix<std::complex<double>, double>*, Plus_U_Base*,
    Setup_DeePKS<std::complex<double>>&, const Input_para&,
    const std::vector<std::string>&, const K_Vectors*,
    HS_Matrix_K<std::complex<double>>*, HContainer<std::complex<double>>*,
    HContainer<std::complex<double>>*);

template LcaoOpsBundle<std::complex<double>, std::complex<double>>
build_multik_ops<std::complex<double>, std::complex<double>>(
    const UnitCell&, const Grid_Driver&, const Parallel_Orbitals*,
    elecstate::Potential*, const TwoCenterBundle&, const LCAO_Orbitals&,
    elecstate::DensityMatrix<std::complex<double>, double>*, Plus_U_Base*,
    Setup_DeePKS<std::complex<double>>&, const Input_para&,
    const std::vector<std::string>&, const K_Vectors*,
    HS_Matrix_K<std::complex<double>>*, HContainer<std::complex<double>>*,
    HContainer<std::complex<double>>*);

} // namespace hamilt
