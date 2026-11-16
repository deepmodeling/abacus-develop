#include "read_input.h"
#include "read_input_tool.h"
#include "source_base/tool_quit.h"

#include <cmath>
#include <set>
#include <utility>

namespace ModuleIO
{
void ReadInput::item_lr_nacv()
{
    {
        Input_Item item("lr_nac_pairs");
        item.annotation = "selected analytic auxiliary-TDA derivative-coupling surface pairs";
        item.category = "Linear Response TDDFT";
        item.type = "Vector of Integer";
        item.default_value = "";
        item.unit = "";
        item.set_availability(
            "esolver_type==ks-lr and basis_type==lcao and gamma_only==true and cal_force==true and dft_plus_u==0");
        item.description =
            R"(An even list of surface indices, for example `0 1 1 2`, requests S0-S1 and S1-S2 analytic derivative couplings after each LR geometry. Zero is the KS reference determinant; positive indices are LR roots counted from one within `lr_nac_spin`. Each pair must satisfy 0 <= bra < ket <= lr_nstates and may appear only once. The empty default disables this output and adds no response solves.

Requires a converged real Gamma `ks-lr` calculation, `cal_force=true`, `kpar=1`, `dft_plus_u=0`, and an LR eigensolver other than spectrum. Occupied bands must be fully occupied and virtual bands empty within floating-point roundoff; fractional-occupation response is not implemented. Hubbard response and projector derivatives are not implemented in the analytic LR force/NAC path. The implementation differentiates the spin-adapted or collinear unrestricted auxiliary TDA determinant wavefunction defined by the LR orbitals and amplitudes; it is not full quadratic-response TDDFT. Use `lr_nac_spin=updown` for an unrestricted reference with nonempty occupied and virtual windows in both spin channels. Closed-shell singlet/triplet adaptations use the corresponding spin label; their S0-triplet spin-free couplings are exactly zero without SOC. Complex-k and spin-flip calculations are rejected.

The running log contains NACV_BEGIN/COMPONENT/END records: raw, ETF and moving-AO vectors in bohr^-1, then raw/ETF gap-weighted numerators in Hartree/bohr. The energy gap is printed in Hartree. An unresolved gap retains its numerators, sets the resolution flag false and leaves the divided vectors zero; no energy-gap clipping is used. Phase and state ordering follow the current LR eigenvectors. Compare degenerate states as subspaces, not individual vector components.

The occupied and virtual LR windows must include every occupied and virtual orbital of the basis. Truncated windows are rejected because padding amplitudes with zeros does not restore omitted orbital-response terms. Both the resolved KS functional and the LR kernel must provide the implemented analytic response derivatives; a KS-supported name alone is insufficient. Resolved functional names are checked at the NAC entry after the KS reference is available. Built-in KS aliases are resolved for ground-state response. Verify the chosen functional and AO basis against the same-definition overlap finite difference before quantitative dynamics.)";
        item.read_value = [](const Input_Item &item, Parameter &para)
        {
            std::vector<int> pairs;
            for (const std::string &value : item.str_values)
            {
                std::size_t consumed = 0;
                int index = 0;
                try
                {
                    index = std::stoi(value, &consumed);
                }
                catch (const std::exception &)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_pairs requires integer surface indices");
                }
                if (consumed != value.size())
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_pairs requires integer surface indices");
                }
                pairs.push_back(index);
            }
            para.input.lr_nac_pairs = std::move(pairs);
        };
        item.check_value = [](const Input_Item &, const Parameter &para)
        {
            const Input_para &input = para.input;
            const std::vector<int> &pairs = input.lr_nac_pairs;
            if (pairs.empty())
            {
                return;
            }
            if (input.dft_plus_u != 0)
            {
                ModuleBase::WARNING_QUIT(
                    "ReadInput", "lr_nac_pairs requires dft_plus_u=0; Hubbard response and projector derivatives are "
                                 "not implemented");
            }
            if (pairs.size() % 2 != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_pairs requires an even number of surface indices");
            }
            std::set<std::pair<int, int>> seen;
            for (std::size_t index = 0; index < pairs.size(); index += 2)
            {
                const int bra = pairs[index];
                const int ket = pairs[index + 1];
                if (bra < 0 || ket <= bra || (input.lr_nstates > 0 && ket > input.lr_nstates))
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_pairs requires 0 <= bra < ket <= lr_nstates");
                }
                if (!seen.insert(std::make_pair(bra, ket)).second)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_pairs contains a duplicate surface pair");
                }
            }
            if (input.esolver_type != "ks-lr" || input.basis_type != "lcao" || !input.gamma_only || !input.cal_force ||
                input.kpar != 1 || para.sys.kpar_lcao != 1 || (input.lr_unrestricted && input.nspin != 2) ||
                (input.nspin != 1 && input.nspin != 2) || input.lr_solver == "spectrum" ||
                input.calculation == "cell-relax" || input.ri_hartree_benchmark != "none")
            {
                ModuleBase::WARNING_QUIT(
                    "ReadInput", "lr_nac_pairs requires real Gamma ks-lr, cal_force=true, kpar=1 and solved LR states");
            }
            if (input.lr_nac_spin != "singlet" && input.lr_nac_spin != "triplet" && input.lr_nac_spin != "updown")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_spin must be singlet, triplet or updown");
            }
            if ((input.lr_nac_spin == "triplet" || input.lr_nac_spin == "updown") && input.nspin != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "lr_nac_spin=triplet/updown requires nspin=2");
            }
        };
        sync_intvec(input.lr_nac_pairs, para.input.lr_nac_pairs.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("lr_nac_spin");
        item.domain.set_choices({"singlet", "triplet", "updown"});
        item.annotation = "spin adaptation or unrestricted channel for selected analytic NAC pairs";
        item.category = "Linear Response TDDFT";
        item.type = "String";
        item.default_value = "singlet";
        item.unit = "";
        item.set_availability(
            "esolver_type==ks-lr and basis_type==lcao and gamma_only==true and cal_force==true and dft_plus_u==0");
        item.description =
            "Select singlet/triplet for a closed-shell reference or updown for a collinear unrestricted "
            "reference with two nonempty occupied spin spaces. Triplet and updown require nspin=2. "
            "S0-triplet spin-free couplings of a closed-shell reference are zero; updown keeps the "
            "reference M_S and does not enforce a pure total-spin eigenstate. This does not enable SOC or ISC. "
            "The default is singlet and an empty lr_nac_pairs list disables NAC output.";
        read_sync_string(input.lr_nac_spin);
        this->add_item(item);
    }
}

void ReadInput::item_lr_nacv_response()
{
    {
        Input_Item item("lr_grad_thr");
        item.annotation = "absolute true-residual norm tolerance for Z-vector CG";
        item.set_availability("lr_grad_solver==cg");
        item.category = "Linear Response TDDFT";
        item.type = "Real";
        item.default_value = "1e-6";
        item.unit = "Ry";
        item.description =
            R"(Positive finite absolute global Euclidean norm bound on $(A+B)Z-R$ for `lr_grad_solver cg`. The norm combines all distributed entries and all right-hand sides in a solve. It is recomputed with the actual response operator after iteration, rather than accepted from the recursive CG residual alone. Applies to LR analytic forces and auxiliary TDA ground-excited and excited-excited NAC adjoints. Ignored by direct Z-vector solvers. It is independent of `lr_thr`, which controls the excitation eigensolver. Small excitation gaps can amplify response errors in divided NAC vectors; check convergence of the requested observable explicitly.)";
        item.check_value = [](const Input_Item &item, const Parameter &para)
        {
            const double tolerance = para.input.lr_grad_thr;
            if (!std::isfinite(tolerance) || tolerance <= 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "lr_grad_thr must be positive and finite");
            }
        };
        read_sync_double(input.lr_grad_thr);
        this->add_item(item);
    }
    {
        Input_Item item("lr_grad_maxiter");
        item.annotation = "maximum number of NAC CG updates";
        item.set_availability("lr_grad_solver==cg");
        item.category = "Linear Response TDDFT";
        item.type = "Integer";
        item.default_value = "100";
        item.unit = "";
        item.description =
            R"(Positive maximum total number of updates for `lr_grad_solver cg` in auxiliary TDA NAC response solves, including corrections restarted from the true residual. Ignored by direct Z-vector solvers. Reaching the limit succeeds only if the explicitly recomputed residual satisfies `lr_grad_thr`; otherwise the calculation fails rather than returning an uncertified NAC.)";
        item.check_value = [](const Input_Item &item, const Parameter &para)
        {
            if (para.input.lr_grad_maxiter < 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "lr_grad_maxiter must be positive");
            }
        };
        read_sync_int(input.lr_grad_maxiter);
        this->add_item(item);
    }
}

} // namespace ModuleIO
