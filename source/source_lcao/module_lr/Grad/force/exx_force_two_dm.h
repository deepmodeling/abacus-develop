#pragma once
#ifdef __EXX
#include <RI/physics/Exx.h>

#include <array>
#include <map>
#include <string>
#include <utility>

namespace LR
{
    /// @brief `RI::Exx` plus a `cal_force` overload taking two *different* density matrices.
    ///
    /// LibRI used to ship this as `RI::LR` (a subclass of `RI::Exx` in
    /// `RI/physics/LR.h`). LibRI commit 5c6c262 repurposed the name `RI::LR` for an
    /// unrelated k-space CVCX/Hartree helper built on `LRI_k`, so the piece the
    /// LR-TDDFT analytical gradient needs is kept here instead.
    ///
    /// `RI::Exx::cal_force` contracts the 3-center derivative dH with whatever sits in
    /// `post_2D.saves["Ds_"+suffix]`, while the density matrix handed to `set_Ds` is the
    /// one consumed inside the loop-3 contraction. Overwriting the former therefore lets
    /// the two sides of Tr[D_IJ * dH[D_KL]] differ, which is what the
    /// Pulay / Hellmann-Feynman split of the EXX gradient requires.
    template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
    class ExxForceTwoDM : public RI::Exx<TA, Tcell, Ndim, Tdata>
    {
        using Base = RI::Exx<TA, Tcell, Ndim, Tdata>;

    public:
        using TC = std::array<Tcell, Ndim>;
        using TAC = std::pair<TA, TC>;

        ExxForceTwoDM() = default;
        /// take over an existing Exx kernel (Cs/Vs/dCs/dVs already set up)
        ExxForceTwoDM(Base&& exx) : Base(std::move(exx)) {}

        using Base::cal_force;

        /// @param Ds_left  the D_IJ contracted with dH; if empty, behaves as `Exx::cal_force`
        /// @param save_names_suffix  "Cs", "Vs", "Ds", "dCs", "dVs"
        void cal_force(const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Ds_left,
            const std::array<std::string, 5>& save_names_suffix = { "","","","","" })
        {
            if (!Ds_left.empty())
            {
                this->post_2D.saves["Ds_" + save_names_suffix[2]]
                    = this->post_2D.set_tensors_map2(Ds_left);
            }
            this->cal_force(save_names_suffix);
        }
    };
}
#endif
