#pragma once

#include "../symm_rotation.h"

namespace ModuleSymmetry
{
    // Keeps synthetic symmetry state confined to focused rotation unit tests.
    class SymmetryRotationTestAccess
    {
    public:
        static void set_operation_counts(Symmetry_rotation& rotation, const int nsym, const int nanti)
        {
            rotation.nsym_ = nsym;
            rotation.nanti_ = nanti;
        }

        static std::vector<SpinRotation::Su2>& spin_rotations(Symmetry_rotation& rotation)
        {
            return rotation.spin_U_;
        }

        static std::map<TapR, std::map<int, TapR>>& sector_stars(Symmetry_rotation& rotation)
        {
            return rotation.irs_.sector_stars_;
        }

        template <typename Tdata>
        static RI::Tensor<Tdata> rotate_atompair_abf(const Symmetry_rotation& rotation,
                                                     const RI::Tensor<Tdata>& input,
                                                     const int isym,
                                                     const int type1,
                                                     const int type2)
        {
            return rotation.rotate_atompair_serial_abf(input, isym, type1, type2);
        }
    };
}
