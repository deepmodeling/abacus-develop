#include "source_lcao/module_rt/boundary_fix.h"

#include "source_base/basis_flattened_indexing.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/parallel_common.h"
#include "source_cell/unitcell.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace
{
using Complex = std::complex<double>;

ModuleBase::Element_Basis_Index::IndexLNM single_orbital_hierarchy()
{
    ModuleBase::Element_Basis_Index::Range range(1);
    range[0].resize(1);
    range[0][0].N = 1;
    range[0][0].M = 1;
    return ModuleBase::Element_Basis_Index::construct_index(range);
}
} // namespace

TEST(BoundaryFixTest, ApplyAtomPhaseFromFlattenedIndexing)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = 2;
    ucell.atoms[0].nw = 1;
    ucell.atoms[0].boundary_shift.resize(2, {0, 0, 0});
    ucell.atoms[0].boundary_shift[1] = {1, 0, 0};
    ucell.set_iat2itia();
    ucell.set_iat2iwt(1);

    const std::vector<int> atom_types = {0, 0};
    const BasisFlattenedIndexing basis(single_orbital_hierarchy(), atom_types, 1);
    Parallel_Blkcyclic layout;
#ifdef __MPI
    const int process_count = Parallel_Common::comm_size(MPI_COMM_WORLD);
    if (process_count == 1)
    {
        layout.set_serial(2, 2);
    }
    else
    {
        layout.init(2, 2, 1, MPI_COMM_WORLD, false);
    }
#else
    layout.set_serial(2, 2);
#endif
    SCOPED_TRACE(testing::Message() << "process (" << layout.irproc() << ", " << layout.icproc()
                                    << "), grid (" << layout.nrproc() << ", " << layout.ncproc()
                                    << "), local rows " << layout.nrloc() << ", local columns "
                                    << layout.ncloc());

    const std::vector<int> ngk = {layout.nrloc()};
    psi::Psi<Complex> psi_prev(1, 1, layout.nrloc(), ngk, true);
    psi_prev.fix_k(0);
    std::fill(psi_prev.get_pointer(), psi_prev.get_pointer() + layout.nrloc(), Complex(1.0, 0.0));

    const int64_t local_size = static_cast<int64_t>(layout.size());
    ct::Tensor hamilt_prev(ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({local_size}));
    ct::Tensor ovlp_prev(ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({local_size}));
    std::fill(hamilt_prev.data<Complex>(), hamilt_prev.data<Complex>() + layout.size(), Complex(1.0, 0.0));
    std::fill(ovlp_prev.data<Complex>(), ovlp_prev.data<Complex>() + layout.size(), Complex(1.0, 0.0));

    const std::vector<ModuleBase::Vector3<double>> kvec_d = {{0.25, 0.0, 0.0}};
    const ModuleBase::Vector3<double> cart_At_prev = {0.0, 0.0, 0.0};
    module_rt::reset_matrix_boundary(ucell,
                                     kvec_d,
                                     cart_At_prev,
                                     basis,
                                     layout,
                                     hamilt_prev,
                                     ovlp_prev,
                                     &psi_prev);

    const Complex* hamilt = hamilt_prev.data<Complex>();
    const Complex* ovlp = ovlp_prev.data<Complex>();
    const Complex expected_matrix[2][2] = {{Complex(1.0, 0.0), Complex(0.0, -1.0)},
                                           {Complex(0.0, 1.0), Complex(1.0, 0.0)}};
    for (int icglb = 0; icglb < 2; ++icglb)
    {
        const int icloc = layout.icloc(icglb);
        if (icloc < 0)
        {
            continue;
        }
        for (int irglb = 0; irglb < 2; ++irglb)
        {
            const int irloc = layout.irloc(irglb);
            if (irloc < 0)
            {
                continue;
            }
            const int index = irloc + icloc * layout.nrloc();
            SCOPED_TRACE(testing::Message() << "global row " << irglb << ", global column "
                                            << icglb << ", local index " << index);
            EXPECT_NEAR(hamilt[index].real(), expected_matrix[irglb][icglb].real(), 1e-12);
            EXPECT_NEAR(hamilt[index].imag(), expected_matrix[irglb][icglb].imag(), 1e-12);
            EXPECT_NEAR(ovlp[index].real(), hamilt[index].real(), 1e-12);
            EXPECT_NEAR(ovlp[index].imag(), hamilt[index].imag(), 1e-12);
        }
    }

    psi_prev.fix_k(0);
    const Complex expected_psi[2] = {Complex(1.0, 0.0), Complex(0.0, 1.0)};
    for (int irglb = 0; irglb < 2; ++irglb)
    {
        const int irloc = layout.irloc(irglb);
        if (irloc >= 0)
        {
            EXPECT_NEAR(psi_prev.get_pointer()[irloc].real(), expected_psi[irglb].real(), 1e-12);
            EXPECT_NEAR(psi_prev.get_pointer()[irloc].imag(), expected_psi[irglb].imag(), 1e-12);
        }
    }
}

TEST(BoundaryFixTest, NoSavedStateIsNoOp)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nw = 1;
    ucell.atoms[0].boundary_shift.resize(1, {0, 0, 0});
    ucell.set_iat2itia();
    ucell.set_iat2iwt(1);

    const std::vector<int> atom_types = {0};
    const BasisFlattenedIndexing basis(single_orbital_hierarchy(), atom_types, 1);
    Parallel_Blkcyclic layout;
    layout.set_serial(1, 1);

    ct::Tensor hamilt_last(ct::DataType::DT_COMPLEX_DOUBLE);
    ct::Tensor ovlp_last(ct::DataType::DT_COMPLEX_DOUBLE);
    const std::vector<ModuleBase::Vector3<double>> kvec_d = {{0.0, 0.0, 0.0}};
    const ModuleBase::Vector3<double> cart_At_last = {0.0, 0.0, 0.0};

    EXPECT_NO_FATAL_FAILURE(module_rt::reset_matrix_boundary(ucell,
                                                              kvec_d,
                                                              cart_At_last,
                                                              basis,
                                                              layout,
                                                              hamilt_last,
                                                              ovlp_last,
                                                              nullptr));
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
