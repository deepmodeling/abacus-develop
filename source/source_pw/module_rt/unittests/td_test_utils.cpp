#include "source_base/parallel_global.h"
#include "source_cell/read_pp.h"
#include "source_pw/module_rt/unittests/td_test_utils.h"

#include <cmath>
#include <utility>

static_assert(noexcept(std::swap(std::declval<Input_para&>(), std::declval<Input_para&>())),
              "Test input installation and restoration must not throw.");
static_assert(noexcept(std::swap(std::declval<System_para&>(), std::declval<System_para&>())),
              "Test system installation and restoration must not throw.");

TestParameters::TestParameters(Parameter& target) : target_(target)
{
    input_.nspin = 4;
    input_.noncolin = true;
    input_.lspinorb = true;
    input_.device = "cpu";
    input_.basis_type = "pw";
    input_.precision = "double";
    system_.nqx = 600;
    system_.nqxq = 600;
    system_.has_double_data = true;
    system_.has_float_data = false;
    // Prepare the complete configuration before installing it; retain the original state in these members.
    std::swap(target_.input, input_);
    std::swap(target_.sys, system_);
}

TestParameters::~TestParameters() noexcept
{
    std::swap(target_.input, input_);
    std::swap(target_.sys, system_);
}

namespace td_test
{
Observables::Observables() : parameters_(PARAM)
{
}

void Observables::SetUp()
{
    cell_.ntype = 1;
    cell_.nat = 2;
    cell_.atoms = new Atom[1];
    cell_.set_atom_flag = true;
    cell_.atoms[0].na = 2;
    cell_.atoms[0].label = "Fe";
    cell_.atoms[0].taud = {{0.11, 0.21, 0.31}, {0.59, 0.43, 0.77}};
    cell_.lat0 = 6.0;
    cell_.latvec = ModuleBase::Matrix3(1, 0, 0, 0.4, 0.8, 0, 0.1, 0.2, 1);
    cell_.G = cell_.latvec.Inverse().Transpose();
    cell_.omega = std::pow(cell_.lat0, 3) * cell_.latvec.Det();
    cell_.tpiba = ModuleBase::TWO_PI / cell_.lat0;
    cell_.tpiba2 = cell_.tpiba * cell_.tpiba;
    for (const ModuleBase::Vector3<double>& position: cell_.atoms[0].taud)
    {
        cell_.atoms[0].tau.push_back(position * cell_.latvec);
    }
    cell_.set_iat2itia();
    cell_.set_iat2iwt(2);

    Pseudopot_upf reader;
    std::string type = "upf201";
    const std::string filename = std::string(TD_TEST_PP_DIR) + "/Fe.upf";
    ASSERT_EQ(reader.init_pseudo_reader(filename, type, cell_.atoms[0].ncpp), 0);
    Atom_pseudo& atom_pp = cell_.atoms[0].ncpp;
    atom_pp.nh = 0;
    for (const int l: atom_pp.lll)
    {
        atom_pp.nh += 2 * l + 1;
    }

#ifdef __MPI
    basis_.initmpi(1, 0, POOL_WORLD);
#endif
    basis_.initgrids(cell_.lat0, cell_.latvec, 12.0);
    const ModuleBase::Vector3<double> k(0.17, 0.11, 0.07);
    basis_.initparameters(false, 3.0, 1, &k);
    basis_.setuptransform();
    basis_.collect_local_pw();
    kv_.set_nkstot(1);
    kv_.set_nks(1);
    kv_.isk = {0};
    kv_.ik2iktot = {0};
    sf_.set(&basis_, 0);
    Parallel_Grid grid;
    sf_.setup(&cell_, grid, &basis_, false);
    pp_.init(cell_, &sf_, &basis_, true);
    pp_.init_vnl(cell_, &basis_);
    charge_.set_rhopw(&basis_);
    charge_.allocate(4, false, false, 0);
    state_.charge = &charge_;
    state_.pot = new elecstate::Potential;
    state_.pot->get_eff_v().create(4, basis_.nrxx);
    pp_.cal_effective_D(state_.pot->get_eff_v(), &basis_, cell_);
    state_.wg.create(1, 2);
    state_.wg(0, 0) = 0.8;
    state_.wg(0, 1) = 0.3;
    const int n = 2 * basis_.npwk_max;
    const std::vector<int> sizes{basis_.npwk[0]};
    psi_.reset(new psi::Psi<Complex>(1, 2, n, sizes, true));
    for (int band = 0; band < 2; ++band)
    {
        double norm = 0;
        for (int i = 0; i < n; ++i)
        {
            const Complex value(std::sin(0.37 * (i + 1) * (band + 1)), std::cos(0.23 * (i + 2) * (band + 2)));
            (*psi_)(band, i) = value;
            norm += std::norm(value);
        }
        for (int i = 0; i < n; ++i)
        {
            (*psi_)(band, i) /= std::sqrt(norm);
        }
    }
}

void Observables::TearDown()
{
    pp_.release_memory();
}

Vector Observables::projectors(const ModuleBase::Vector3<double>& a)
{
    Vector beta(pp_.nkb * basis_.npwk_max);
    pp_.getvnl<double, Device>(nullptr, cell_, 0, a, beta.data());
    return beta;
}

Vector Observables::nonlocal(const Vector& beta, const Vector& x, int atom) const
{
    // Explicit beta D beta-adjoint action in Hartree units, independent of the production contraction kernels.
    const int npw = basis_.npwk[0];
    const int nh = cell_.atoms[0].ncpp.nh;
    Vector b(2 * nh, 0.0);
    for (int s = 0; s < 2; ++s)
    {
        for (int p = 0; p < nh; ++p)
        {
            for (int g = 0; g < npw; ++g)
            {
                b[s * nh + p] += std::conj(beta[(atom * nh + p) * npw + g]) * x[s * npw + g];
            }
        }
    }
    Vector y(x.size(), 0.0);
    for (int s = 0; s < 2; ++s)
    {
        for (int t = 0; t < 2; ++t)
        {
            for (int p = 0; p < nh; ++p)
            {
                for (int q = 0; q < nh; ++q)
                {
                    const Complex db = 0.5 * pp_.deeq_nc(2 * s + t, atom, p, q) * b[t * nh + q];
                    for (int g = 0; g < npw; ++g)
                    {
                        y[s * npw + g] += beta[(atom * nh + p) * npw + g] * db;
                    }
                }
            }
        }
    }
    return y;
}

Vector Observables::hamiltonian(const ModuleBase::Vector3<double>& a, const Vector& x)
{
    const Vector beta = projectors(a);
    Vector y(x.size(), 0.0);
    for (int atom = 0; atom < cell_.nat; ++atom)
    {
        const Vector term = nonlocal(beta, x, atom);
        for (size_t i = 0; i < x.size(); ++i)
        {
            y[i] += term[i];
        }
    }
    const int npw = basis_.npwk[0];
    for (int g = 0; g < npw; ++g)
    {
        const ModuleBase::Vector3<double> q = basis_.getgpluskcar(0, g) * cell_.tpiba + a;
        for (int s = 0; s < 2; ++s)
        {
            y[s * npw + g] += 0.5 * q.norm2() * x[s * npw + g];
        }
    }
    return y;
}

Vector pauli(const Vector& x, int axis)
{
    const Complex matrices[3][4] = {{0.0, 1.0, 1.0, 0.0}, {0.0, Complex(0, -1), Complex(0, 1), 0.0}, {1.0, 0.0, 0.0, -1.0}};
    const int npw = x.size() / 2;
    Vector y(x.size(), 0.0);
    for (int s = 0; s < 2; ++s)
    {
        for (int t = 0; t < 2; ++t)
        {
            for (int g = 0; g < npw; ++g)
            {
                y[s * npw + g] += matrices[axis][2 * s + t] * x[t * npw + g];
            }
        }
    }
    return y;
}

Complex inner(const Vector& x, const Vector& y)
{
    Complex value(0);
    for (size_t i = 0; i < x.size(); ++i)
    {
        value += std::conj(x[i]) * y[i];
    }
    return value;
}
} // namespace td_test

int main(int argc, char** argv)
{
    int nproc = 1;
    int nthread = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread, rank);
#ifdef __MPI
    POOL_WORLD = MPI_COMM_SELF;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
