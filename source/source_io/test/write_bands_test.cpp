#include "source_io/module_energy/write_bands.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_variable.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/module_parallel/para_world.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

/** Save and restore only the global settings consumed by band output. */
class TestParameters
{
  public:
    explicit TestParameters(Parameter& parameter) : parameter_(parameter)
    {
        append_ = parameter_.input.out_app_flag;
        directory_ = parameter_.sys.global_out_dir;
    }

    ~TestParameters()
    {
        parameter_.input.out_app_flag = append_;
        parameter_.sys.global_out_dir = directory_;
    }

    void set_output(const std::string& directory, bool append)
    {
        parameter_.sys.global_out_dir = directory;
        parameter_.input.out_app_flag = append;
    }

  private:
    Parameter& parameter_;
    bool append_ = false;
    std::string directory_;
};

namespace
{

using BandTopology = std::pair<int, int>;

// GoogleTest evaluates this generator during InitGoogleTest, after MPI_Init.
std::vector<BandTopology> band_topologies()
{
#ifdef __MPI
    const Parallel::ParaWorld world = Parallel::ParaWorld::make_mpi("test", MPI_COMM_WORLD);
    if (world.size() == 4)
    {
        return {{2, 1}, {1, 2}, {2, 2}};
    }
#endif
    return {{1, 1}};
}

std::string topology_name(const testing::TestParamInfo<BandTopology>& info)
{
    return "Kpar" + std::to_string(info.param.first) + "Bndpar" + std::to_string(info.param.second);
}

class BandOutputTest : public testing::TestWithParam<BandTopology>
{
  protected:
    std::unique_ptr<Parallel::ParaWorld> world_;
    std::unique_ptr<TestParameters> parameters_;
    K_Vectors kv_;
    ModuleBase::matrix ekb_;
    Input_para input_;
    std::string directory_;
    int kpar_ = 1;
    int bndpar_ = 1;
    int my_pool_ = 0;
    int band_group_ = 0;
    int rank_in_pool_ = 0;
    bool owns_pools_ = false;

    void SetUp() override
    {
#ifdef __MPI
        world_.reset(new Parallel::ParaWorld(Parallel::ParaWorld::make_mpi("test", MPI_COMM_WORLD)));
#else
        world_.reset(new Parallel::ParaWorld(Parallel::ParaWorld::serial("test")));
#endif
        kpar_ = GetParam().first;
        bndpar_ = GetParam().second;
        ASSERT_GT(kpar_, 0);
        ASSERT_GT(bndpar_, 0);
        ASSERT_EQ(world_->size() % (kpar_ * bndpar_), 0);
#ifdef __MPI
        int nproc_in_band_group = 0;
        int rank_in_band_group = 0;
        int nproc_in_pool = 0;
        Parallel_Global::init_pools(world_->size(), world_->rank(), bndpar_, kpar_,
                                    nproc_in_band_group, rank_in_band_group, band_group_,
                                    nproc_in_pool, rank_in_pool_, my_pool_);
        owns_pools_ = true;
#endif
        directory_ = "write_bands_p" + std::to_string(world_->size()) + "_"
                     + testing::UnitTest::GetInstance()->current_test_info()->name();
        std::replace(directory_.begin(), directory_.end(), '/', '_');
        directory_ += "/";
        bool ready = true;
        if (world_->rank() == 0)
        {
#ifdef _WIN32
            ready = _mkdir(directory_.c_str()) == 0 || errno == EEXIST;
#else
            ready = mkdir(directory_.c_str(), 0700) == 0 || errno == EEXIST;
#endif
            if (ready)
            {
                remove_outputs();
            }
        }
        Parallel::bcast_bool(ready, *world_, 0);
        ASSERT_TRUE(ready);
        parameters_.reset(new TestParameters(PARAM));
        parameters_->set_output(directory_, false);
    }

    void TearDown() override
    {
        parameters_.reset();
        if (world_)
        {
            Parallel::barrier(*world_);
            if (world_->rank() == 0 && !directory_.empty())
            {
                remove_outputs();
#ifdef _WIN32
                EXPECT_EQ(_rmdir(directory_.c_str()), 0);
#else
                EXPECT_EQ(rmdir(directory_.c_str()), 0);
#endif
            }
            Parallel::barrier(*world_);
        }
#ifdef __MPI
        if (owns_pools_)
        {
            MPI_Comm_free(&POOL_WORLD);
            if (KP_WORLD != MPI_COMM_NULL)
            {
                MPI_Comm_free(&KP_WORLD);
            }
            MPI_Comm_free(&INT_BGROUP);
            MPI_Comm_free(&BP_WORLD);
        }
#endif
    }

    void remove_outputs() const
    {
        for (const std::string& name : {"direct.txt", "up.txt", "down.txt", "band.txt", "bands1.txt", "bands2.txt"})
        {
            std::remove((directory_ + name).c_str());
        }
    }

    static double eigenvalue(int spin, int global_k, int global_band)
    {
        return 100.0 * spin + 10.0 * global_k + global_band - 2.0;
    }

    void prepare(int nspin, int nbands, bool distributed)
    {
        const int nk_per_spin = 5;
        const int spin_channels = nspin == 2 ? 2 : 1;
        int global_k_count = nk_per_spin;
        kv_.para_k.kinfo(global_k_count, kpar_, my_pool_, rank_in_pool_, world_->size(), spin_channels);
        const int local_k_count = kv_.para_k.nks_np;
        const int k_offset = my_pool_ * (nk_per_spin / kpar_) + std::min(my_pool_, nk_per_spin % kpar_);
        const int stored_rows = local_k_count * spin_channels;
        kv_.set_nks(stored_rows);
        kv_.set_nkstot(nk_per_spin * spin_channels);
        kv_.set_nkstot_nospin(nk_per_spin);
        kv_.kvec_c.resize(stored_rows);
        kv_.isk.resize(stored_rows);
        kv_.ik2iktot.resize(stored_rows);
        // Segment ids are global, as in the production k-point distribution.
        kv_.kl_segids = {0, 0, 1, 1, 1};
        const std::vector<ModuleBase::Vector3<double>> points{
            {0.0, 0.0, 0.0}, {0.3, 0.4, 0.0}, {1.3, 0.4, 0.0},
            {1.3, 0.4, 1.0}, {1.3, 0.4, 2.0}};
        const int local_bands = distributed ? nbands / bndpar_ + (band_group_ < nbands % bndpar_ ? 1 : 0) : nbands;
        const int band_offset = distributed ? band_group_ * (nbands / bndpar_) + std::min(band_group_, nbands % bndpar_) : 0;
        ekb_.create(stored_rows, local_bands);
        for (int spin = 0; spin < spin_channels; ++spin)
        {
            for (int ik = 0; ik < local_k_count; ++ik)
            {
                const int row = spin * local_k_count + ik;
                const int global_k = k_offset + ik;
                kv_.kvec_c[row] = points[global_k];
                kv_.isk[row] = spin;
                kv_.ik2iktot[row] = spin * nk_per_spin + global_k;
                for (int band = 0; band < local_bands; ++band)
                {
                    ekb_(row, band) = eigenvalue(spin, global_k, band_offset + band);
                }
            }
        }
        input_.nspin = nspin;
        input_.nbands = nbands;
        input_.out_band = {1, 8};
    }

    void expect_output(const std::string& name, int spin, int nbands,
                       double fermi, int precision, int copies) const
    {
        if (world_->rank() != 0)
        {
            return;
        }
        SCOPED_TRACE(name);
        std::ifstream file(directory_ + name);
        ASSERT_TRUE(file.is_open());
        const std::vector<double> distances{0.0, 0.5, 0.5, 1.5, 2.5};
        const double tolerance = 0.51 * std::pow(10.0, -precision);
        std::string line;
        for (int copy = 0; copy < copies; ++copy)
        {
            for (int ik = 0; ik < 5; ++ik)
            {
                ASSERT_TRUE(static_cast<bool>(std::getline(file, line))) << "Missing k-point " << ik;
                std::istringstream row(line);
                int index = 0;
                ASSERT_TRUE(static_cast<bool>(row >> index));
                EXPECT_EQ(index, ik + 1);
                std::vector<double> expected{distances[ik]};
                for (int band = 0; band < nbands; ++band)
                {
                    expected.push_back((eigenvalue(spin, ik, band) - fermi) * 13.605698);
                }
                for (const double value : expected)
                {
                    std::string token;
                    ASSERT_TRUE(static_cast<bool>(row >> token)) << line;
                    std::istringstream number(token);
                    double actual = 0.0;
                    ASSERT_TRUE(static_cast<bool>(number >> actual));
                    EXPECT_NEAR(actual, value, tolerance);
                    const std::size_t dot = token.find('.');
                    ASSERT_NE(dot, std::string::npos);
                    EXPECT_EQ(token.size() - dot - 1, static_cast<std::size_t>(precision));
                    EXPECT_TRUE(number.eof());
                }
                std::string extra;
                EXPECT_FALSE(static_cast<bool>(row >> extra)) << "Unexpected column: " << extra;
            }
        }
        EXPECT_FALSE(static_cast<bool>(std::getline(file, line))) << "Unexpected line: " << line;
    }

    std::string read_text(const std::string& name) const
    {
        std::ifstream file(directory_ + name);
        return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
    }

    void seed_file(const std::string& name, const std::string& content) const
    {
        if (world_->rank() == 0)
        {
            std::ofstream file(directory_ + name);
            ASSERT_TRUE(file.is_open());
            file << content;
        }
    }

    void expect_missing(const std::string& name) const
    {
        if (world_->rank() == 0)
        {
            EXPECT_FALSE(std::ifstream(directory_ + name).is_open()) << name;
        }
    }

    void check_entrypoint(int nspin)
    {
        for (const bool distributed : {false, true})
        {
            SCOPED_TRACE(distributed ? "distributed" : "replicated");
            prepare(nspin, 3, distributed);
            ModuleIO::write_bands(input_, ekb_, kv_);
            if (nspin == 2)
            {
                expect_output("bands1.txt", 0, 3, 0.0, 8, 1);
                expect_output("bands2.txt", 1, 3, 0.0, 8, 1);
                expect_missing("band.txt");
            }
            else
            {
                expect_output("band.txt", 0, 3, 0.0, 8, 1);
                expect_missing("bands1.txt");
                expect_missing("bands2.txt");
            }
            // Root must finish reading before another collective call overwrites output.
            Parallel::barrier(*world_);
        }
    }
};

TEST_P(BandOutputTest, DirectSingleSpin)
{
    prepare(1, 3, false);
    ModuleIO::nscf_bands(0, directory_ + "direct.txt", 3, 0.0, 8, ekb_, kv_);
    expect_output("direct.txt", 0, 3, 0.0, 8, 1);
    if (world_->rank() == 0)
    {
        const std::string text = read_text("direct.txt");
        EXPECT_EQ(text.substr(0, text.find('\n')),
                  "   1 0.00000000 -27.21139600 -13.60569800 0.00000000");
    }
}

TEST_P(BandOutputTest, DirectBothSpins)
{
    prepare(2, 3, false);
    ModuleIO::nscf_bands(0, directory_ + "up.txt", 3, 0.0, 8, ekb_, kv_);
    ModuleIO::nscf_bands(1, directory_ + "down.txt", 3, 0.0, 8, ekb_, kv_);
    expect_output("up.txt", 0, 3, 0.0, 8, 1);
    expect_output("down.txt", 1, 3, 0.0, 8, 1);
}

TEST_P(BandOutputTest, PrecisionAndFermiShift)
{
    prepare(1, 3, false);
    for (const int precision : {4, 8})
    {
        ModuleIO::nscf_bands(0, directory_ + "direct.txt", 3, 0.25, precision, ekb_, kv_);
        expect_output("direct.txt", 0, 3, 0.25, precision, 1);
        if (world_->rank() == 0 && precision == 4)
        {
            const std::string text = read_text("direct.txt");
            EXPECT_EQ(text.substr(0, text.find('\n')), "   1 0.0000 -30.6128 -17.0071 -3.4014");
        }
        Parallel::barrier(*world_);
    }
}

TEST_P(BandOutputTest, OverwriteAndAppend)
{
    prepare(1, 3, true);
    seed_file("band.txt", "stale output\n");
    Parallel::barrier(*world_);
    ModuleIO::write_bands(input_, ekb_, kv_);
    expect_output("band.txt", 0, 3, 0.0, 8, 1);
    Parallel::barrier(*world_);
    parameters_->set_output(directory_, true);
    ModuleIO::write_bands(input_, ekb_, kv_);
    expect_output("band.txt", 0, 3, 0.0, 8, 2);
}

TEST_P(BandOutputTest, OutputDisabled)
{
    prepare(1, 3, true);
    input_.out_band[0] = 0;
    ModuleIO::write_bands(input_, ekb_, kv_);
    expect_missing("band.txt");
    seed_file("band.txt", "keep this content\n");
    Parallel::barrier(*world_);
    ModuleIO::write_bands(input_, ekb_, kv_);
    if (world_->rank() == 0)
    {
        EXPECT_EQ(read_text("band.txt"), "keep this content\n");
    }
    expect_missing("bands1.txt");
    expect_missing("bands2.txt");
}

TEST_P(BandOutputTest, SingleSpinFiles)
{
    check_entrypoint(1);
}

TEST_P(BandOutputTest, CollinearSpinFiles)
{
    check_entrypoint(2);
}

TEST_P(BandOutputTest, SpinorFile)
{
    check_entrypoint(4);
}

TEST_P(BandOutputTest, EmptyBandShard)
{
    prepare(2, 1, true);
    if (bndpar_ > 1 && band_group_ == bndpar_ - 1)
    {
        EXPECT_EQ(ekb_.nc, 0);
    }
    ModuleIO::write_bands(input_, ekb_, kv_);
    expect_output("bands1.txt", 0, 1, 0.0, 8, 1);
    expect_output("bands2.txt", 1, 1, 0.0, 8, 1);
}

INSTANTIATE_TEST_SUITE_P(Topology, BandOutputTest, testing::ValuesIn(band_topologies()), topology_name);

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
    int nproc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (nproc != 1 && nproc != 4)
    {
        std::fprintf(stderr, "Band output tests require 1 or 4 MPI processes.\n");
        MPI_Finalize();
        return 1;
    }
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
