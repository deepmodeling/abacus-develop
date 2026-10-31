#include "dfpt_serial_fixture.h"
#include "source_pw/module_dfpt/dfpt_write_dyn.h"

#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/************************************************
 *  unit test of dfpt_write_dyn.cpp and the
 *  Gamma-first q-ordering validation
 *  (validate_q_ordering in dfpt_pw_data.cpp)
 ***********************************************/

/// Writer-focused fixture: only the two-atom cell geometry and the
/// qlist/data binding are needed; the plane-wave bases of the base
/// SetUp are not on the writer path.
class DFPTWriteDynTest : public DFPTSerialBase
{
  protected:
    void SetUp() override
    {
        SetUpCell();
        MakeTwoAtomCell();
        // bind data to qlist exactly like SetupBases does (nat = 2)
        data_.init(&qlist_, 1, 1, 1, 1, 1, 2, nullptr);
    }

    /// install an explicit q list (Gamma-first or not) into qlist_
    void set_qvecs(const std::vector<ModuleBase::Vector3<double>>& qvecs)
    {
        qlist_.kvec_d = qvecs;
        qlist_.nkstot = static_cast<int>(qvecs.size());
    }

    /// a distinguishable 6x6 raw matrix: dm(row, col) = (row + 1) + col
    ModuleBase::ComplexMatrix make_raw() const
    {
        ModuleBase::ComplexMatrix dm(6, 6);
        for (int r = 0; r < 6; ++r)
        {
            for (int c = 0; c < 6; ++c)
            {
                dm(r, c) = std::complex<double>(r + 1, c);
            }
        }
        return dm;
    }

    void set_polar_tensors()
    {
        ModuleBase::matrix eps(3, 3);
        for (int i = 0; i < 3; ++i)
        {
            eps(i, i) = 2.5;
        }
        data_.set_dielectric(eps);
        for (int ia = 0; ia < 2; ++ia)
        {
            ModuleBase::matrix zstar(3, 3);
            zstar(0, 0) = 1.0 + ia;
            zstar(1, 1) = 1.0 + ia;
            zstar(2, 2) = 1.0 + ia;
            data_.set_born(ia, zstar);
        }
    }

    std::string read_file(const std::string& path) const
    {
        std::ifstream in(path);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }
};

TEST_F(DFPTWriteDynTest, GammaOnlyPolarRecordsAndInjectedLog)
{
    // Gamma-first list [Gamma, X]; compute_q0 on, LO-TO not consulted
    set_qvecs({ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.5, 0.0, 0.0)});
    data_.set_dynmat_raw(0, make_raw());
    data_.set_dynmat_raw(1, make_raw());
    set_polar_tensors();
    data_.set_compute_q0(true);

    const std::string dir = testing::TempDir();
    std::ostringstream running_log;
    DFPT_Write_Dyn(ucell_, data_, 2, dir, 0, &running_log);

    const std::string dyn0 = read_file(dir + "/dyn.q0");
    const std::string dyn1 = read_file(dir + "/dyn.q1");

    // polar records only on the Gamma entry, Dielectric before Born (QE order)
    EXPECT_NE(dyn0.find("Dielectric Tensor"), std::string::npos);
    EXPECT_NE(dyn0.find("Effective Charges E-U"), std::string::npos);
    EXPECT_LT(dyn0.find("Dielectric Tensor"), dyn0.find("Effective Charges E-U"));
    EXPECT_EQ(dyn1.find("Dielectric Tensor"), std::string::npos);
    EXPECT_EQ(dyn1.find("Effective Charges E-U"), std::string::npos);

    // the raw snapshot is exported (dm(0,0) = 1+0i, masses 12 and 12 for
    // the (1,1) block: Phi = sqrt(12*12) * 1 = 12)
    EXPECT_NE(dyn0.find("1.2000000000e+01"), std::string::npos);

    // both created-file confirmations land in the injected log stream
    EXPECT_NE(running_log.str().find("DFPT_Write_Dyn: wrote " + dir + "/dyn.q0"), std::string::npos);
    EXPECT_NE(running_log.str().find("DFPT_Write_Dyn: wrote " + dir + "/dyn.q1"), std::string::npos);

    EXPECT_EQ(std::remove((dir + "/dyn.q0").c_str()), 0);
    EXPECT_EQ(std::remove((dir + "/dyn.q1").c_str()), 0);
}

TEST_F(DFPTWriteDynTest, NonzeroRankWritesNothing)
{
    set_qvecs({ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    data_.set_dynmat_raw(0, make_raw());
    data_.set_compute_q0(true);

    const std::string dir = testing::TempDir();
    std::ostringstream running_log;
    DFPT_Write_Dyn(ucell_, data_, 1, dir, 1, &running_log);
    EXPECT_TRUE(running_log.str().empty());

    std::ifstream in(dir + "/dyn.q0");
    EXPECT_FALSE(in.is_open());
}

TEST_F(DFPTWriteDynTest, GammaFirstOrderingAccepted)
{
    set_qvecs({ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.5, 0.5),
               ModuleBase::Vector3<double>(0.25, 0.25, 0.25)});
    EXPECT_NO_THROW(ModuleDFPT::validate_q_ordering(qlist_, true));
}

TEST_F(DFPTWriteDynTest, NonFirstGammaRejectedWhenComputeQ0)
{
    // [X, Gamma]: the q0 pipeline would run at the nonzero-q index 0
    set_qvecs({ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleDFPT::validate_q_ordering(qlist_, true), ::testing::ExitedWithCode(1), "");
    const std::string first_diag = testing::internal::GetCapturedStdout();
    EXPECT_NE(first_diag.find("FIRST q point"), std::string::npos);

    // no Gamma at all: same rejection branch
    set_qvecs({ModuleBase::Vector3<double>(0.5, 0.0, 0.0)});
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleDFPT::validate_q_ordering(qlist_, true), ::testing::ExitedWithCode(1), "");
    const std::string no_gamma_diag = testing::internal::GetCapturedStdout();
    EXPECT_NE(no_gamma_diag.find("FIRST q point"), std::string::npos);
}

TEST_F(DFPTWriteDynTest, DuplicateGammaRejectedWhenComputeQ0)
{
    // [Gamma, X, Gamma]: the duplicate at index 2 would silently miss the
    // polar tensors while the first entry carries them
    set_qvecs({ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleDFPT::validate_q_ordering(qlist_, true), ::testing::ExitedWithCode(1), "");
    const std::string dup_diag = testing::internal::GetCapturedStdout();
    EXPECT_NE(dup_diag.find("duplicate Gamma q point at index 2"), std::string::npos);
}

TEST_F(DFPTWriteDynTest, ArbitraryOrderingAllowedWhenComputeQ0Off)
{
    // plain phonon runs keep full freedom: [X, Gamma] and duplicates are
    // legal when the polar pipeline is disabled
    set_qvecs({ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    EXPECT_NO_THROW(ModuleDFPT::validate_q_ordering(qlist_, false));

    set_qvecs({ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
               ModuleBase::Vector3<double>(0.0, 0.0, 0.0)});
    EXPECT_NO_THROW(ModuleDFPT::validate_q_ordering(qlist_, false));
}
