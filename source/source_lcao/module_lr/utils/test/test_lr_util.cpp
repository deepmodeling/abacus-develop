#include "../lr_util.h"
#include "source_cell/unitcell.h"
#include <gtest/gtest.h>
#include <vector>

TEST(LRForcePP, AllowsForcesAndRelaxationWithoutNLCC)
{
    UnitCell cell;
    cell.ntype = 1;
    cell.atoms = new Atom[1];
    cell.set_atom_flag = true;
    EXPECT_NO_THROW(LR_Util::check_force_pp(cell, true, "scf"));
    EXPECT_NO_THROW(LR_Util::check_force_pp(cell, true, "relax"));
}

TEST(LRForcePP, AllowsNLCCTDDFTSpectra)
{
    UnitCell cell;
    cell.ntype = 2;
    cell.atoms = new Atom[2];
    cell.set_atom_flag = true;
    cell.pseudo_fn = {"H.upf", "Si_NLCC.upf"};
    cell.atoms[1].ncpp.nlcc = true;
    EXPECT_NO_THROW(LR_Util::check_force_pp(cell, false, "scf"));
    EXPECT_NO_THROW(LR_Util::check_force_pp(cell, false, "nscf"));
}

TEST(LRForcePP, RejectsNLCCForcesAndIdentifiesThePseudopotential)
{
    UnitCell cell;
    cell.ntype = 2;
    cell.atoms = new Atom[2];
    cell.set_atom_flag = true;
    cell.pseudo_fn = {"H.upf", "Si_NLCC.upf"};
    cell.atoms[1].ncpp.nlcc = true;
    cell.atoms[1].label = "Si";
    try
    {
        LR_Util::check_force_pp(cell, true, "scf");
        FAIL() << "NLCC forces were accepted";
    }
    catch (const std::invalid_argument& error)
    {
        const std::string message = error.what();
        EXPECT_NE(message.find("element Si"), std::string::npos);
        EXPECT_NE(message.find("Si_NLCC.upf"), std::string::npos);
        EXPECT_NE(message.find("not implemented"), std::string::npos);
    }
}

TEST(LRForcePP, RejectsNLCCRelaxationEvenWithoutExplicitForceFlag)
{
    UnitCell cell;
    cell.ntype = 2;
    cell.atoms = new Atom[2];
    cell.set_atom_flag = true;
    cell.pseudo_fn = {"H.upf", "Si_NLCC.upf"};
    cell.atoms[1].ncpp.nlcc = true;
    EXPECT_THROW(LR_Util::check_force_pp(cell, false, "relax"), std::invalid_argument);
    EXPECT_THROW(LR_Util::check_force_pp(cell, true, "relax"), std::invalid_argument);
}

TEST(LR_XCSigma, PreservesAntiparallelSpinGradientsWithoutHSE)
{
    // grad(up)=(2,0,0), grad(down)=(-1,0,0): a valid Gram matrix with sigma_ud < -1.
    const std::vector<double> sigma{4.0, -2.0, 1.0};
    std::vector<double> hse_buffer;
    const bool is_hse06 = false;
    const auto& input = LR_Util::prepare_xc_sigma(sigma, is_hse06, hse_buffer);
    EXPECT_EQ(input.data(), sigma.data());
    EXPECT_EQ(input[1], -2.0);
    EXPECT_TRUE(hse_buffer.empty());
}

TEST(LR_XCSigma, RetainsExistingHSEStabilizationWithoutMutatingTheDensity)
{
    const std::vector<double> sigma{0.01, -0.008, 0.0064};
    std::vector<double> hse_buffer;
    const bool is_hse06 = true;
    const auto& input = LR_Util::prepare_xc_sigma(sigma, is_hse06, hse_buffer);
    EXPECT_NE(input.data(), sigma.data());
    EXPECT_EQ(sigma[1], -0.008);
    EXPECT_EQ(input[0], sigma[0]);
    EXPECT_EQ(input[1], 1e-6);
    EXPECT_EQ(input[2], sigma[2]);
}
