#include "source_estate/module_pot/sccs_output.h"
#include "source_estate/module_pot/pot_sccs.h"
#include "source_estate/module_pot/pot_pcc.h"
#include "source_hamilt/module_sccs/test/sccs_test.h"
#include "source_estate/module_pot/sol_grid_field.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_base/parallel_reduce.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const elecstate::SccsResume no_resume = elecstate::SccsResume();
}

using SccsOutputTest = SccsTest::PwTest;

TEST_F(SccsOutputTest, OutputLevelsSelectTheReportsAndPreserveStreamState)
{
    elecstate::SccsOutput result;
    result.valid = true;
    result.iterations = 4;
    result.reaction_energy = -0.1;
    std::ostringstream output;
    output.setf(std::ios::scientific, std::ios::floatfield);
    output.precision(3);
    const std::ios::fmtflags flags = output.flags();
    elecstate::write_sccs_out(output, result, 0, -0.2, true, 0.05);
    EXPECT_TRUE(output.str().empty());
    elecstate::write_sccs_out(output, result, 1, -0.2, true, 0.05);
    // Cells of 16 characters, each starting on a column edge.
    const std::string iteration_line = " SCCS_ITER       CG_STEPS        4               E_SOL/Ry        "
                                       "-2.00000000e-01 E_PCC/Ry        5.00000000e-02\n";
    EXPECT_EQ(output.str(), iteration_line);
    EXPECT_EQ(output.str().find("SCCS_RESIDUAL"), std::string::npos);
    output.str("");
    elecstate::write_sccs_out(output, result, 2, -0.2, true, 0.05);
    // Every cell starts on the grid of 16 columns: no cell runs into the next,
    // every cell begins on its column edge, and so do the table headers.
    const std::size_t cell = 16;
    std::istringstream lines(output.str());
    std::string line;
    int table_rows = 0;
    while (std::getline(lines, line))
    {
        for (std::size_t start = 1 + cell; start < line.size(); start += cell)
        {
            EXPECT_EQ(line[start - 1], ' ') << line;
            EXPECT_NE(line[start], ' ') << line;
        }
        if (line.find("PCC0D_MOMENTS") == std::string::npos)
        {
            continue;
        }
        ++table_rows;
        if (line.find("KIND") != std::string::npos)
        {
            EXPECT_EQ(line.find("Q/e"), 1 + 2 * cell);
            EXPECT_EQ(line.find("QRR/eBohr2"), 1 + 6 * cell);
        }
    }
    EXPECT_EQ(table_rows, 5);
    EXPECT_NE(output.str().find(" PCC0D_MOMENTS   SCREENED"), std::string::npos);
    EXPECT_NE(output.str().find("CACHED_SOURCES  0"), std::string::npos);
    EXPECT_EQ(output.str().find("V_FILLED"), std::string::npos);
    EXPECT_EQ(output.precision(), 3);
    EXPECT_EQ(output.flags(), flags);
    output.str("");
    result.pcc.slab = true;
    result.pcc.normal.z = 1.0;
    elecstate::write_sccs_out(output, result, 2, -0.2, true, 0.05);
    EXPECT_NE(output.str().find("PCC2D_ORIGIN"), std::string::npos);
    EXPECT_NE(output.str().find(" PCC2D_MOMENTS   SCREENED"), std::string::npos);
    EXPECT_EQ(output.precision(), 3);
    EXPECT_EQ(output.flags(), flags);
}

TEST_F(SccsOutputTest, CellAsWideAsTheColumnKeepsOneBlank)
{
    elecstate::SccsOutput result;
    result.valid = true;
    result.iterations = 4;
    std::ostringstream output;
    // A three-digit exponent fills the whole 16-character cell.
    elecstate::write_sccs_out(output, result, 1, -1.5e-100, true, 0.05);
    EXPECT_NE(output.str().find("-1.50000000e-100 E_PCC/Ry"), std::string::npos) << output.str();
}

TEST_F(SccsOutputTest, AdvancedDiagnosticsPreservePccPotentialEnergyAndForce)
{
    UnitCell cell;
    Atom atom;
    cell.lat0 = length;
    cell.tpiba = tpiba;
    cell.omega = basis.omega;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.mass = 1.0;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    std::vector<double> density(basis.nrxx, 1e-6);
    const double electron_count = 1e-6 * basis.omega;
    double* channel = density.data();
    Charge charge;
    charge.nspin = 1;
    charge.rho = &channel;
    const char* boundaries[] = {"pcc_0d", "pcc_2d"};
    for (const char* boundary : boundaries)
    {
        Input_para input;
        input.assume_isolated = boundary;
        input.sccs_epsilon = 5.0;
        ModuleSccs::SccsConfig config;
        ModuleSccs::PolarizationSolverParameters solver;
        elecstate::sccs_cfg_input(input, config, solver);
        elecstate::PotSccs plain(&basis, config, solver, electron_count, no_resume);
        config.collect_report = true;
        elecstate::PotSccs checked(&basis, config, solver, electron_count, no_resume);
        elecstate::PotPcc::Dimension dimension = elecstate::PotPcc::Dimension::molecule;
        if (config.boundary == ModuleSccs::Boundary::Pcc2d)
        {
            dimension = elecstate::PotPcc::Dimension::slab;
        }
        elecstate::PotPcc pcc(&basis, dimension, 2, electron_count);
        ModuleBase::matrix vacuum(1, basis.nrxx);
        pcc.cal_v_eff(&charge, &cell, vacuum);
        ModuleBase::matrix first(1, basis.nrxx);
        ModuleBase::matrix second(1, basis.nrxx);
        plain.cal_v_eff(&charge, &cell, first);
        checked.cal_v_eff(&charge, &cell, second);
        EXPECT_DOUBLE_EQ(plain.get_energy(), checked.get_energy());
        for (int i = 0; i < basis.nrxx; ++i)
        {
            EXPECT_DOUBLE_EQ(first(0, i), second(0, i));
        }
        ModuleBase::matrix force1(1, 3);
        ModuleBase::matrix force2(1, 3);
        plain.add_sol_force(cell, force1);
        checked.add_sol_force(cell, force2);
        for (int i = 0; i < 3; ++i)
        {
            EXPECT_DOUBLE_EQ(force1(0, i), force2(0, i));
        }
        const double pcc_energy = pcc.get_energy();
        std::ostringstream summary;
        checked.write_iter(summary, 2, 1e-6, pcc_energy);
        EXPECT_NE(summary.str().find("SCCS_GAUSS"), std::string::npos);
        EXPECT_TRUE(checked.has_output());
        config.start_drho = 1e-4;
        elecstate::PotSccs deferred(&basis, config, solver, electron_count, no_resume);
        EXPECT_FALSE(deferred.is_active());
        EXPECT_FALSE(deferred.has_output());
        std::ostringstream waiting;
        deferred.write_iter(waiting, 1, 1e-3, pcc_energy);
        EXPECT_NE(waiting.str().find("SCCS_DEFERRED"), std::string::npos);
    }
}

// out_sol fields of PotSccs, written by ctrl_output_fp with collective calls.
using SccsSolventFieldTest = SccsTest::PwTest;

TEST_F(SccsSolventFieldTest, SolventFieldsFollowTheLastActiveEvaluation)
{
    UnitCell cell;
    Atom atom;
    cell.lat0 = length;
    cell.tpiba = tpiba;
    cell.omega = basis.omega;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = &atom;
    atom.na = 1;
    atom.ncpp.zv = 1.0;
    atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.25)};
    // A density across the switching range gives a nontrivial cavity.
    std::vector<double> density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double step = ir % 7;
        density[ir] = 1e-5 * (1.0 + step);
    }
    double* channels[] = {density.data()};
    Charge charge;
    charge.nspin = 1;
    charge.rho = channels;
    ModuleSccs::SccsConfig config;
    ModuleSccs::PolarizationSolverParameters solver;
    Input_para input;
    elecstate::sccs_cfg_input(input, config, solver);
    config.cavity.epsilon_bulk = 5.0;
    config.cavity.density_min = 2e-5;
    config.cavity.density_max = 5e-5;
    config.start_drho = 1e-3;
    elecstate::PotSccs delayed(&basis, config, solver, 1.0, no_resume);
    std::vector<elecstate::SolventGridField> fields;
    delayed.add_sol_fields(fields);
    EXPECT_TRUE(fields.empty());
    config.start_drho = 0.0;
    elecstate::PotSccs active(&basis, config, solver, 1.0, no_resume);
    active.add_sol_fields(fields);
    EXPECT_TRUE(fields.empty());
    ModuleBase::matrix potential(1, basis.nrxx);
    active.cal_v_eff(&charge, &cell, potential);
    active.add_sol_fields(fields);
    ASSERT_EQ(fields.size(), 2u);
    EXPECT_EQ(fields[0].name, "eps");
    EXPECT_EQ(fields[1].name, "cavity");
    ASSERT_EQ(fields[0].values.size(), static_cast<std::size_t>(basis.nrxx));
    const double log_bulk = std::log(config.cavity.epsilon_bulk);
    bool transition = false;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double solute = fields[1].values[ir];
        EXPECT_GE(solute, 0.0);
        EXPECT_LE(solute, 1.0);
        const double expected = std::exp(log_bulk * (1.0 - solute));
        EXPECT_NEAR(fields[0].values[ir], expected, 1e-12);
        transition = transition || (solute > 0.0 && solute < 1.0);
    }
    EXPECT_TRUE(transition);
}

// With more pool ranks than z planes a rank owns no real-space grid points.
class PotSccsEmptyPartitionTest : public SccsTest::EmptyPartitionTest
{
protected:
    // One active evaluation, then the out_sol field list. The pool sum over the
    // field count is collective, so a rank with a different list hangs.
    void evaluate_and_collect(const ModuleSccs::SccsConfig& config,
                            double& energy,
                            std::vector<std::string>& names,
                            std::vector<double>& sums)
    {
        double empty_ranks = basis.nrxx == 0 ? 1.0 : 0.0;
        Parallel_Reduce::reduce_pool(empty_ranks);
        if (SccsTest::pool_size >= 3)
        {
            EXPECT_GE(empty_ranks, 1.0);
        }
        UnitCell cell;
        Atom atom;
        cell.lat0 = length;
        cell.tpiba = tpiba;
        cell.omega = basis.omega;
        cell.latvec = lattice;
        cell.ntype = 1;
        cell.nat = 1;
        cell.atoms = &atom;
        atom.label = "H";
        atom.na = 1;
        atom.ncpp.zv = 1.0;
        atom.tau = {ModuleBase::Vector3<double>(0.25, 0.25, 0.1)};
        const std::vector<double> mode = cosine_mode(0);
        const double mean_density = 3.5e-5;
        std::vector<double> density(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            density[ir] = mean_density * (1.0 + 0.8 * mode[ir]);
        }
        double* channels[] = {density.data()};
        Charge charge;
        charge.nspin = 1;
        charge.rho = channels;
        ModuleSccs::PolarizationSolverParameters solver;
        const double electron_count = mean_density * basis.omega;
        elecstate::PotSccs component(&basis, config, solver, electron_count, no_resume);
        ModuleBase::matrix potential(1, basis.nrxx);
        component.cal_v_eff(&charge, &cell, potential);
        energy = component.get_energy();

        std::vector<elecstate::SolventGridField> fields;
        component.add_sol_fields(fields);
        const std::size_t field_count = fields.size();
        names.clear();
        sums.assign(field_count, 0.0);
        for (std::size_t field = 0; field < field_count; ++field)
        {
            names.push_back(fields[field].name);
            ASSERT_EQ(fields[field].values.size(), static_cast<std::size_t>(basis.nrxx));
            for (double value : fields[field].values)
            {
                sums[field] += value;
            }
        }
        // Collective: a rank with a different field count would hang here.
        Parallel_Reduce::reduce_pool(sums.data(), static_cast<int>(field_count));
    }

    ModuleSccs::SccsConfig active_config() const
    {
        ModuleSccs::SccsConfig config;
        ModuleSccs::PolarizationSolverParameters solver;
        Input_para input;
        elecstate::sccs_cfg_input(input, config, solver);
        config.cavity.epsilon_bulk = 5.0;
        config.cavity.density_min = 2e-5;
        config.cavity.density_max = 5e-5;
        config.start_drho = 0.0;
        return config;
    }
};

// Every rank lists eps and cavity, including one without grid points, and the
// energy and fields match the serial run.
TEST_F(PotSccsEmptyPartitionTest, SolventOutputIsCollectiveOnEmptyGridPartitions)
{
    double energy = 0.0;
    std::vector<std::string> names;
    std::vector<double> sums;
    evaluate_and_collect(active_config(), energy, names, sums);
    const std::vector<std::string> expected_names = {"eps", "cavity"};
    ASSERT_EQ(names, expected_names);
    // Serial reference values.
    EXPECT_NEAR(energy, -0.748979082383938, 1e-10);
    EXPECT_NEAR(sums[0], 795.834560078087, 1e-8);
    EXPECT_NEAR(sums[1], 154.256465345260, 1e-8);
}
