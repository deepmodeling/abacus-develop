#include "source_cell/mdcell_cif_io.h"

#include "source_base/parallel_cell.h"
#include "source_cell/cif_io.h"
#include "source_cell/mdcell.h"

#include "gtest/gtest.h"

#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

TEST(MDCellCifIOTest, WritesAtomsInStableIdentityOrder)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 4.0;
    lattice.e22 = 4.0;
    lattice.e33 = 4.0;
    const ModuleBase::Matrix3 reciprocal = lattice.Inverse();
    std::vector<LocalAtom> atoms;
    for (int index = 2; index >= 0; --index)
    {
        if (index % domain.size() == domain.rank())
        {
            LocalAtom atom;
            atom.type = 0;
            atom.type_index = index;
            atom.frac.set(0.1 * (index + 1), 0.2, 0.3);
            atoms.push_back(atom);
        }
    }
    MDCell cell;
    const std::vector<std::string> labels(1, "He");
    const std::vector<double> masses(1, 4.0);
    const std::vector<std::int64_t> counts(1, 3);
    cell.initialize_from_owned_atoms(lattice, reciprocal, 1.0, 64.0, 3, atoms, labels, masses, counts, 0.0, domain);
    const std::string filename = testing::TempDir() + "mdcell_cif_order.cif";
    mdcell::print_cif_file(cell, filename, "# MD step 1");
    if (domain.rank() == 0)
    {
        std::map<std::string, std::vector<std::string> > data;
        ModuleIO::CifParser::read(filename, data);
        ASSERT_EQ(data["_atom_site_fract_x"].size(), 3U);
        EXPECT_NEAR(std::stod(data["_atom_site_fract_x"][0]), 0.1, 1e-8);
        EXPECT_NEAR(std::stod(data["_atom_site_fract_x"][1]), 0.2, 1e-8);
        EXPECT_NEAR(std::stod(data["_atom_site_fract_x"][2]), 0.3, 1e-8);
        EXPECT_NEAR(std::stod(data["_cell_length_a"][0]), 4.0 * 0.52917721067, 1e-7);
        std::remove(filename.c_str());
    }
}

TEST(MDCellCifIOTest, RejectsDuplicateIdentity)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 1.0;
    lattice.e22 = 1.0;
    lattice.e33 = 1.0;
    const ModuleBase::Matrix3 reciprocal = lattice.Inverse();
    std::vector<LocalAtom> atoms;
    if (domain.rank() == 0)
    {
        LocalAtom atom;
        atoms.push_back(atom);
        atoms.push_back(atom);
    }
    MDCell cell;
    const std::vector<std::string> labels(1, "He");
    const std::vector<double> masses(1, 4.0);
    const std::vector<std::int64_t> counts(1, 1);
    cell.initialize_from_owned_atoms(lattice, reciprocal, 1.0, 1.0, 1, atoms, labels, masses, counts, 0.0, domain);
    const std::string filename = testing::TempDir() + "mdcell_cif_duplicate.cif";
    EXPECT_THROW(mdcell::print_cif_file(cell, filename, "# MD step 1"), std::runtime_error);
}
