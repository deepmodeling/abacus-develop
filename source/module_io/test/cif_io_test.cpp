#include <gtest/gtest.h>
#include "module_io/cif_io.h"
#include <cmath>
#include <random>
#include "module_base/formatter.h"
#include <fstream>

/**
 * this is the unittest for ABACUS i/o interface with Crystal Information File (CIF) format.
 * 
 * Intergrately, there are two examples, one easier, one more complicated.
 */
const std::string mp2516584 = ""
"# generated using pymatgen\n"
"data_C\n"
"_symmetry_space_group_name_H-M   'P 1'\n"
"_cell_length_a   2.46772428\n"
"_cell_length_b   2.46772428\n"
"_cell_length_c   8.68503800\n"
"_cell_angle_alpha   90.00000000\n"
"_cell_angle_beta   90.00000000\n"
"_cell_angle_gamma   120.00000758\n"
"_symmetry_Int_Tables_number   1\n"
"_chemical_formula_structural   C\n"
"_chemical_formula_sum   C4\n"
"_cell_volume   45.80317575\n"
"_cell_formula_units_Z   4\n"
"loop_\n"
" _symmetry_equiv_pos_site_id\n"
" _symmetry_equiv_pos_as_xyz\n"
"  1  'x, y, z'\n"
"loop_\n"
" _atom_site_type_symbol\n"
" _atom_site_label\n"
" _atom_site_symmetry_multiplicity\n"
" _atom_site_fract_x\n"
" _atom_site_fract_y\n"
" _atom_site_fract_z\n"
" _atom_site_occupancy\n"
"  C  C0  1  0.00000000  0.00000000  0.75000000  1\n"
"  C  C1  1  0.00000000  0.00000000  0.25000000  1\n"
"  C  C2  1  0.33333300  0.66666700  0.75000000  1\n"
"  C  C3  1  0.66666700  0.33333300  0.25000000  1\n";

const std::string cod1000065 = ""
"#------------------------------------------------------------------------------\n"
"#$Date: 2016-02-18 17:37:37 +0200 (Thu, 18 Feb 2016) $\n"
"#$Revision: 176729 $\n"
"#$URL: svn://www.crystallography.net/cod/cif/1/00/00/1000065.cif $\n"
"#------------------------------------------------------------------------------\n"
"#\n"
"# This file is available in the Crystallography Open Database (COD),\n"
"# http://www.crystallography.net/\n"
"#\n"
"# All data on this site have been placed in the public domain by the\n"
"# contributors.\n"
"#\n"
"data_1000065\n"
"loop_\n"
"_publ_author_name\n"
"'Nixon, D E'\n"
"'Parry, G S'\n"
"'Ubbelohde, A R'\n"
"_publ_section_title\n"
";\n"
"Order-disorder transformations in graphite nitrates\n"
";\n"
"_journal_coden_ASTM              PRLAAZ\n"
"_journal_name_full\n"
";\n"
"Proceedings of the Royal Society of London, Series A: Mathematical and\n" 
"Physical Sciences (76,1906-)\n"
";\n"
"_journal_page_first              324\n"
"_journal_page_last               339\n"
"_journal_paper_doi               10.1098/rspa.1966.0098\n"
"_journal_volume                  291\n"
"_journal_year                    1966\n"
"_chemical_formula_analytical     'C (H N O3)'\n"
"_chemical_formula_structural     C\n"
"_chemical_formula_sum            C\n"
"_chemical_name_common            'Graphite nitrate'\n"
"_chemical_name_systematic        Carbon\n"
"_space_group_IT_number           166\n"
"_symmetry_cell_setting           trigonal\n"
"_symmetry_space_group_name_Hall  '-R 3 2\"'\n"
"_symmetry_space_group_name_H-M   'R -3 m :H'\n"
"_cell_angle_alpha                90\n"
"_cell_angle_beta                 90\n"
"_cell_angle_gamma                120\n"
"_cell_formula_units_Z            12\n"
"_cell_length_a                   2.46\n"
"_cell_length_b                   2.46\n"
"_cell_length_c                   33.45\n"
"_cell_volume                     175.3\n"
"_cod_original_sg_symbol_H-M      'R -3 m H'\n"
"_cod_database_code               1000065\n"
"loop_\n"
"_symmetry_equiv_pos_as_xyz\n"
"x,y,z\n"
"-y,x-y,z\n"
"y-x,-x,z\n"
"-y,-x,z\n"
"x,x-y,z\n"
"y-x,y,z\n"
"-x,-y,-z\n"
"y,y-x,-z\n"
"x-y,x,-z\n"
"y,x,-z\n"
"-x,y-x,-z\n"
"x-y,-y,-z\n"
"1/3+x,2/3+y,2/3+z\n"
"2/3+x,1/3+y,1/3+z\n"
"1/3-y,2/3+x-y,2/3+z\n"
"2/3-y,1/3+x-y,1/3+z\n"
"1/3-x+y,2/3-x,2/3+z\n"
"2/3-x+y,1/3-x,1/3+z\n"
"1/3-y,2/3-x,2/3+z\n"
"2/3-y,1/3-x,1/3+z\n"
"1/3+x,2/3+x-y,2/3+z\n"
"2/3+x,1/3+x-y,1/3+z\n"
"1/3-x+y,2/3+y,2/3+z\n"
"2/3-x+y,1/3+y,1/3+z\n"
"1/3-x,2/3-y,2/3-z\n"
"2/3-x,1/3-y,1/3-z\n"
"1/3+y,2/3-x+y,2/3-z\n"
"2/3+y,1/3-x+y,1/3-z\n"
"1/3+x-y,2/3+x,2/3-z\n"
"2/3+x-y,1/3+x,1/3-z\n"
"1/3+y,2/3+x,2/3-z\n"
"2/3+y,1/3+x,1/3-z\n"
"1/3-x,2/3-x+y,2/3-z\n"
"2/3-x,1/3-x+y,1/3-z\n"
"1/3+x-y,2/3-y,2/3-z\n"
"2/3+x-y,1/3-y,1/3-z\n"
"loop_\n"
"_atom_site_label\n"
"_atom_site_type_symbol\n"
"_atom_site_symmetry_multiplicity\n"
"_atom_site_Wyckoff_symbol\n"
"_atom_site_fract_x\n"
"_atom_site_fract_y\n"
"_atom_site_fract_z\n"
"_atom_site_occupancy\n"
"_atom_site_attached_hydrogens\n"
"_atom_site_calc_flag\n"
"C1 C0 6 c 0. 0. 0.05 1. 0 d\n"
"C2 C0 6 c 0. 0. 0.283 1. 0 d\n"
"loop_\n"
"_atom_type_symbol\n"
"_atom_type_oxidation_number\n"
"C0 0.000\n";

TEST(CifParserTest, ReadSimpleTest)
{
    int rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef __MPI
    if (rank == 0)
    {
#endif
    std::ofstream ofs("mp-2516584.cif");
    ofs << mp2516584;
    ofs.close();
#ifdef __MPI
    }
    MPI_Barrier(MPI_COMM_WORLD); // make sure the file is written
#endif
    std::map<std::string, std::vector<std::string>> data;
    ModuleIO::CifParser::read("mp-2516584.cif", data);
    // delete the file
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD); // make sure the file is already read
    if (rank == 0)
    {
#endif
    std::remove("mp-2516584.cif");
#ifdef __MPI
    }
#endif

    EXPECT_EQ(data.size(), 23);
    EXPECT_EQ(data["_symmetry_space_group_name_H-M"][0], "'P 1'");
    EXPECT_EQ(data["_cell_length_a"][0], "2.46772428");
    EXPECT_EQ(data["_cell_length_b"][0], "2.46772428");
    EXPECT_EQ(data["_cell_length_c"][0], "8.68503800");
    EXPECT_EQ(data["_cell_angle_alpha"][0], "90.00000000");
    EXPECT_EQ(data["_cell_angle_beta"][0], "90.00000000");
    EXPECT_EQ(data["_cell_angle_gamma"][0], "120.00000758");
    EXPECT_EQ(data["_symmetry_Int_Tables_number"][0], "1");
    EXPECT_EQ(data["_chemical_formula_structural"][0], "C");
    EXPECT_EQ(data["_chemical_formula_sum"][0], "C4");
    EXPECT_EQ(data["_cell_volume"][0], "45.80317575");
    EXPECT_EQ(data["_cell_formula_units_Z"][0], "4");
    EXPECT_EQ(data["_symmetry_equiv_pos_site_id"][0], "1");
    EXPECT_EQ(data["_symmetry_equiv_pos_as_xyz"][0], "'x, y, z'");
    EXPECT_EQ(data["_atom_site_type_symbol"][0], "C");
    EXPECT_EQ(data["_atom_site_label"][0], "C0");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][0], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][0], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_y"][0], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_z"][0], "0.75000000");
    EXPECT_EQ(data["_atom_site_occupancy"][0], "1");
    EXPECT_EQ(data["_atom_site_type_symbol"][1], "C");
    EXPECT_EQ(data["_atom_site_label"][1], "C1");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][1], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][1], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_y"][1], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_z"][1], "0.25000000");
    EXPECT_EQ(data["_atom_site_occupancy"][1], "1");
    EXPECT_EQ(data["_atom_site_type_symbol"][2], "C");
    EXPECT_EQ(data["_atom_site_label"][2], "C2");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][2], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][2], "0.33333300");
    EXPECT_EQ(data["_atom_site_fract_y"][2], "0.66666700");
    EXPECT_EQ(data["_atom_site_fract_z"][2], "0.75000000");
    EXPECT_EQ(data["_atom_site_occupancy"][2], "1");
    EXPECT_EQ(data["_atom_site_type_symbol"][3], "C");
    EXPECT_EQ(data["_atom_site_label"][3], "C3");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][3], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][3], "0.66666700");
    EXPECT_EQ(data["_atom_site_fract_y"][3], "0.33333300");
    EXPECT_EQ(data["_atom_site_fract_z"][3], "0.25000000");
    EXPECT_EQ(data["_atom_site_occupancy"][3], "1");
}

TEST(CifParserTest, ReadMediumTest)
{
    int rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef __MPI
    if (rank == 0)
    {
#endif
    std::ofstream ofs("cod-1000065.cif");
    ofs << cod1000065;
    ofs.close();
#ifdef __MPI
    }
    MPI_Barrier(MPI_COMM_WORLD); // make sure the file is written
#endif
    std::map<std::string, std::vector<std::string>> data;
    ModuleIO::CifParser::read("cod-1000065.cif", data);
    // delete the file
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD); // make sure the file is already read
    if (rank == 0)
    {
#endif
    std::remove("cod-1000065.cif");
#ifdef __MPI
    }
#endif

    EXPECT_EQ(data.size(), 43);
    EXPECT_EQ(data["_publ_author_name"][0], "'Nixon, D E' 'Parry, G S' 'Ubbelohde, A R'");
    EXPECT_EQ(data["_publ_section_title"][0], "; Order-disorder transformations in graphite nitrates ;");
    EXPECT_EQ(data["_journal_coden_ASTM"][0], "PRLAAZ");
    EXPECT_EQ(data["_journal_name_full"][0], "; Proceedings of the Royal Society of London, Series A: Mathematical and Physical Sciences (76,1906-) ;");
    EXPECT_EQ(data["_journal_page_first"][0], "324");
    EXPECT_EQ(data["_journal_page_last"][0], "339");
    EXPECT_EQ(data["_journal_paper_doi"][0], "10.1098/rspa.1966.0098");
    EXPECT_EQ(data["_journal_volume"][0], "291");
    EXPECT_EQ(data["_journal_year"][0], "1966");
    EXPECT_EQ(data["_chemical_formula_analytical"][0], "'C (H N O3)'");
    EXPECT_EQ(data["_chemical_formula_structural"][0], "C");
    EXPECT_EQ(data["_chemical_formula_sum"][0], "C");
    EXPECT_EQ(data["_chemical_name_common"][0], "'Graphite nitrate'");
    EXPECT_EQ(data["_chemical_name_systematic"][0], "Carbon");
    EXPECT_EQ(data["_space_group_IT_number"][0], "166");
    EXPECT_EQ(data["_symmetry_cell_setting"][0], "trigonal");
    EXPECT_EQ(data["_symmetry_space_group_name_Hall"][0], "'-R 3 2\"'");
    EXPECT_EQ(data["_symmetry_space_group_name_H-M"][0], "'R -3 m :H'");
    EXPECT_EQ(data["_cell_angle_alpha"][0], "90");
    EXPECT_EQ(data["_cell_angle_beta"][0], "90");
    EXPECT_EQ(data["_cell_angle_gamma"][0], "120");
    EXPECT_EQ(data["_cell_formula_units_Z"][0], "12");
    EXPECT_EQ(data["_cell_length_a"][0], "2.46");
    EXPECT_EQ(data["_cell_length_b"][0], "2.46");
    EXPECT_EQ(data["_cell_length_c"][0], "33.45");
    EXPECT_EQ(data["_cell_volume"][0], "175.3");
    EXPECT_EQ(data["_cod_original_sg_symbol_H-M"][0], "'R -3 m H'");
    EXPECT_EQ(data["_cod_database_code"][0], "1000065");
    EXPECT_EQ(data["_symmetry_equiv_pos_as_xyz"][0], "x,y,z -y,x-y,z y-x,-x,z -y,-x,z x,x-y,z y-x,y,z -x,-y,-z y,y-x,-z x-y,x,-z y,x,-z -x,y-x,-z x-y,-y,-z 1/3+x,2/3+y,2/3+z 2/3+x,1/3+y,1/3+z 1/3-y,2/3+x-y,2/3+z 2/3-y,1/3+x-y,1/3+z 1/3-x+y,2/3-x,2/3+z 2/3-x+y,1/3-x,1/3+z 1/3-y,2/3-x,2/3+z 2/3-y,1/3-x,1/3+z 1/3+x,2/3+x-y,2/3+z 2/3+x,1/3+x-y,1/3+z 1/3-x+y,2/3+y,2/3+z 2/3-x+y,1/3+y,1/3+z 1/3-x,2/3-y,2/3-z 2/3-x,1/3-y,1/3-z 1/3+y,2/3-x+y,2/3-z 2/3+y,1/3-x+y,1/3-z 1/3+x-y,2/3+x,2/3-z 2/3+x-y,1/3+x,1/3-z 1/3+y,2/3+x,2/3-z 2/3+y,1/3+x,1/3-z 1/3-x,2/3-x+y,2/3-z 2/3-x,1/3-x+y,1/3-z 1/3+x-y,2/3-y,2/3-z 2/3+x-y,1/3-y,1/3-z");
    EXPECT_EQ(data["_atom_site_label"][0], "C1");
    EXPECT_EQ(data["_atom_site_type_symbol"][0], "C0");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][0], "6");
    EXPECT_EQ(data["_atom_site_Wyckoff_symbol"][0], "c");
    EXPECT_EQ(data["_atom_site_fract_x"][0], "0.");
    EXPECT_EQ(data["_atom_site_fract_y"][0], "0.");
    EXPECT_EQ(data["_atom_site_fract_z"][0], "0.05");
    EXPECT_EQ(data["_atom_site_occupancy"][0], "1.");
    EXPECT_EQ(data["_atom_site_attached_hydrogens"][0], "0");
    EXPECT_EQ(data["_atom_site_calc_flag"][0], "d");
    EXPECT_EQ(data["_atom_site_label"][1], "C2");
    EXPECT_EQ(data["_atom_site_type_symbol"][1], "C0");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][1], "6");
    EXPECT_EQ(data["_atom_site_Wyckoff_symbol"][1], "c");
    EXPECT_EQ(data["_atom_site_fract_x"][1], "0.");
    EXPECT_EQ(data["_atom_site_fract_y"][1], "0.");
    EXPECT_EQ(data["_atom_site_fract_z"][1], "0.283");
    EXPECT_EQ(data["_atom_site_occupancy"][1], "1.");
    EXPECT_EQ(data["_atom_site_attached_hydrogens"][1], "0");
    EXPECT_EQ(data["_atom_site_calc_flag"][1], "d");
    EXPECT_EQ(data["_atom_type_symbol"][0], "C0");
    EXPECT_EQ(data["_atom_type_oxidation_number"][0], "0.000");
}
// because it is relatively hard to define loop_ by ABACUS itself, here the cooperative test
// will be performed by write-read manner.
TEST(CifParserTest, WriteTest)
{
    int rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    const std::string fcif = "test.cif";
    const std::vector<double> abc_angles = {2.46637620, 2.46637620, 24.84784531, 90.0, 90.0, 120.0};
    const int natom = 4;
    const std::vector<std::string> atom_site_labels = {"C", "C", "C", "C"};
    const std::vector<double> atom_site_fract = {0.0, 0.0, 0.75, 
                                                 0.0, 0.0, 0.25, 
                                                 0.333333, 0.666667, 0.75, 
                                                 0.666667, 0.333333, 0.25};
    ModuleIO::CifParser::write(fcif, 
                               abc_angles.data(), 
                               natom, 
                               atom_site_labels.data(), 
                               atom_site_fract.data(),
                               "# Generated during unittest of function ModuleIO::CifParser::write",
                               "data_test",
                               rank);
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD); // make sure the file is written
#endif
    std::map<std::string, std::vector<std::string>> data;
    ModuleIO::CifParser::read(fcif, data, rank);
    // delete the file
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD); // make sure the file is already read
    if (rank == 0)
    {
#endif
    std::remove(fcif.c_str());
#ifdef __MPI
    }
#endif

    EXPECT_EQ(data.size(), 23);
    EXPECT_EQ(data["_symmetry_space_group_name_H-M"][0], "'P 1'");
    EXPECT_EQ(data["_cell_length_a"][0], "2.46637620");
    EXPECT_EQ(data["_cell_length_b"][0], "2.46637620");
    EXPECT_EQ(data["_cell_length_c"][0], "24.84784531");
    EXPECT_EQ(data["_cell_angle_alpha"][0], "90.00000000");
    EXPECT_EQ(data["_cell_angle_beta"][0], "90.00000000");
    EXPECT_EQ(data["_cell_angle_gamma"][0], "120.00000000");
    EXPECT_EQ(data["_symmetry_Int_Tables_number"][0], "1");
    EXPECT_EQ(data["_chemical_formula_structural"][0], "C");
    EXPECT_EQ(data["_chemical_formula_sum"][0], "C4");
    EXPECT_EQ(data["_cell_volume"][0], "130.89950618");
    EXPECT_EQ(data["_cell_formula_units_Z"][0], "1");
    EXPECT_EQ(data["_symmetry_equiv_pos_site_id"][0], "1");
    EXPECT_EQ(data["_symmetry_equiv_pos_as_xyz"][0], "'x, y, z'");
    EXPECT_EQ(data["_atom_site_type_symbol"][0], "C");
    EXPECT_EQ(data["_atom_site_label"][0], "C0");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][0], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][0], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_y"][0], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_z"][0], "0.75000000");
    EXPECT_EQ(data["_atom_site_occupancy"][0], "1.0");
    EXPECT_EQ(data["_atom_site_type_symbol"][1], "C");
    EXPECT_EQ(data["_atom_site_label"][1], "C1");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][1], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][1], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_y"][1], "0.00000000");
    EXPECT_EQ(data["_atom_site_fract_z"][1], "0.25000000");
    EXPECT_EQ(data["_atom_site_occupancy"][1], "1.0");
    EXPECT_EQ(data["_atom_site_type_symbol"][2], "C");
    EXPECT_EQ(data["_atom_site_label"][2], "C2");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][2], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][2], "0.33333300");
    EXPECT_EQ(data["_atom_site_fract_y"][2], "0.66666700");
    EXPECT_EQ(data["_atom_site_fract_z"][2], "0.75000000");
    EXPECT_EQ(data["_atom_site_occupancy"][2], "1.0");
    EXPECT_EQ(data["_atom_site_type_symbol"][3], "C");
    EXPECT_EQ(data["_atom_site_label"][3], "C3");
    EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][3], "1");
    EXPECT_EQ(data["_atom_site_fract_x"][3], "0.66666700");
    EXPECT_EQ(data["_atom_site_fract_y"][3], "0.33333300");
    EXPECT_EQ(data["_atom_site_fract_z"][3], "0.25000000");
    EXPECT_EQ(data["_atom_site_occupancy"][3], "1.0");
}


int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return ret;
}