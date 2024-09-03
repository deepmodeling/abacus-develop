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


// TEST(CifParserTest, VecABCAnglesTest)
// {
//     // this test will test the conversion between vector and abc angles
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_real_distribution<double> dis(0.0, 1.0);
//     std::vector<double> vec(9);
//     // perform test for 10 times
//     for (int it = 0; it < 10; ++it)
//     {
//         std::for_each(vec.begin(), vec.end(), [&dis, &gen](double& d) { d = dis(gen); });
//         std::vector<double> abc_angles(6); // norma, normb, normc, alpha, beta, gamma
//         ModuleIO::CifParser::vec_to_abc_angles(vec.data(), abc_angles.data());
//         std::vector<double> vec_out(9);
//         ModuleIO::CifParser::abc_angles_to_vec(abc_angles.data(), vec_out.data());
//         // however, the function above will assume a as (a, 0, 0), so we convert again.
//         std::vector<double> abc_angles_out(6);
//         ModuleIO::CifParser::vec_to_abc_angles(vec_out.data(), abc_angles_out.data());
//         for (int i = 0; i < 6; ++i)
//         {
//             EXPECT_NEAR(abc_angles[i], abc_angles_out[i], 1e-6);
//         }
//     }
//     // above test the consistency between vec_to_abc_angles and abc_angles_to_vec
//     // the following test the correctness of vec_to_abc_angles. If the conversion is correct,
//     // then both are correct.
//     vec = {0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5}; // FCC
//     std::vector<double> abc_angles(6);
//     ModuleIO::CifParser::vec_to_abc_angles(vec.data(), abc_angles.data());
//     // the norm of a, b, c should be 0.5*sqrt(2)
//     EXPECT_NEAR(abc_angles[0], 0.5 * std::sqrt(2.0), 1e-6);
//     EXPECT_NEAR(abc_angles[1], 0.5 * std::sqrt(2.0), 1e-6);
//     EXPECT_NEAR(abc_angles[2], 0.5 * std::sqrt(2.0), 1e-6);
//     // the angle between a and b should be 60 degree
//     EXPECT_NEAR(abc_angles[3], 60.0, 1e-6);
//     // the angle between a and c should be 60 degree
//     EXPECT_NEAR(abc_angles[4], 60.0, 1e-6);
//     // the angle between b and c should be 60 degree
//     EXPECT_NEAR(abc_angles[5], 60.0, 1e-6);
// }

// TEST(CifParserTest, ToVolumeTest)
// {
//     // here we test two pratical system.
//     // the first is cif file Materials Project, the mp-2516584
//     std::vector<double> abc_angles = {2.46637620, 2.46637620, 24.84784531, 90.0, 90.0, 120.0};
//     double volume = ModuleIO::CifParser::abc_angles_to_volume(abc_angles.data());
//     EXPECT_NEAR(volume, 130.89950564, 1e-6);
//     // the second is also from Materials Project, the mp-1274279
//     abc_angles = {3.10229376, 3.10695673, 5.31894279, 74.11536155, 73.57430467, 61.15267638};
//     volume = ModuleIO::CifParser::abc_angles_to_volume(abc_angles.data());
//     EXPECT_NEAR(volume, 42.49423842, 1e-6);
// }

// TEST(CifParserTest, BuildChemFormulaTest)
// {
//     // this test will test the function _build_chem_formula
//     std::string sum, structural;
//     std::string atom_site_labels[4] = {"C", "C", "C", "C"};
//     ModuleIO::CifParser::_build_chem_formula(4, atom_site_labels, sum, structural);
//     EXPECT_EQ(sum, "C4");
//     EXPECT_EQ(structural, "C");

//     // the second test is for a more complicated system
//     std::string atom_site_labels2[4] = {"C", "C", "C", "O"};
//     ModuleIO::CifParser::_build_chem_formula(4, atom_site_labels2, sum, structural);
//     EXPECT_EQ(sum, "C3O");
//     EXPECT_EQ(structural, "CO");
// }

// TEST(CifParserTest, SplitOutsideEncloseTest)
// {
//     // test the function that can avoid split inside '
//     std::string in = "'P 1'";
//     std::vector<std::string> out = ModuleIO::CifParser::_split_outside_enclose(in, " ", {"'", "'"});
//     EXPECT_EQ(out.size(), 1);
//     EXPECT_EQ(out[0], "'P 1'");
//     in = "C1 C2 'C3 C4' C5 'C6 C7' C8";
//     out = ModuleIO::CifParser::_split_outside_enclose(in, " ", {"'", "'"});
//     EXPECT_EQ(out.size(), 6);
//     EXPECT_EQ(out[0], "C1");
//     EXPECT_EQ(out[1], "C2");
//     EXPECT_EQ(out[2], "'C3 C4'");
//     EXPECT_EQ(out[3], "C5");
//     EXPECT_EQ(out[4], "'C6 C7'");
//     EXPECT_EQ(out[5], "C8");

// }

// TEST(CifParserTest, FromCifBasicUtilsTest)
// {
//     // this test will test the function _split_loop_block
//     const std::string block1 = ""
//     "_symmetry_space_group_name_H-M   'P 1'\n"
//     "_cell_length_a   2.46772428\n"
//     "_cell_length_b   2.46772428\n"
//     "_cell_length_c   8.68503800\n"
//     "_cell_angle_alpha   90.00000000\n"
//     "_cell_angle_beta   90.00000000\n"
//     "_cell_angle_gamma   120.00000758\n"
//     "_symmetry_Int_Tables_number   1\n"
//     "_chemical_formula_structural   C\n"
//     "_chemical_formula_sum   C4\n"
//     "_cell_volume   45.80317575\n"
//     "_cell_formula_units_Z   4\n";
//     std::vector<std::string> words = ModuleIO::CifParser::_split_loop_block(block1);
//     EXPECT_EQ(words.size(), 24);
//     std::map<std::string, std::vector<std::string>> data = ModuleIO::CifParser::_build_block_data(words);
//     EXPECT_EQ(data.size(), 12);
//     EXPECT_EQ(data["_symmetry_space_group_name_H-M"][0], "'P 1'");
//     EXPECT_EQ(data["_cell_length_a"][0], "2.46772428");
//     EXPECT_EQ(data["_cell_length_b"][0], "2.46772428");
//     EXPECT_EQ(data["_cell_length_c"][0], "8.68503800");
//     EXPECT_EQ(data["_cell_angle_alpha"][0], "90.00000000");
//     EXPECT_EQ(data["_cell_angle_beta"][0], "90.00000000");
//     EXPECT_EQ(data["_cell_angle_gamma"][0], "120.00000758");
//     EXPECT_EQ(data["_symmetry_Int_Tables_number"][0], "1");
//     EXPECT_EQ(data["_chemical_formula_structural"][0], "C");
//     EXPECT_EQ(data["_chemical_formula_sum"][0], "C4");
//     EXPECT_EQ(data["_cell_volume"][0], "45.80317575");
//     EXPECT_EQ(data["_cell_formula_units_Z"][0], "4");
    
//     // from Crystallography Open Database (COD), cod 4128192
//     const std::string block2 = ""
//     "_publ_author_name\n"
//     "'Galley, Shane S.'\n"
//     "'Pattenaude, Scott A.'\n"
//     "'Gaggioli, Carlo Alberto'\n"
//     "'Qiao, Yusen'\n"
//     "'Sperling, Joseph M.'\n"
//     "'Zeller, Matthias'\n"
//     "'Pakhira, Srimanta'\n"
//     "'Mendoza-Cortes, Jose L'\n"
//     "'Schelter, Eric J.'\n"
//     "'Albrecht-Schmitt, Thomas E'\n"
//     "'Gagliardi, Laura'\n"
//     "'Bart, Suzanne C.'\n"
//     "_publ_section_title\n"
//     ";\n"
//     " Synthesis and Characterization of Tris-chelate Complexes for\n"
//     " Understanding f-Orbital Bonding in Later Actinides.\n"
//     ";\n"
//     "_journal_issue                   6\n"
//     "_journal_name_full               'Journal of the American Chemical Society'\n";
//     words = ModuleIO::CifParser::_split_loop_block(block2);
//     EXPECT_EQ(words.size(), 8);

//     const std::string block3 = ""
//     " _symmetry_equiv_pos_site_id\n"
//     " _symmetry_equiv_pos_as_xyz\n"
//     "  1  'x, y, z'\n";
//     words = ModuleIO::CifParser::_split_loop_block(block3);
//     EXPECT_EQ(words.size(), 3);
//     data = ModuleIO::CifParser::_build_block_data(words);
//     EXPECT_EQ(data.size(), 2);
//     EXPECT_EQ(data["_symmetry_equiv_pos_site_id"][0], "1");
//     EXPECT_EQ(data["_symmetry_equiv_pos_as_xyz"][0], "'x, y, z'");

//     const std::string block4 = ""
//     "_atom_site_type_symbol\n"
//     "_atom_site_label\n"
//     "_atom_site_symmetry_multiplicity\n"
//     "_atom_site_fract_x\n"
//     "_atom_site_fract_y\n"
//     "_atom_site_fract_z\n"
//     "_atom_site_occupancy\n"
//     " Fe2+  Fe0  1  0.00000000  0.00000000  0.00000000  1\n"
//     " Fe2+  Fe1  1  0.50000000  0.50000000  0.50000000  1\n"
//     " O2-  O2  1  0.77376500  0.76754050  0.74832550  1\n"
//     " O2-  O3  1  0.22623500  0.23245950  0.25167450  1\n";
//     words = ModuleIO::CifParser::_split_loop_block(block4);
//     EXPECT_EQ(words.size(), 8);
//     data = ModuleIO::CifParser::_build_block_data(words);
//     EXPECT_EQ(data.size(), 7);
//     std::vector<std::string> col1ref = {"Fe2+", "Fe2+", "O2-", "O2-"};
//     std::vector<std::string> col2ref = {"Fe0", "Fe1", "O2", "O3"};
//     std::vector<std::string> col3ref = {"1", "1", "1", "1"};
//     std::vector<std::string> col4ref = {"0.00000000", "0.50000000", "0.77376500", "0.22623500"};
//     std::vector<std::string> col5ref = {"0.00000000", "0.50000000", "0.76754050", "0.23245950"};
//     std::vector<std::string> col6ref = {"0.00000000", "0.50000000", "0.74832550", "0.25167450"};
//     std::vector<std::string> col7ref = {"1", "1", "1", "1"};
//     for (int i = 0; i < 4; ++i)
//     {
//         EXPECT_EQ(data["_atom_site_type_symbol"][i], col1ref[i]);
//         EXPECT_EQ(data["_atom_site_label"][i], col2ref[i]);
//         EXPECT_EQ(data["_atom_site_symmetry_multiplicity"][i], col3ref[i]);
//         EXPECT_EQ(data["_atom_site_fract_x"][i], col4ref[i]);
//         EXPECT_EQ(data["_atom_site_fract_y"][i], col5ref[i]);
//         EXPECT_EQ(data["_atom_site_fract_z"][i], col6ref[i]);
//         EXPECT_EQ(data["_atom_site_occupancy"][i], col7ref[i]);
//     }
// }

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
#endif
    std::map<std::string, std::vector<std::string>> data;
    ModuleIO::CifParser::read("mp-2516584.cif", data, rank);
    // delete the file
#ifdef __MPI
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
#endif
    std::map<std::string, std::vector<std::string>> data;
    ModuleIO::CifParser::read("cod-1000065.cif", data, rank);
    // delete the file
#ifdef __MPI
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
    std::ofstream ofs(fcif);
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
    std::map<std::string, std::vector<std::string>> data;
    ModuleIO::CifParser::read(fcif, data, rank);
    // delete the file
#ifdef __MPI
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