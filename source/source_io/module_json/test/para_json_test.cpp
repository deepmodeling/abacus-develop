#include "gtest/gtest.h"

#include <nlohmann/json.hpp>

#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_global.h"
#include "source_base/vector3.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/unitcell.h"
#include "source_io/module_json/abacusjson.h"
#include "source_io/module_json/general_info.h"
#include "source_io/module_json/init_info.h"
#include "source_io/module_json/output_info.h"
#include "source_io/module_parameter/parameter.h"
#include "source_main/version.h"

class AbacusJsonTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        Json::AbacusJson::document() = Json::jsonValue::object();
    }

    void TearDown() override
    {
        std::remove("test.json");
        std::remove("json-output-not-a-directory");
    }

    const Json::jsonValue& document() const
    {
        return Json::AbacusJson::document();
    }
};

TEST_F(AbacusJsonTest, OutputJson)
{
    // Exercise our writer, including escaping, number types and insertion order.
    Json::AbacusJson::document() = {
        {"z", "quote: \"; slash: \\; newline: \n; UTF-8: \xCE\xB1"},
        {"a", std::string("a\0b", 3)},
        {"nested", {{"int", 1}, {"float", 0.1}, {"bool", true}, {"null", nullptr},
                    {"array", Json::jsonValue::array({1, 2, 3})}}}};
    Json::AbacusJson::write_to_json("test.json");
    std::ifstream file("test.json");
    ASSERT_TRUE(file.is_open());
    const Json::jsonValue result = Json::jsonValue::parse(file);
    EXPECT_EQ(result, document());
    EXPECT_EQ(result.dump(), document().dump());
}

TEST_F(AbacusJsonTest, GeneralInfo)
{
    Parameter param;
    Json::gen_general_info(param);

    const Json::jsonValue& info = document().at("general_info");
    EXPECT_EQ(info["version"], VERSION);
    EXPECT_EQ(info["device"], param.inp.device);
#ifdef __MPI
    EXPECT_EQ(info["mpi_num"], Parallel_Global::mpi_number);
    EXPECT_EQ(info["omp_num"], Parallel_Global::omp_number);
#else
    EXPECT_EQ(info["mpi_num"], 1);
    EXPECT_EQ(info["omp_num"], 1);
#endif
    EXPECT_EQ(info["orbital_dir"], param.inp.orbital_dir);
    EXPECT_EQ(info["pseudo_dir"], param.inp.pseudo_dir);
    EXPECT_EQ(info["stru_file"], param.globalv.global_in_stru);
    EXPECT_EQ(info["kpt_file"], param.inp.kpoint_file);
    EXPECT_TRUE(info["start_time"].is_string());
    EXPECT_TRUE(info["end_time"].is_string());
    std::vector<std::string> keys;
    for (Json::jsonValue::const_iterator field = info.begin(); field != info.end(); ++field)
    {
        keys.push_back(field.key());
    }
    EXPECT_EQ(keys, (std::vector<std::string>{"version", "commit", "device", "mpi_num", "omp_num",
                                            "pseudo_dir", "orbital_dir", "stru_file", "kpt_file",
                                            "start_time", "end_time"}));
    Json::add_nkstot(2);
    Json::gen_general_info(param);
    EXPECT_EQ(document()["init"]["nkstot"], 2);
    EXPECT_EQ(document()["general_info"].size(), keys.size());
}

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}

Magnetism::~Magnetism()
{
}

TEST_F(AbacusJsonTest, InitInfo)
{
    UnitCell ucell;
    Atom atomlist[3];

    ucell.symm.pgname = "T_d";
    ucell.symm.spgname = "O_h";
    ucell.atoms = atomlist;
    ucell.ntype = 3;

    Input_para inp;
    inp.nbands = 10;
    inp.ecutwfc = 50.0;
    inp.smearing_method = "gauss";
    inp.smearing_sigma = 0.015;
    inp.kspacing = {0.04, 0.04, 0.04};
    inp.koffset = {0.0, 0.0, 0.0};
    inp.kmesh_type = "gamma";

    ucell.atoms[0].label = "Si";
    ucell.atoms[0].ncpp.zv = 3;
    ucell.atoms[0].na = 1;
    ucell.atoms[1].label = "C";
    ucell.atoms[1].ncpp.zv = 4;
    ucell.atoms[1].na = 2;
    ucell.atoms[2].label = "O";
    ucell.atoms[2].ncpp.zv = 5;
    ucell.atoms[2].na = 3;

    ucell.nat = 6;

    Json::add_nkstot(1);
    Json::gen_init(&ucell, inp);

    const Json::jsonValue& init = document().at("init");
    EXPECT_EQ(init["nkstot"], 1);
    EXPECT_EQ(init["natom"], 6);
    EXPECT_EQ(init["nband"], 10);
    EXPECT_EQ(init["point_group"], "T_d");
    EXPECT_EQ(init["point_group_in_space"], "O_h");
    EXPECT_EQ(init.at("nelectron_each_type"), (Json::jsonValue{{"Si", 3.0}, {"C", 4.0}, {"O", 5.0}}));
    EXPECT_EQ(init.at("natom_each_type"), (Json::jsonValue{{"Si", 1}, {"C", 2}, {"O", 3}}));
    EXPECT_EQ(init.at("nelectron"), 26);
    EXPECT_TRUE(init.at("nelectron").is_number_integer());
    EXPECT_TRUE(init.at("nelectron_each_type").at("C").is_number_float());
    EXPECT_EQ(init["ecutwfc"], 50.0);
    EXPECT_EQ(init["ecutwfc_unit"], "Ry");
    EXPECT_EQ(init["smearing_method"], "gauss");
    EXPECT_EQ(init["smearing_sigma"], 0.015);
    EXPECT_EQ(init["smearing_sigma_unit"], "Ry");
    EXPECT_EQ(init["kmesh_type"], "gamma");
    EXPECT_EQ(init["kspacing"], Json::jsonValue::array({0.04, 0.04, 0.04}));
    EXPECT_EQ(init["koffset"], Json::jsonValue::array({0.0, 0.0, 0.0}));

    // Rebuild the per-species maps rather than retaining entries from a previous call.
    ucell.ntype = 2;
    ucell.nat = 3;
    Json::gen_init(&ucell, inp);
    EXPECT_EQ(init.at("natom_each_type"), (Json::jsonValue{{"Si", 1}, {"C", 2}}));
    EXPECT_EQ(init.at("nelectron_each_type"), (Json::jsonValue{{"Si", 3.0}, {"C", 4.0}}));
    EXPECT_EQ(init.at("nkstot"), 1);
}

TEST_F(AbacusJsonTest, InitStructure)
{
    UnitCell ucell;
    Atom atom;

    ucell.latvec = ModuleBase::Matrix3(0.1, 0.1, 0.1,
                                       0.2, 0.2, 0.2,
                                       0.3, 0.3, 0.3);
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.pseudo_fn = {"si.ufp"};
    ucell.orbital_fn = {""};
    ucell.atoms = &atom;
    ucell.lat0 = 10.0;
    ucell.lat0_angstrom = ucell.lat0 * ModuleBase::BOHR_TO_A;

    atom.label = "Si";
    atom.na = 2;
    atom.ncpp.zv = 4.0;
    atom.tau = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                ModuleBase::Vector3<double>(0.1, 0.1, 0.1)};
    atom.mag = {0.0, 131.0};

    Input_para inp;
    Json::gen_stru(&ucell, inp);

    const Json::jsonValue& init = document().at("init");
    EXPECT_EQ(init["mag"], Json::jsonValue::array({0.0, 131.0}));
    EXPECT_EQ(init["pp"]["Si"], "si.ufp");
    EXPECT_TRUE(init["orb"]["Si"].is_null());
    EXPECT_EQ(init["label"][0], "Si");
    EXPECT_EQ(init["element"]["Si"], "");
    EXPECT_EQ(init["coordinate"][0], Json::jsonValue::array({0.0, 0.0, 0.0}));
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(init["coordinate"][1][i].get<double>(), ModuleBase::BOHR_TO_A, 1.0e-12);
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(init["cell"][i][j].get<double>(), (i + 1) * ModuleBase::BOHR_TO_A, 1.0e-12);
        }
    }

    // Reuse the cell to check shared init fields and repeated generation.
    ucell.orbital_fn[0] = "Si.orb";
    inp.orbital_dir = "orbitals/";
    inp.kspacing = {0.1, 0.2, 0.3};
    inp.koffset = {0.0, 0.5, 0.0};
    Json::gen_stru(&ucell, inp);
    Json::add_nkstot(3);
    Json::gen_init(&ucell, inp);
    const Json::jsonValue first = document();
    EXPECT_EQ(first["init"]["orb"]["Si"], "orbitals/Si.orb");
    EXPECT_EQ(first["init"]["nkstot"], 3);
    EXPECT_EQ(first["init"]["natom"], 2);
    Json::gen_init(&ucell, inp);
    Json::gen_stru(&ucell, inp);
    EXPECT_EQ(document(), first);
    EXPECT_EQ(document().dump(), first.dump()); // Preserve key order, too.
}

TEST_F(AbacusJsonTest, OutputRecords)
{
    EXPECT_THROW(Json::add_output_energy(-1.0), std::invalid_argument);
    EXPECT_FALSE(document().contains("output"));
    Json::AbacusJson::document()["output"] = Json::jsonValue::array();
    EXPECT_THROW(Json::add_output_energy(-1.0), std::out_of_range);
    Json::init_output_array_obj();
    ASSERT_EQ(document().at("output").size(), 1u);
    const Json::jsonValue initial = document()["output"][0];
    for (const char* key : {"e_fermi", "energy", "scf_converge", "force", "stress"})
    {
        EXPECT_TRUE(initial[key].is_null());
    }
    for (const char* key : {"coordinate", "mag", "cell"})
    {
        EXPECT_TRUE(initial[key].is_array());
        EXPECT_TRUE(initial[key].empty());
    }

    Json::add_output_efermi_converge(1.5, true);
    Json::add_output_energy(-10.0);
    Json::add_output_scf_mag(1.0, 2.0, -9.0, -0.2, 1.0e-3, 0.5);
    Json::add_output_scf_mag(1.0, 2.0, -10.0, -1.0, 1.0e-5, 0.6);

    const Json::jsonValue first = document()["output"][0];
    EXPECT_EQ(first["e_fermi"], 1.5);
    EXPECT_EQ(first["energy"], -10.0);
    EXPECT_EQ(first["scf_converge"], true);
    EXPECT_EQ(first["total_mag"], 1.0);
    EXPECT_EQ(first["absolute_mag"], 2.0);
    ASSERT_EQ(first["scf"].size(), 2u);
    EXPECT_EQ(first["scf"][1]["energy"], -10.0);
    EXPECT_EQ(first["scf"][1]["ediff"], -1.0);
    EXPECT_EQ(first["scf"][1]["drho"], 1.0e-5);
    EXPECT_EQ(first["scf"][1]["time"], 0.6);

    Json::init_output_array_obj();
    Json::add_output_energy(-11.0);
    ASSERT_EQ(document()["output"].size(), 2u);
    EXPECT_EQ(document()["output"][0], first);
    EXPECT_EQ(document()["output"][1]["energy"], -11.0);
}

TEST_F(AbacusJsonTest, RejectsInvalidSections)
{
    for (const Json::jsonValue& invalid : {Json::jsonValue(nullptr), Json::jsonValue(1),
                                           Json::jsonValue("invalid"), Json::jsonValue::array()})
    {
        Json::AbacusJson::document()["init"] = invalid;
        EXPECT_THROW(Json::add_nkstot(1), std::invalid_argument);
        EXPECT_EQ(document().at("init"), invalid);
    }
    for (const Json::jsonValue& invalid : {Json::jsonValue(nullptr), Json::jsonValue(1),
                                           Json::jsonValue("invalid"), Json::jsonValue::object()})
    {
        Json::AbacusJson::document()["output"] = invalid;
        EXPECT_THROW(Json::init_output_array_obj(), std::invalid_argument);
        EXPECT_THROW(Json::add_output_energy(-1.0), std::invalid_argument);
        EXPECT_EQ(document().at("output"), invalid);
    }
}

TEST_F(AbacusJsonTest, RejectsInvalidRecordsAndScfHistory)
{
    for (const Json::jsonValue& invalid : {Json::jsonValue(nullptr), Json::jsonValue(1),
                                           Json::jsonValue::array()})
    {
        Json::AbacusJson::document()["output"] = Json::jsonValue::array({invalid});
        EXPECT_THROW(Json::add_output_energy(-1.0), std::invalid_argument);
        EXPECT_EQ(document().at("output").at(0), invalid);
    }
    Json::init_output_array_obj();
    for (const Json::jsonValue& invalid : {Json::jsonValue(nullptr), Json::jsonValue(1),
                                           Json::jsonValue::object()})
    {
        Json::AbacusJson::document()["output"].back()["scf"] = invalid;
        EXPECT_THROW(Json::add_output_scf_mag(0.0, 0.0, -1.0, 0.0, 0.1, 0.1), std::invalid_argument);
        EXPECT_EQ(document().at("output").back().at("scf"), invalid);
    }
}

TEST_F(AbacusJsonTest, OutputStructureForceAndStress)
{
    UnitCell ucell;
    Atom atom;
    ucell.atoms = &atom;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.lat0_angstrom = 2.0;
    ucell.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                       0.0, 2.0, 0.0,
                                       0.0, 0.0, 3.0);
    atom.na = 1;
    atom.tau = {ModuleBase::Vector3<double>(0.25, -0.5, 0.75)};
    atom.mag = {1.5};

    ModuleBase::matrix force(1, 3);
    force(0, 0) = 0.5e-8;
    force(0, 1) = 2.0;
    force(0, 2) = -3.0;

    ModuleBase::matrix stress(3, 3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            stress(i, j) = 3 * i + j + 1;
        }
    }

    Json::init_output_array_obj();
    Json::add_output_cell_coo_stress_force(ucell, force, 2.0, stress, 0.5, true, true);
    const Json::jsonValue first = document()["output"][0];
    ASSERT_EQ(first["force"].size(), 1u);
    EXPECT_EQ(first["force"][0], Json::jsonValue::array({0.0, 4.0, -6.0}));
    EXPECT_EQ(first["coordinate"][0], Json::jsonValue::array({0.5, -1.0, 1.5}));
    EXPECT_EQ(first["mag"], Json::jsonValue::array({1.5}));
    ASSERT_EQ(first["stress"].size(), 3u);
    ASSERT_EQ(first["cell"].size(), 3u);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_DOUBLE_EQ(first["stress"][i][j].get<double>(), stress(i, j) * 0.5);
            EXPECT_DOUBLE_EQ(first["cell"][i][j].get<double>(), i == j ? 2.0 * (i + 1) : 0.0);
        }
    }

    // Replacing the same step's data must not append extra rows or nested arrays.
    Json::add_output_cell_coo_stress_force(ucell, force, 2.0, stress, 0.5, true, true);
    EXPECT_EQ(document()["output"][0], first);

    Json::init_output_array_obj();
    Json::add_output_cell_coo_stress_force(ucell, force, 2.0, stress, 0.5, false, false);
    EXPECT_EQ(document()["output"][0], first);
    const Json::jsonValue& second = document()["output"][1];
    EXPECT_TRUE(second["force"].is_null());
    EXPECT_TRUE(second["stress"].is_null());
    EXPECT_EQ(second["coordinate"], first["coordinate"]);
    EXPECT_EQ(second["cell"], first["cell"]);
}

TEST_F(AbacusJsonTest, NonFiniteNumbersSerializeAsNull)
{
    Json::init_output_array_obj();
    Json::add_output_energy(std::numeric_limits<double>::quiet_NaN());
    Json::add_output_efermi_converge(std::numeric_limits<double>::infinity(), false);
    Json::AbacusJson::write_to_json("test.json");
    std::ifstream file("test.json");
    ASSERT_TRUE(file.is_open());
    const Json::jsonValue result = Json::jsonValue::parse(file);
    EXPECT_TRUE(result.at("output").at(0).at("energy").is_null());
    EXPECT_TRUE(result.at("output").at(0).at("e_fermi").is_null());
}

TEST_F(AbacusJsonTest, FileOpenFailureIsReported)
{
    const std::string blocker = "json-output-not-a-directory";
    {
        std::ofstream file(blocker);
        ASSERT_TRUE(file.is_open());
    }
    EXPECT_THROW(Json::AbacusJson::write_to_json(blocker + "/abacus.json"), std::runtime_error);
    EXPECT_EQ(std::remove(blocker.c_str()), 0);
}
