#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "module_base/formatter.h"
#include "module_io/cif_io.h"
#include <regex>
#include <queue> // first push, first pop
#include <cassert>

void ModuleIO::CifParser::_build_chem_formula(const int natom, 
                  const std::string* atom_site_labels,
                  std::string& sum,
                  std::string& structural)
{
    std::vector<std::string> kinds;
    std::vector<std::string> labels(natom);
    std::copy(atom_site_labels, atom_site_labels + natom, labels.begin());
    for (int i = 0; i < natom; ++i)
    {
        if (std::find(kinds.begin(), kinds.end(), labels[i]) == kinds.end())
        {
            kinds.push_back(labels[i]);
        }
    }
    std::vector<size_t> counts(kinds.size());
    std::transform(kinds.begin(), kinds.end(), counts.begin(), [&labels](const std::string& kind) {
        return std::count(labels.begin(), labels.end(), kind);
    });
    for (size_t i = 0; i < kinds.size(); ++i)
    {
        sum += kinds[i];
        structural += kinds[i];
        if (counts[i] > 1)
        {
            sum += std::to_string(counts[i]);
        }
    }
}

void ModuleIO::CifParser::vec_to_abc_angles(const double* vec, double* abc_angles)
{
    const std::vector<double> a = {vec[0], vec[1], vec[2]};
    const std::vector<double> b = {vec[3], vec[4], vec[5]};
    const std::vector<double> c = {vec[6], vec[7], vec[8]};
    const double anorm = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    const double bnorm = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    const double cnorm = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    const double alpha = std::acos((b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) / (bnorm * cnorm));
    const double beta = std::acos((a[0] * c[0] + a[1] * c[1] + a[2] * c[2]) / (anorm * cnorm));
    const double gamma = std::acos((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (anorm * bnorm));
    abc_angles[0] = anorm;
    abc_angles[1] = bnorm;
    abc_angles[2] = cnorm;
    abc_angles[3] = rad2deg(alpha);
    abc_angles[4] = rad2deg(beta);
    abc_angles[5] = rad2deg(gamma);
}
void ModuleIO::CifParser::abc_angles_to_vec(const double* abc_angles, double* vec)
{
    const double a = abc_angles[0];
    const double b = abc_angles[1];
    const double c = abc_angles[2];
    const double alpha = abc_angles[3];
    const double beta = abc_angles[4];
    const double gamma = abc_angles[5];
    vec[0] = a;
    vec[1] = 0.0;
    vec[2] = 0.0;
    vec[3] = b * std::cos(deg2rad(gamma));
    vec[4] = b * std::sin(deg2rad(gamma));
    vec[5] = 0.0;
    vec[6] = c * std::cos(deg2rad(beta));
    vec[7] = c * (std::cos(deg2rad(alpha)) - std::cos(deg2rad(beta)) * std::cos(deg2rad(gamma))) / std::sin(deg2rad(gamma));
    vec[8] = std::sqrt(c * c - vec[6] * vec[6] - vec[7] * vec[7]);
}
double ModuleIO::CifParser::vec_to_volume(const double* vec)
{
    // vector's mixed product
    return vec[0] * (vec[4] * vec[8] - vec[5] * vec[7]) - vec[1] * (vec[3] * vec[8] - vec[5] * vec[6]) + vec[2] * (vec[3] * vec[7] - vec[4] * vec[6]);
}
double ModuleIO::CifParser::abc_angles_to_volume(const double* abc_angles)
{
    std::vector<double> vec(9);
    abc_angles_to_vec(abc_angles, vec.data());
    return vec_to_volume(vec.data());
}

void ModuleIO::CifParser::_unpack_ucell(const UnitCell& ucell,
                                       std::vector<double>& veca,
                                       std::vector<double>& vecb,
                                       std::vector<double>& vecc,
                                       int& natom,
                                       std::vector<std::string>& atom_site_labels,
                                       std::vector<double>& atom_site_fract_coords)
{
    const double bohr2angstrom = 0.52917721067;
    const double lat0 = ucell.lat.lat0;
    veca.resize(3);
    vecb.resize(3);
    vecc.resize(3);
    veca[0] = ucell.a1.x;
    veca[1] = ucell.a1.y;
    veca[2] = ucell.a1.z;
    vecb[0] = ucell.a2.x;
    vecb[1] = ucell.a2.y;
    vecb[2] = ucell.a2.z;
    vecc[0] = ucell.a3.x;
    vecc[1] = ucell.a3.y;
    vecc[2] = ucell.a3.z;
    std::for_each(veca.begin(), veca.end(), [lat0, bohr2angstrom](double& x) { x *= lat0 * bohr2angstrom; });
    std::for_each(vecb.begin(), vecb.end(), [lat0, bohr2angstrom](double& x) { x *= lat0 * bohr2angstrom; });
    std::for_each(vecc.begin(), vecc.end(), [lat0, bohr2angstrom](double& x) { x *= lat0 * bohr2angstrom; });
    natom = ucell.nat;
    atom_site_labels.resize(natom);
    atom_site_fract_coords.resize(3 * natom);
    for (int i = 0; i < natom; ++i)
    {
        atom_site_labels[i] = ucell.atoms[ucell.iat2it[i]].ncpp.psd; // this is purely the element symbol
        atom_site_fract_coords[3 * i] = ucell.atoms[ucell.iat2it[i]].taud[ucell.iat2ia[i]].x;
        atom_site_fract_coords[3 * i + 1] = ucell.atoms[ucell.iat2it[i]].taud[ucell.iat2ia[i]].y;
        atom_site_fract_coords[3 * i + 2] = ucell.atoms[ucell.iat2it[i]].taud[ucell.iat2ia[i]].z;
    }
}

void ModuleIO::CifParser::to_cif(const std::string& fcif,
                                 const double* abc_angles,
                                 const int natom,
                                 const std::string* atom_site_labels,
                                 const double* atom_site_fract_coords,
                                 const std::string& title,
                                 const std::string& data_tag,
                                 const double* atom_site_occups,
                                 const std::string& cell_formula_units_z)
{
    std::ofstream ofs(fcif);
    if (!ofs)
    {
        std::cerr << "Error: failed to open file " << fcif << std::endl;
        return;
    }
    ofs << title << std::endl;
    ofs << data_tag << std::endl;
    ofs << "_symmetry_space_group_name_H-M   'P 1'" << std::endl;
    ofs << "_cell_length_a   " << abc_angles[0] << std::endl;
    ofs << "_cell_length_b   " << abc_angles[1] << std::endl;
    ofs << "_cell_length_c   " << abc_angles[2] << std::endl;
    ofs << "_cell_angle_alpha   " << abc_angles[3] << std::endl;
    ofs << "_cell_angle_beta   " << abc_angles[4] << std::endl;
    ofs << "_cell_angle_gamma   " << abc_angles[5] << std::endl;
    ofs << "_symmetry_Int_Tables_number   1" << std::endl;
    std::string chem_sum, chem_structural;
    _build_chem_formula(natom, atom_site_labels, chem_sum, chem_structural);
    ofs << "_chemical_formula_structural   " << chem_structural << std::endl;
    ofs << "_chemical_formula_sum   " << chem_sum << std::endl;
    ofs << "_cell_volume   " << abc_angles_to_volume(abc_angles) << std::endl;
    ofs << "_cell_formula_units_Z   " << cell_formula_units_z << std::endl;
    ofs << "loop_" << std::endl;
    ofs << " _symmetry_equiv_pos_site_id" << std::endl;
    ofs << " _symmetry_equiv_pos_as_xyz" << std::endl;
    ofs << "  1  'x, y, z'" << std::endl;
    ofs << "loop_" << std::endl;
    ofs << " _atom_site_type_symbol" << std::endl;
    ofs << " _atom_site_label" << std::endl;
    ofs << " _atom_site_symmetry_multiplicity" << std::endl;
    ofs << " _atom_site_fract_x" << std::endl;
    ofs << " _atom_site_fract_y" << std::endl;
    ofs << " _atom_site_fract_z" << std::endl;
    ofs << " _atom_site_occupancy" << std::endl;
    std::vector<double> occups(natom, 1.0);
    if (atom_site_occups != nullptr)
    {// overwrite the default occupancies
        std::copy(atom_site_occups, atom_site_occups + natom, occups.begin());
    }
    // then output atomic information with format: %3s%4s%3d%12.8f%12.8f%12.8f%3d
    FmtCore fmt("%3s%4s%3d%12.8f%12.8f%12.8f%3d");
    int j = 0;
    std::string numbered_label;
    std::string cache;
    for (int i = 0; i < natom; ++i)
    {
        const std::string label = atom_site_labels[i];
        numbered_label = label + std::to_string(i);
        cache = fmt.format(label.c_str(), numbered_label.c_str(), 1, 
        atom_site_fract_coords[j], atom_site_fract_coords[j + 1], atom_site_fract_coords[j + 2], occups[i]);
        ofs << cache << std::endl;
        j += 3;
    }
    ofs << std::endl;
    ofs.close();
}

void ModuleIO::CifParser::to_cif(const std::string& fcif,
                                 const UnitCell& ucell,
                                 const std::string& title,
                                 const std::string& data_tag)
{
    std::vector<double> veca, vecb, vecc;
    int natom;
    std::vector<std::string> atom_site_labels;
    std::vector<double> atom_site_fract_coords;
    _unpack_ucell(ucell, veca, vecb, vecc, natom, atom_site_labels, atom_site_fract_coords);
    double abc_angles[6];
    vec_to_abc_angles(veca.data(), abc_angles);
    to_cif(fcif.c_str(), abc_angles, natom, atom_site_labels.data(), atom_site_fract_coords.data(), title, data_tag);
}

// reading cif is another hard (physically) and laborious work. The cif sometimes can be easily read line by line,
// sometimes word by word. The structure of cif file is, except the line startswith "#", all other lines can be split
// by blocks leading by "loop_", then in each "loop_", there are contents can be split by those keywords that 
// startswith "_". There is also another exception that, if there is no space between keywords, then their values will
// appear after all keywords are listed. In this case, all values actually form a table, which is needed to be 
// furtherly formatted (rows are memory-contiguous).
// Thus the reading strategy are, 
// 1. first split the file into blocks by "loop_"
// 2. in each block, split with words starting with "_"
// 3. scan the splited words

void ModuleIO::CifParser::from_cif(const std::string& fcif,
                                   std::map<std::string, std::vector<std::string>>& out)
{
    std::ifstream ifs(fcif);
    std::string cache; // first read all lines into cache
    while (ifs.good())
    {
        std::string line;
        std::getline(ifs, line);
        cache += line + " ";
    }
    std::vector<std::string> blocks = FmtCore::split(cache, "loop_");
    for (auto& block: blocks)
    {
        std::vector<std::string> words = _split_loop_block(block);
        std::map<std::string, std::vector<std::string>> data = _build_block_data(words);
        out.insert(data.begin(), data.end());
    }
    // then print
    for (auto& kv: out)
    {
        std::cout << kv.first << ": ";
        for (auto& s: kv.second)
        {
            std::cout << s << " ";
        }
        std::cout << std::endl;
    }
}

// the second step, for each block, split with words starting with "_"
std::vector<std::string> ModuleIO::CifParser::_split_loop_block(const std::string& block)
{
    std::vector<std::string> out;
    std::string word, cache;
    std::stringstream ss(block);
    while (ss.good())
    {
        ss >> word;
        if (FmtCore::startswith(word, "_"))
        {
            if (!cache.empty())
            {
                out.push_back(cache);
                cache.clear();
            }
            out.push_back(word);
        }
        else
        {
            cache += word + " ";
        }
    }
    return out;
}

std::map<std::string, std::vector<std::string>> ModuleIO::CifParser::_build_block_data(const std::vector<std::string>& block)
{
    // the strategy is, read key(s) and cache, then read value(s) and cache, till meet the next key(s).
    // If there are only one key is stored in cache, then the number of values should also be 1.
    // , if there are more than only key, then a table must be made, and the number of columns,
    // say the number of keys should be provided.
    std::map<std::string, std::vector<std::string>> out;
    std::vector<std::string> kcache;
    std::vector<std::string> vcache; // actually it should be only one value
    for (auto& str: block)
    {
        if (FmtCore::startswith(str, "_")) // it is a new key read
        {
            if (!kcache.empty()) // there are data should be stored into "out"
            {
                if (kcache.size() == 1) // only one key, then only one value
                {
                    out[kcache[0]] = vcache; // store
                }
                else // a more complicated case, a table should be made.
                {
                    // split the value with white space
                    assert(vcache.size() == 1);
                    std::vector<std::string> words = FmtCore::split(vcache[0], " ");
                    const size_t ncols = kcache.size();
                    assert(words.size() % ncols == 0);
                    const size_t nrows = words.size() / ncols;
                    for (size_t i = 0; i < ncols; i++)
                    {
                        std::vector<std::string> col(nrows);
                        for (size_t j = 0; j < nrows; j++)
                        {
                            col[j] = words[j * ncols + i];
                        }
                        out[kcache[i]] = col;
                    }
                }
                kcache.clear();
                vcache.clear();
            }
            kcache.push_back(str);
        }
        else
        {
            vcache.push_back(str);
        }
    }
    return out;
}