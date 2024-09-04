#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "module_base/formatter.h"
#include "module_io/cif_io.h"
#include <regex>
#include <cassert>
#include "module_base/tool_quit.h"
#ifdef __MPI
#include "module_base/parallel_common.h"
#endif

double deg2rad(double deg) { return deg * M_PI / 180.0; }
double rad2deg(double rad) { return rad * 180.0 / M_PI; }

void _build_chem_formula(const int natom, 
                         const std::string* atom_site_labels,
                         std::string& sum,
                         std::string& structural)
{
    sum.clear();
    structural.clear();
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

void vec_to_abc_angles(const double* vec, double* abc_angles)
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
void abc_angles_to_vec(const double* abc_angles, double* vec)
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
double vec_to_volume(const double* vec)
{
    // vector's mixed product
    return vec[0] * (vec[4] * vec[8] - vec[5] * vec[7]) - vec[1] * (vec[3] * vec[8] - vec[5] * vec[6]) + vec[2] * (vec[3] * vec[7] - vec[4] * vec[6]);
}
double abc_angles_to_volume(const double* abc_angles)
{
    std::vector<double> vec(9);
    abc_angles_to_vec(abc_angles, vec.data());
    return vec_to_volume(vec.data());
}
std::vector<std::string> _split_outside_enclose(const std::string& in, 
                                                const std::string& delim,
                                                const std::vector<std::string>& enclose)
{
    // a very naive impl. for only CIF possible cases
    assert(enclose.size() == 2); // other complicated case not implemented yet
    // first split with delim. then scan all fragments, if there are enclose symbol, then will first meet
    // a fragment startswith the opening, then after fragments, there will be a one ends with closing.
    // between them, fragments will be concatenated with delim.
    std::vector<std::string> out;
    std::string cache;
    std::vector<std::string> words = FmtCore::split(in, delim);
    bool in_enclose = false;
    for (auto& word: words)
    {
        if (FmtCore::startswith(word, enclose[0]))
        {
            in_enclose = true;
            cache += word + delim;
        }
        else if (FmtCore::endswith(word, enclose[1]))
        {
            in_enclose = false;
            cache += word;
            out.push_back(cache);
            cache.clear();
        }
        else
        {
            if (in_enclose)
            {
                cache += word + delim;
            }
            else
            {
                out.push_back(word);
            }
        }
    }
    return out;
}

// the second step, for each block, split with words starting with "_"
std::vector<std::string> _split_loop_block(const std::string& block)
{
    std::vector<std::string> out;
    std::string word, cache;
    std::stringstream ss(FmtCore::strip(FmtCore::strip(block, "\n")));
    while (ss.good())
    {
        ss >> word;
        if (FmtCore::startswith(word, "_"))
        {
            if (!cache.empty())
            {
                out.push_back(FmtCore::strip(cache));
                cache.clear();
            }
            out.push_back(FmtCore::strip(word));
        }
        else
        {
            cache += word + " ";
        }
    }
    // the last word
    if (!cache.empty())
    {
        out.push_back(FmtCore::strip(cache));
    }
    return out;
}

std::map<std::string, std::vector<std::string>> _build_table(const std::vector<std::string>& keys,
                                                             const std::vector<std::string>& values)
{
    std::map<std::string, std::vector<std::string>> out;
    const size_t ncols = keys.size();
    assert(values.size() % ncols == 0);
    const size_t nrows = values.size() / ncols;
    for (size_t i = 0; i < ncols; i++)
    {
        std::vector<std::string> col(nrows);
        for (size_t j = 0; j < nrows; j++)
        {
            col[j] = values[j * ncols + i];
        }
        out[keys[i]] = col;
    }
    return out;
}

std::map<std::string, std::vector<std::string>> _build_block_data(const std::vector<std::string>& block)
{
    // after calling the _split_loop_block, the data now composed of elements that either startswith "_"
    // or not. Between elements startswith "_", there is at most one element that does not startswith "_".
    // a scan can be performed to group those keys.
    std::vector<std::vector<std::string>> keys;
    std::vector<std::string> kcache;
    std::vector<std::string> values;
    // first drop all elements that does not startswith "_" before the first element that startswith "_"
    std::vector<std::string> block_ = block;
    auto it = std::find_if(block.begin(), block.end(), [](const std::string& s) { return FmtCore::startswith(s, "_"); });
    if (it != block.begin())
    {
        block_.erase(block_.begin(), it);
    }
    for (auto& elem: block_)
    {
        if (FmtCore::startswith(elem, "_"))
        {
            kcache.push_back(elem);
        }
        else
        {
            keys.push_back(kcache);
            values.push_back(elem);
            kcache.clear();
        }
    }
    assert(keys.size() == values.size()); // ensure the number of keys and values are the same
    // then for each elem in keys, if there are more than one element, then it is a table. Make it a table
    // , otherwise it is a simple key-value pair, directly add it to the output.
    std::map<std::string, std::vector<std::string>> out;
    for (size_t i = 0; i < keys.size(); i++)
    {
        if (keys[i].size() > 1)
        {
            const std::vector<std::string> words = _split_outside_enclose(values[i], " ", {"'", "'"});
            std::map<std::string, std::vector<std::string>> table = _build_table(keys[i], words);
            out.insert(table.begin(), table.end());
        }
        else
        {
            out[keys[i][0]] = {values[i]};
        }
    }
    
    return out;
}

void bcast_cifmap(std::map<std::string, std::vector<std::string>>& map, // the map to be broadcasted 
                  const int rank = 0)   // source rank: from which rank to broadcast
{
#ifdef __MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // use Parallel_Common namespace bcast_int and bcast_string to broadcast the size of map and key, value pairs
    int size = map.size();
    Parallel_Common::bcast_int(size); // seems it can only broadcast from rank 0, so presently no need to specify
                                      // the rank to broadcast
    std::vector<std::string> keys(size);
    std::vector<std::vector<std::string>> values(size);
    int i = 0;
    if (myrank == rank) // if the rank is the source rank, then pack the map to keys and values
    {
        for (auto& elem: map)
        {
            keys[i] = elem.first;
            values[i] = elem.second;
            i++;
        }
    }
    for (int i = 0; i < size; i++)
    {
        Parallel_Common::bcast_string(keys[i]);
        int valsize = values[i].size();
        Parallel_Common::bcast_int(valsize);
        values[i].resize(valsize);
        Parallel_Common::bcast_string(values[i].data(), valsize);
    }
    if (myrank != rank) // if the rank is not the source rank, then unpack the keys and values to map
    {
        map.clear();
        for (int i = 0; i < size; i++)
        {
            map[keys[i]] = values[i];
        }
    }
#endif
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
    assert(natom > 0); // ensure the number of atoms is positive
    atom_site_labels.resize(natom);
    atom_site_fract_coords.resize(3 * natom);
    for (int i = 0; i < natom; ++i)
    {
        atom_site_labels[i] = ucell.atoms[ucell.iat2it[i]].ncpp.psd; // the most standard label
        atom_site_labels[i] = atom_site_labels[i].empty() ? ucell.atom_label[ucell.iat2it[i]]: atom_site_labels[i];
        atom_site_labels[i] = atom_site_labels[i].empty() ? ucell.atoms[ucell.iat2it[i]].label: atom_site_labels[i];
        assert(!atom_site_labels[i].empty()); // ensure the label is not empty
        atom_site_fract_coords[3 * i] = ucell.atoms[ucell.iat2it[i]].taud[ucell.iat2ia[i]].x;
        atom_site_fract_coords[3 * i + 1] = ucell.atoms[ucell.iat2it[i]].taud[ucell.iat2ia[i]].y;
        atom_site_fract_coords[3 * i + 2] = ucell.atoms[ucell.iat2it[i]].taud[ucell.iat2ia[i]].z;
    }
}

void ModuleIO::CifParser::write(const std::string& fcif,
                                const double* abc_angles,
                                const int natom,
                                const std::string* atom_site_labels,
                                const double* atom_site_fract_coords,
                                const std::string& title,
                                const std::string& data_tag,
                                const int rank,
                                const double* atom_site_occups,
                                const std::string& cell_formula_units_z)
{
#ifdef __MPI // well...very simple...
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank != rank) // if present rank is not the rank assigned to write the cif file, then return
    {
        return;
    }
#endif
    std::ofstream ofs(fcif);
    if (!ofs)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::CifParser::write", "Cannot open file " + fcif);
    }
    ofs << title << std::endl;
    ofs << data_tag << std::endl;
    ofs << "_symmetry_space_group_name_H-M   'P 1'" << std::endl;
    ofs << "_cell_length_a " << FmtCore::format("%12.8f", abc_angles[0]) << std::endl;
    ofs << "_cell_length_b " << FmtCore::format("%12.8f", abc_angles[1]) << std::endl;
    ofs << "_cell_length_c " << FmtCore::format("%12.8f", abc_angles[2]) << std::endl;
    ofs << "_cell_angle_alpha " << FmtCore::format("%12.8f", abc_angles[3]) << std::endl;
    ofs << "_cell_angle_beta " << FmtCore::format("%12.8f", abc_angles[4]) << std::endl;
    ofs << "_cell_angle_gamma " << FmtCore::format("%12.8f", abc_angles[5]) << std::endl;
    ofs << "_symmetry_Int_Tables_number   1" << std::endl;
    std::string chem_sum, chem_structural;
    _build_chem_formula(natom, atom_site_labels, chem_sum, chem_structural);
    ofs << "_chemical_formula_structural   " << chem_structural << std::endl;
    ofs << "_chemical_formula_sum   " << chem_sum << std::endl;
    ofs << "_cell_volume " << FmtCore::format("%14.8f", abc_angles_to_volume(abc_angles)) << std::endl;
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
    FmtCore fmt("%3s %4s %3d %12.8f %12.8f %12.8f %2.1f");
    int j = 0;
    std::string numbered_label;
    std::string cache;
    for (int i = 0; i < natom; ++i)
    {
        const std::string label = atom_site_labels[i];
        numbered_label = label + std::to_string(i);
        cache = fmt.format(label, numbered_label, 1, 
        atom_site_fract_coords[j], atom_site_fract_coords[j + 1], atom_site_fract_coords[j + 2], occups[i]);
        ofs << cache << std::endl;
        j += 3;
    }
    ofs << std::endl;
    ofs.close();
}

void ModuleIO::CifParser::write(const std::string& fcif,
                                const std::vector<double>& abc_angles,
                                const std::vector<std::string>& atom_site_labels,
                                const std::vector<double>& atom_site_fract_coords,
                                const std::string& title,
                                const std::string& data_tag,
                                const int rank,
                                const std::vector<double>& atom_site_occups,
                                const std::string& cell_formula_units_z)
{
#ifdef __MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank != rank) // if present rank is not the rank assigned to write the cif file, then return
    {
        return;
    }
#endif
    const double* occups = atom_site_occups.empty() ? nullptr : atom_site_occups.data();
    write(fcif.c_str(), 
          abc_angles.data(), 
          atom_site_labels.size(), 
          atom_site_labels.data(), 
          atom_site_fract_coords.data(), 
          title, 
          data_tag, 
          rank,
          occups, 
          cell_formula_units_z);
}

void ModuleIO::CifParser::write(const std::string& fcif,
                                const UnitCell& ucell,
                                const std::string& title,
                                const std::string& data_tag,
                                const int rank)
{
#ifdef __MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank != rank)  // if present rank is not the rank assigned to write the cif file, then return
    {
        return;
    }
#endif
    std::vector<double> veca, vecb, vecc;
    int natom;
    std::vector<std::string> atom_site_labels;
    std::vector<double> atom_site_fract_coords;
    _unpack_ucell(ucell, veca, vecb, vecc, natom, atom_site_labels, atom_site_fract_coords);
    std::vector<double> vec(9);
    std::copy(veca.begin(), veca.end(), vec.begin());
    std::copy(vecb.begin(), vecb.end(), vec.begin() + 3);
    std::copy(vecc.begin(), vecc.end(), vec.begin() + 6);
    std::vector<double> abc_angles(6);
    vec_to_abc_angles(vec.data(), abc_angles.data());
    write(fcif.c_str(), abc_angles.data(), natom, atom_site_labels.data(), atom_site_fract_coords.data(), title, data_tag);
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

void ModuleIO::CifParser::read(const std::string& fcif,
                               std::map<std::string, std::vector<std::string>>& out,
                               const int rank)
{
    // okey for read, cannot just use if rank != 0 then return, because need to broadcast the map
    out.clear();
#ifdef __MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == rank) // only the rank assigned to read the cif file will read the file
    {
#endif
    std::ifstream ifs(fcif);
    if (!ifs)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::CifParser::read", "Cannot open file " + fcif);
    }
    std::string cache; // first read all lines into cache
    while (ifs.good())
    {
        std::string line;
        std::getline(ifs, line);
        if (FmtCore::startswith(FmtCore::strip(line), "#"))
        {
            out["comment"].push_back(line);
        }
        else if (FmtCore::startswith(FmtCore::strip(line), "data_"))
        {
            out["data_tag"].push_back(line);
        }
        else
        {
            cache += line + " ";
        }
    }
    std::vector<std::string> blocks = FmtCore::split(FmtCore::strip(cache), "loop_");
    blocks.erase(std::remove_if(blocks.begin(), blocks.end(), [](const std::string& s) { return s.empty(); }), blocks.end());
    
    for (auto& block: blocks)
    {
        std::vector<std::string> words = _split_loop_block(block);
        std::map<std::string, std::vector<std::string>> data = _build_block_data(words);
        out.insert(data.begin(), data.end());
    }
#ifdef __MPI
    }
    bcast_cifmap(out, rank);
#endif
}

ModuleIO::CifParser::CifParser(const std::string& fcif)
{
    read(fcif, raw_);
}

std::vector<std::string> ModuleIO::CifParser::get(const std::string& key)
{
    if (raw_.find(key) != raw_.end())
    {
        return raw_[key];
    }
    else
    {
        return std::vector<std::string>();
    }
}
