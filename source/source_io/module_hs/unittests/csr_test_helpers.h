#ifndef CSR_TEST_HELPERS_H
#define CSR_TEST_HELPERS_H

/**
 * @file csr_test_helpers.h
 * @brief Shared helper utilities for the CSR writer unit tests in this
 *        directory (hsr_writer / hs_sparse_io / dhs_sparse_writer).
 *
 * Each test executable includes this header from exactly one translation
 * unit, so the anonymous-namespace helpers below are instantiated once per
 * test binary.
 */

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_hcontainer/atom_pair.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
std::string read_file(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

std::vector<int> read_binary_ints(const std::string& filename, const size_t count)
{
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    std::vector<int> values(count, 0);
    for (size_t i = 0; i < count; ++i)
    {
        ifs.read(reinterpret_cast<char*>(&values[i]), sizeof(int));
    }
    return values;
}

template <typename T>
T read_binary_value(std::ifstream& ifs)
{
    T value{};
    ifs.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}

struct NativeDoubleRBlock
{
    int rx = 0;
    int ry = 0;
    int rz = 0;
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<long long> row_ptr;
};

struct NativeDoubleRecord
{
    int step = 0;
    int nbasis = 0;
    std::vector<NativeDoubleRBlock> blocks;
};

NativeDoubleRecord read_native_double_record(std::ifstream& ifs)
{
    NativeDoubleRecord record;
    record.step = read_binary_value<int>(ifs);
    record.nbasis = read_binary_value<int>(ifs);
    const int nR = read_binary_value<int>(ifs);
    record.blocks.resize(nR);
    for (NativeDoubleRBlock& block: record.blocks)
    {
        block.rx = read_binary_value<int>(ifs);
        block.ry = read_binary_value<int>(ifs);
        block.rz = read_binary_value<int>(ifs);
        const int nnz = read_binary_value<int>(ifs);
        block.values.resize(nnz);
        block.columns.resize(nnz);
        block.row_ptr.resize(record.nbasis + 1);
        for (double& value: block.values)
        {
            value = read_binary_value<double>(ifs);
        }
        for (int& column: block.columns)
        {
            column = read_binary_value<int>(ifs);
        }
        for (long long& pointer: block.row_ptr)
        {
            pointer = read_binary_value<long long>(ifs);
        }
    }
    return record;
}

std::vector<std::string> read_lines(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    return lines;
}

int count_substr(const std::string& text, const std::string& pattern)
{
    int count = 0;
    std::string::size_type pos = 0;
    while ((pos = text.find(pattern, pos)) != std::string::npos)
    {
        ++count;
        pos += pattern.size();
    }
    return count;
}

void init_unitcell(UnitCell& ucell)
{
    ucell.latName = "user_defined_lattice";
    ucell.lat0 = 10.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e33 = 1.0;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Si";
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nw = 2;
    ucell.atoms[0].taud.resize(1);
    ucell.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.25, 0.5);
}

void init_serial_orbitals(Parallel_Orbitals& pv)
{
    pv.atom_begin_row.resize(2);
    pv.atom_begin_col.resize(2);
    pv.atom_begin_row[0] = 0;
    pv.atom_begin_row[1] = 2;
    pv.atom_begin_col[0] = 0;
    pv.atom_begin_col[1] = 2;
    pv.nrow = 2;
    pv.ncol = 2;
    pv.set_serial(2, 2);
}

void fill_matrix(hamilt::HContainer<double>& matrix, Parallel_Orbitals& pv, double* values)
{
    hamilt::AtomPair<double> pair(0, 0, 0, 0, 0, &pv, values);
    matrix.insert_pair(pair);
}

template <typename T>
void fill_matrix_at_R(hamilt::HContainer<T>& matrix,
                      Parallel_Orbitals& pv,
                      const int rx,
                      const int ry,
                      const int rz,
                      T* values)
{
    hamilt::AtomPair<T> pair(0, 0, rx, ry, rz, &pv, values);
    matrix.insert_pair(pair);
}

void init_sparse_output_globals()
{
    GlobalV::DRANK = 0;
}

void remove_derivative_files(const std::string& fileflag, int step = -1)
{
    std::vector<std::string> filenames;
    if (step >= 0)
    {
        const std::string suffix = "g" + std::to_string(step + 1) + "_nao.csr";
        for (const char axis : {'x', 'y', 'z'})
        {
            for (int ispin = 1; ispin <= 2; ++ispin)
            {
                filenames.push_back("d" + fileflag + "r" + axis + "s" + std::to_string(ispin) + suffix);
            }
        }
    }
    else
    {
        for (const char axis : {'x', 'y', 'z'})
        {
            for (int ispin = 1; ispin <= 2; ++ispin)
            {
                filenames.push_back("d" + fileflag + "r" + axis + "s" + std::to_string(ispin) + "_nao.csr");
            }
        }
    }
    for (const std::string& filename: filenames)
    {
        std::remove(filename.c_str());
    }
}

bool starts_with(const std::string& text, const std::string& prefix)
{
    return text.find(prefix) == 0;
}
} // namespace

#endif // CSR_TEST_HELPERS_H
