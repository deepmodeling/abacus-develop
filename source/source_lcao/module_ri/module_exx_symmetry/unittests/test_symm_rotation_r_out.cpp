// Regression for the symrot_R.txt writers: the spatial-only variant used by
// nspin=4, the antiunitary trailer, and the failure paths.
#include "source_cell/module_symmetry/symmetry.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::fprintf(stderr, "%s\n", message);
        std::exit(EXIT_FAILURE);
    }
}

std::string make_temp_directory()
{
    char pattern[] = "/tmp/abacus_symrot_writer_XXXXXX";
    const char* created = mkdtemp(pattern);
    require(created != nullptr, "could not create a temporary directory");
    return std::string(created) + "/";
}

std::vector<std::string> read_lines(const std::string& path)
{
    std::ifstream input(path);
    require(input.good(), "the writer produced no file");
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        if (!line.empty())
        {
            lines.push_back(line);
        }
    }
    return lines;
}

// One unitary operation (identity) and, when magnetic, one antiunitary one.
// Symmetry stores the operations in fixed-size arrays.
void fill_symmetry(ModuleSymmetry::Symmetry& symmetry, const bool magnetic)
{
    symmetry.nrotk = 1;
    symmetry.gmatrix[0] = ModuleBase::Matrix3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    symmetry.gtrans[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    symmetry.magnetic_nspin4 = magnetic;
    symmetry.nrotk_anti = magnetic ? 1 : 0;
    if (magnetic)
    {
        symmetry.gmatrix_anti[0] = ModuleBase::Matrix3(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0);
        symmetry.gtrans_anti[0] = ModuleBase::Vector3<double>(0.5, 0.5, 0.5);
    }
}
} // namespace

int main()
{
    const std::string directory = make_temp_directory();

    // Spatial-only export, non-magnetic: a marker line and one operation.
    {
        ModuleSymmetry::Symmetry symmetry;
        fill_symmetry(symmetry, false);
        std::ostringstream log;
        ModuleSymmetry::print_symrot_info_R_spatial(symmetry, 2, 3, directory, log);
        const auto lines = read_lines(directory + "symrot_R.txt");
        bool has_marker = false;
        for (const std::string& line : lines)
        {
            if (line.find("Spatial-only export") != std::string::npos)
            {
                has_marker = true;
            }
        }
        require(has_marker, "the spatial export did not mark itself");
        for (const std::string& line : lines)
        {
            require(line.rfind("Antiunitary operations:", 0) != 0,
                    "a non-magnetic export must not carry an antiunitary trailer");
        }
        require(log.str().find("symrot_R.txt") != std::string::npos,
                "the writer did not announce its path in the log");
    }

    // Magnetic: the antiunitary trailer carries the operation.
    {
        ModuleSymmetry::Symmetry symmetry;
        fill_symmetry(symmetry, true);
        std::ostringstream log;
        ModuleSymmetry::print_symrot_info_R_spatial(symmetry, 2, 3, directory, log);
        const auto lines = read_lines(directory + "symrot_R.txt");
        bool has_trailer = false;
        int trailer_index = -1;
        for (int index = 0; index < static_cast<int>(lines.size()); ++index)
        {
            if (lines[static_cast<std::size_t>(index)].rfind("Antiunitary operations:", 0) == 0)
            {
                has_trailer = true;
                trailer_index = index;
            }
        }
        require(has_trailer, "the magnetic export lost its antiunitary trailer");
        require(lines[static_cast<std::size_t>(trailer_index)].find(" 1") != std::string::npos,
                "the magnetic export did not record one antiunitary operation");
        require(static_cast<int>(lines.size()) > trailer_index + 4,
                "the antiunitary operation body is missing");
    }

    // An unwritable directory must throw a path-bearing error, not return.
    {
        ModuleSymmetry::Symmetry symmetry;
        fill_symmetry(symmetry, false);
        std::ostringstream log;
        bool threw = false;
        std::string message;
        try
        {
            ModuleSymmetry::print_symrot_info_R_spatial(symmetry, 2, 3,
                                                        "/nonexistent-directory-for-symrot-test/",
                                                        log);
        }
        catch (const std::exception& error)
        {
            threw = true;
            message = error.what();
        }
        require(threw, "an unwritable output directory did not throw");
        require(message.find("symrot_R.txt") != std::string::npos,
                "the error message does not name the file");
    }

    std::remove((directory + "symrot_R.txt").c_str());
    rmdir(directory.c_str());
    std::printf("symrot writer tests passed\n");
    return EXIT_SUCCESS;
}
