#include "../librpa_stru_symmetry.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

std::vector<std::string> split_lines(const std::string& text)
{
    std::vector<std::string> lines;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line))
    {
        lines.push_back(line);
    }
    return lines;
}

void require_operation_row(const std::string& line, const char* message)
{
    std::istringstream input(line);
    std::vector<std::string> fields;
    std::string field;
    while (input >> field)
    {
        fields.push_back(field);
    }
    require(fields.size() == 12, message);
    for (int i = 0; i < 9; ++i)
    {
        std::size_t consumed = 0;
        (void)std::stoi(fields[i], &consumed);
        require(consumed == fields[i].size(), "rotation fields must be integers");
    }
    for (int i = 9; i < 12; ++i)
    {
        std::size_t consumed = 0;
        const double value = std::stod(fields[i], &consumed);
        require(consumed == fields[i].size() && std::isfinite(value),
                "translation fields must be finite floating-point values");
    }
}

void test_common_spatial_and_magnetic_rows()
{
    const RpaLriDetail::LibRpaSymmetryOperation unitary[] = {
        {{1, 0, 0, 0, 1, 0, 0, 0, 1}, {0.0, 0.0, 0.0}},
        {{-1, 0, 0, 0, -1, 0, 0, 0, -1}, {0.5, 0.0, 0.0}}};
    const RpaLriDetail::LibRpaSymmetryOperation antiunitary[] = {
        {{1, 0, 0, 0, -1, 0, 0, 0, -1}, {0.125, 0.25, 0.375}}};
    const std::vector<RpaLriDetail::LibRpaSymmetryOperation> unitary_operations(unitary, unitary + 2);
    const std::vector<RpaLriDetail::LibRpaSymmetryOperation> antiunitary_operations(antiunitary, antiunitary + 1);

    std::ostringstream output;
    RpaLriDetail::write_librpa_symmetry_rows(output, unitary_operations, antiunitary_operations);
    const std::string serialized = output.str();
    const auto lines = split_lines(serialized);

    require(lines.size() == 4, "symmetry output must contain one count line and three rows");
    require(lines[0] == "3 row", "ordinary and antiunitary operations share one row block");
    require_operation_row(lines[1], "identity operation row must contain 12 fields");
    require_operation_row(lines[2], "ordinary inversion row must contain 12 fields");
    require_operation_row(lines[3], "antiunitary spatial row must contain 12 fields");
    require(lines[1].find("1.250000000000000e-01") == std::string::npos,
            "ordinary rows must precede the antiunitary row");
    require(lines[3].find("1.250000000000000e-01") != std::string::npos,
            "antiunitary translation must be preserved");
    require(serialized.find("spin_symmetry") == std::string::npos,
            "legacy spin-specific trailer must not be emitted");
}
} // namespace

int main()
{
    test_common_spatial_and_magnetic_rows();
    std::cout << "LibRPA stru_out symmetry tests passed\n";
    return 0;
}
