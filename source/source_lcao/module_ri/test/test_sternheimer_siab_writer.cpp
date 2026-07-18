#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/read_input.h"
#include "source_lcao/module_ri/sternheimer_siab_writer.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

namespace siab = module_ri::sternheimer_siab;
using Complex = std::complex<double>;

const char* canonical_fixture = R"fixture(<STERNHEIMER_SIAB_HEADER>
format_version 1
n_reference 2
n_primitive 4
n_blocks 2
grid_volume_bohr3 0.125
</STERNHEIMER_SIAB_HEADER>
<PRIMITIVE_BLOCKS>
# element atom_index l m n_primitive offset
H 0 0 0 2 0
H 0 1 0 2 2
</PRIMITIVE_BLOCKS>
<REFERENCE_METADATA>
# occupied_state auxiliary_channel frequency_ha occupation frequency_weight norm
0 0 0.5 2.0 0.3 1.2
0 0 1.5 2.0 0.7 0.8
</REFERENCE_METADATA>
<OVERLAP_Q>
# row-major <Y_rho|B_e>, real imag
0.8 0.0
0.1 0.0
0.2 0.0
0.0 0.0
0.4 0.0
0.3 0.0
0.0 0.0
0.1 0.0
</OVERLAP_Q>
<OVERLAP_S>
# row-major <B_e|B_ep>, real imag
1.0 0.0
0.0 0.0
0.0 0.0
0.0 0.0
0.0 0.0
1.0 0.0
0.0 0.0
0.0 0.0
0.0 0.0
0.0 0.0
1.0 0.0
0.0 0.0
0.0 0.0
0.0 0.0
0.0 0.0
1.0 0.0
</OVERLAP_S>
<PROVENANCE_JSON>
{"abacus_commit":"19ab21e01d02cc805604ed77a6e269af698fdd1d","auxiliary_basis_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","cell_bohr":[20.0,0.0,0.0,0.0,20.0,0.0,0.0,0.0,20.0],"ecut_ry":25.0,"kernel":"full_coulomb","orbital_sha256":"7e398340398306a6baf1c61ea68944d81ed43667473fbcc290d6541c4a661d1c","pseudopotential_sha256":"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb","spin_convention":"occupation_in_metadata"}
</PROVENANCE_JSON>
)fixture";

std::string test_path(const std::string& suffix)
{
    return ::testing::TempDir() + "/sternheimer_siab_writer_" + suffix + ".dat";
}

std::string read_text(const std::string& path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

void write_text(const std::string& path, const std::string& text)
{
    std::ofstream output(path.c_str(), std::ios::binary | std::ios::trunc);
    output << text;
}

std::size_t count_occurrences(const std::string& text, const std::string& needle)
{
    std::size_t count = 0;
    std::size_t position = 0;
    while ((position = text.find(needle, position)) != std::string::npos)
    {
        ++count;
        position += needle.size();
    }
    return count;
}

std::string section_body(const std::string& text, const std::string& section)
{
    const std::string opening = "<" + section + ">\n";
    const std::string closing = "</" + section + ">";
    const std::size_t begin = text.find(opening);
    const std::size_t end = text.find(closing);
    EXPECT_NE(begin, std::string::npos);
    EXPECT_NE(end, std::string::npos);
    if (begin == std::string::npos || end == std::string::npos)
    {
        return "";
    }
    return text.substr(begin + opening.size(), end - begin - opening.size());
}

std::size_t data_line_count(const std::string& body)
{
    std::istringstream input(body);
    std::string line;
    std::size_t count = 0;
    while (std::getline(input, line))
    {
        if (!line.empty() && line[0] != '#')
        {
            ++count;
        }
    }
    return count;
}

std::vector<siab::PrimitiveBlock> canonical_blocks()
{
    return {siab::PrimitiveBlock{"H", 0, 0, 0, 2, 0}, siab::PrimitiveBlock{"H", 0, 1, 0, 2, 2}};
}

std::vector<siab::ReferenceRow> canonical_rows_reversed()
{
    siab::ReferenceRow high;
    high.occupied_state = 0;
    high.auxiliary_channel = 0;
    high.frequency_ha = 1.5;
    high.occupation = 2.0;
    high.frequency_weight = 0.7;
    high.norm = 0.8;
    high.q = {{0.4, 0.0}, {0.3, 0.0}, {0.0, 0.0}, {0.1, 0.0}};
    high.frequency_index = 1;

    siab::ReferenceRow low;
    low.occupied_state = 0;
    low.auxiliary_channel = 0;
    low.frequency_ha = 0.5;
    low.occupation = 2.0;
    low.frequency_weight = 0.3;
    low.norm = 1.2;
    low.q = {{0.8, 0.0}, {0.1, 0.0}, {0.2, 0.0}, {0.0, 0.0}};
    low.frequency_index = 0;
    return {high, low};
}

std::vector<Complex> canonical_overlap_s()
{
    std::vector<Complex> overlap_s(16, Complex(0.0, 0.0));
    for (std::size_t i = 0; i != 4; ++i)
    {
        overlap_s[i * 4 + i] = Complex(1.0, 0.0);
    }
    return overlap_s;
}

siab::Provenance canonical_provenance()
{
    siab::Provenance provenance;
    provenance.abacus_commit = "19ab21e01d02cc805604ed77a6e269af698fdd1d";
    provenance.auxiliary_basis_sha256 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    provenance.cell_bohr = {20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0};
    provenance.ecut_ry = 25.0;
    provenance.kernel = "full_coulomb";
    provenance.orbital_sha256 = "7e398340398306a6baf1c61ea68944d81ed43667473fbcc290d6541c4a661d1c";
    provenance.pseudopotential_sha256 = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
    provenance.spin_convention = "occupation_in_metadata";
    return provenance;
}

void write_canonical(const std::string& path,
                     const std::vector<siab::ReferenceRow>& rows = canonical_rows_reversed(),
                     const std::vector<Complex>& overlap_s = canonical_overlap_s(),
                     const siab::Provenance& provenance = canonical_provenance())
{
    siab::write_v1(path, 0.125, canonical_blocks(), rows, overlap_s, provenance);
}

std::string formatted_q_real(const double value, const std::string& suffix)
{
    const std::string path = test_path("float_" + suffix);
    std::vector<siab::ReferenceRow> rows = canonical_rows_reversed();
    rows[1].q[0] = Complex(value, 0.0);
    write_canonical(path, rows);

    std::istringstream input(section_body(read_text(path), "OVERLAP_Q"));
    std::string line;
    std::string real;
    while (std::getline(input, line))
    {
        if (!line.empty() && line[0] != '#')
        {
            std::istringstream values(line);
            values >> real;
            break;
        }
    }
    std::remove(path.c_str());
    return real;
}

double parse_classic_double(const std::string& text)
{
    std::istringstream input(text);
    input.imbue(std::locale::classic());
    double value = 0.0;
    input >> value;
    const bool exact_subnormal = value != 0.0 && std::abs(value) < std::numeric_limits<double>::min();
    EXPECT_TRUE(!input.fail() || exact_subnormal);
    EXPECT_EQ(input.rdbuf()->in_avail(), 0);
    return value;
}

} // namespace

TEST(SternheimerSIABWriter, MatchesCanonicalFixtureAndSortsRowsWithQ)
{
    const std::string path = test_path("canonical");
    std::remove(path.c_str());
    std::remove((path + ".tmp").c_str());

    write_canonical(path);
    const std::string text = read_text(path);
    EXPECT_EQ(text, canonical_fixture);

    const std::vector<std::string> sections = {"STERNHEIMER_SIAB_HEADER",
                                               "PRIMITIVE_BLOCKS",
                                               "REFERENCE_METADATA",
                                               "OVERLAP_Q",
                                               "OVERLAP_S",
                                               "PROVENANCE_JSON"};
    std::size_t previous = 0;
    for (const std::string& section: sections)
    {
        const std::string opening = "<" + section + ">";
        const std::string closing = "</" + section + ">";
        EXPECT_EQ(count_occurrences(text, opening), 1);
        EXPECT_EQ(count_occurrences(text, closing), 1);
        const std::size_t position = text.find(opening);
        EXPECT_NE(position, std::string::npos);
        EXPECT_GE(position, previous);
        previous = position;
    }

    EXPECT_EQ(data_line_count(section_body(text, "PRIMITIVE_BLOCKS")), 2);
    EXPECT_EQ(data_line_count(section_body(text, "REFERENCE_METADATA")), 2);
    EXPECT_EQ(data_line_count(section_body(text, "OVERLAP_Q")), 8);
    EXPECT_EQ(data_line_count(section_body(text, "OVERLAP_S")), 16);

    const std::string expected_json
        = "{\"abacus_commit\":\"19ab21e01d02cc805604ed77a6e269af698fdd1d\","
          "\"auxiliary_basis_sha256\":\"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\","
          "\"cell_bohr\":[20.0,0.0,0.0,0.0,20.0,0.0,0.0,0.0,20.0],\"ecut_ry\":25.0,"
          "\"kernel\":\"full_coulomb\","
          "\"orbital_sha256\":\"7e398340398306a6baf1c61ea68944d81ed43667473fbcc290d6541c4a661d1c\","
          "\"pseudopotential_sha256\":\"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\","
          "\"spin_convention\":\"occupation_in_metadata\"}\n";
    EXPECT_EQ(section_body(text, "PROVENANCE_JSON"), expected_json);

    std::remove(path.c_str());
}

TEST(SternheimerSIABWriter, PreservesShortestExactDoubleRoundTrip)
{
    const double adjacent = std::nextafter(1.0, 2.0);
    const std::vector<std::pair<double, std::string>> values = {
        {adjacent, "adjacent"},
        {std::numeric_limits<double>::denorm_min(), "denorm_min"},
        {std::numeric_limits<double>::min(), "min"},
        {std::numeric_limits<double>::max(), "max"},
        {-0.0, "negative_zero"},
    };
    for (const std::pair<double, std::string>& value: values)
    {
        const std::string formatted = formatted_q_real(value.first, value.second);
        const double parsed = parse_classic_double(formatted);
        EXPECT_EQ(parsed, value.first);
        if (value.first == 0.0)
        {
            EXPECT_EQ(std::signbit(parsed), std::signbit(value.first));
        }
    }

    EXPECT_EQ(formatted_q_real(adjacent, "adjacent_exact"), "1.0000000000000002");
    EXPECT_EQ(formatted_q_real(123.0, "one_hundred_twenty_three"), "123.0");
    EXPECT_EQ(formatted_q_real(10000.0, "ten_thousand"), "1e4");
    EXPECT_EQ(formatted_q_real(1.0e-5, "one_e_minus_five"), "1e-5");
    EXPECT_EQ(formatted_q_real(-0.0, "negative_zero_exact"), "-0.0");
}

TEST(SternheimerSIABWriter, EscapesJsonStringsWithoutChangingUtf8Bytes)
{
    const std::string path = test_path("json_escape");
    siab::Provenance provenance = canonical_provenance();
    const std::string utf8_two = "\xC2\xA2";
    const std::string utf8_three = "\xE6\xA0\xB8";
    const std::string utf8_four = "\xF0\x9F\x98\x80";
    const std::string utf8 = utf8_two + "_" + utf8_three + "_" + utf8_four;
    provenance.kernel = std::string("full_\"\\\b\f\n\r\t") + static_cast<char>(0x01) + "_" + utf8;

    write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance);
    const std::string json = section_body(read_text(path), "PROVENANCE_JSON");
    const std::string expected_kernel = std::string(R"json("kernel":"full_\"\\\b\f\n\r\t\u0001_)json") + utf8 + "\"";
    EXPECT_NE(json.find(expected_kernel), std::string::npos);
    EXPECT_EQ(json.find(static_cast<char>(0x01)), std::string::npos);

    std::remove(path.c_str());
}

TEST(SternheimerSIABWriter, RejectsMalformedUtf8BeforeTouchingTmp)
{
    const std::string path = test_path("invalid_utf8");
    const std::string tmp = path + ".tmp";
    const std::vector<std::string> malformed = {
        "\x80",             // lone continuation byte
        "\xE2\x28\xA1",     // non-continuation byte inside a sequence
        "\xE2\x82",         // truncated three-byte sequence
        "\xC0\xAF",         // overlong slash
        "\xED\xA0\x80",     // UTF-8 encoded surrogate U+D800
        "\xF4\x90\x80\x80", // code point above U+10FFFF
    };

    write_text(tmp, "sentinel");
    for (const std::string& invalid: malformed)
    {
        siab::Provenance provenance = canonical_provenance();
        provenance.kernel = invalid;
        EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                     std::invalid_argument);
        EXPECT_EQ(read_text(tmp), "sentinel");
    }

    using StringMember = std::string siab::Provenance::*;
    const std::vector<StringMember> json_string_fields = {
        &siab::Provenance::abacus_commit,
        &siab::Provenance::auxiliary_basis_sha256,
        &siab::Provenance::kernel,
        &siab::Provenance::orbital_sha256,
        &siab::Provenance::pseudopotential_sha256,
        &siab::Provenance::spin_convention,
    };
    for (const StringMember field: json_string_fields)
    {
        siab::Provenance provenance = canonical_provenance();
        provenance.*field = "\x80";
        EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                     std::invalid_argument);
        EXPECT_EQ(read_text(tmp), "sentinel");
    }

    std::remove(tmp.c_str());
}

TEST(SternheimerSIABWriter, AcceptsCompleteNonorthogonalCellWithNegativeComponents)
{
    const std::string path = test_path("nonorthogonal_cell");
    siab::Provenance provenance = canonical_provenance();
    provenance.cell_bohr = {20.0, 1.0, -2.0, 0.0, 19.0, 3.0, -1.0, 0.0, 18.0};

    write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance);
    const std::string json = section_body(read_text(path), "PROVENANCE_JSON");
    EXPECT_NE(json.find("\"cell_bohr\":[20.0,1.0,-2.0,0.0,19.0,3.0,-1.0,0.0,18.0]"), std::string::npos);

    std::remove(path.c_str());
}

TEST(SternheimerSIABWriter, RejectsIncompleteSingularAndNonfiniteCells)
{
    const std::string path = test_path("invalid_cell");
    siab::Provenance provenance = canonical_provenance();
    provenance.cell_bohr = {20.0, 20.0, 20.0};
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                 std::invalid_argument);

    provenance = canonical_provenance();
    provenance.cell_bohr = {1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                 std::invalid_argument);

    provenance = canonical_provenance();
    provenance.cell_bohr[4] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                 std::invalid_argument);
}

TEST(SternheimerSIABWriter, RejectsBadDimensionsAndDuplicateRowKeys)
{
    const std::string path = test_path("dimensions");
    std::remove(path.c_str());
    std::remove((path + ".tmp").c_str());
    std::vector<siab::ReferenceRow> rows = canonical_rows_reversed();
    rows[0].q.pop_back();
    EXPECT_THROW(write_canonical(path, rows), std::invalid_argument);

    rows = canonical_rows_reversed();
    rows[0].frequency_index = rows[1].frequency_index;
    EXPECT_THROW(write_canonical(path, rows), std::invalid_argument);

    std::vector<Complex> overlap_s = canonical_overlap_s();
    overlap_s.pop_back();
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), overlap_s), std::invalid_argument);
    EXPECT_FALSE(std::ifstream(path.c_str()).good());
}

TEST(SternheimerSIABWriter, RejectsNonHermitianAndNonfiniteValues)
{
    const std::string path = test_path("values");
    std::vector<Complex> overlap_s = canonical_overlap_s();
    overlap_s[1] = Complex(0.25, 0.5);
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), overlap_s), std::invalid_argument);

    std::vector<siab::ReferenceRow> rows = canonical_rows_reversed();
    rows[0].frequency_weight = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(write_canonical(path, rows), std::invalid_argument);

    rows = canonical_rows_reversed();
    rows[0].q[0] = Complex(std::numeric_limits<double>::infinity(), 0.0);
    EXPECT_THROW(write_canonical(path, rows), std::invalid_argument);

    overlap_s = canonical_overlap_s();
    overlap_s[0] = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), overlap_s), std::invalid_argument);
}

TEST(SternheimerSIABWriter, AcceptsAsciiPrimitiveElementTokens)
{
    const std::string path = test_path("element_tokens");
    const std::vector<std::string> accepted = {"H", "H1", "Si_aux", "H+", "H-", "A+B-C_1"};
    for (const std::string& element: accepted)
    {
        std::vector<siab::PrimitiveBlock> blocks = canonical_blocks();
        blocks[0].element = element;
        siab::write_v1(path, 0.125, blocks, canonical_rows_reversed(), canonical_overlap_s(), canonical_provenance());
        EXPECT_NE(read_text(path).find(element + " 0 0 0 2 0\n"), std::string::npos);
    }
    std::remove(path.c_str());
}

TEST(SternheimerSIABWriter, RejectsUnsafePrimitiveElementTokensBeforeTouchingTmp)
{
    const std::string path = test_path("invalid_element_tokens");
    const std::string tmp = path + ".tmp";
    const std::vector<std::string> rejected = {"", " H", "H H", "H\nX", "#H", "<H>", ">H", "+H", "-H", "\x80", "H\x80"};
    write_text(tmp, "sentinel");
    for (const std::string& element: rejected)
    {
        std::vector<siab::PrimitiveBlock> blocks = canonical_blocks();
        blocks[0].element = element;
        EXPECT_THROW(siab::write_v1(path,
                                    0.125,
                                    blocks,
                                    canonical_rows_reversed(),
                                    canonical_overlap_s(),
                                    canonical_provenance()),
                     std::invalid_argument);
        EXPECT_EQ(read_text(tmp), "sentinel");
    }
    std::remove(tmp.c_str());
}

TEST(SternheimerSIABWriter, RejectsInvalidBlocksAndHashesBeforeTouchingTmp)
{
    const std::string path = test_path("validation");
    const std::string tmp = path + ".tmp";
    write_text(tmp, "sentinel");

    std::vector<siab::PrimitiveBlock> blocks = canonical_blocks();
    blocks[1].offset = 3;
    EXPECT_THROW(
        siab::write_v1(path, 0.125, blocks, canonical_rows_reversed(), canonical_overlap_s(), canonical_provenance()),
        std::invalid_argument);
    EXPECT_EQ(read_text(tmp), "sentinel");

    siab::Provenance provenance = canonical_provenance();
    provenance.orbital_sha256 = "not-a-sha256";
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                 std::invalid_argument);
    EXPECT_EQ(read_text(tmp), "sentinel");

    std::remove(tmp.c_str());
}

TEST(SternheimerSIABWriter, AtomicallyReplacesExistingPathAndRemovesTmp)
{
    const std::string path = test_path("atomic");
    const std::string tmp = path + ".tmp";
    write_text(path, "old contents");
    std::remove(tmp.c_str());

    write_canonical(path);
    EXPECT_EQ(read_text(path), canonical_fixture);
    EXPECT_FALSE(std::ifstream(tmp.c_str()).good());

    std::remove(path.c_str());
}

TEST(SternheimerSIABWriter, AppendsCompleteTask4ProvenanceDeterministically)
{
    const std::string path = test_path("task4_provenance");
    siab::Provenance provenance = canonical_provenance();
    provenance.executable_sha256 = "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc";
    provenance.exx_pca_thr = 1.0e-6;
    provenance.sternheimer_nfreq = 2;
    provenance.frequency_ha = {0.2, 0.8};
    provenance.frequency_weights_ha = {0.3, 0.7};
    provenance.mpi_ranks = 2;
    provenance.omp_threads = 32;

    write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance);
    const std::string json = section_body(read_text(path), "PROVENANCE_JSON");
    EXPECT_NE(json.find(
                  "\"spin_convention\":\"occupation_in_metadata\","
                  "\"executable_sha256\":\"cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\","
                  "\"exx_pca_thr\":1e-6,\"frequency_ha\":[0.2,0.8],"
                  "\"frequency_weights_ha\":[0.3,0.7],\"mpi_ranks\":2,\"omp_threads\":32,"
                  "\"sternheimer_nfreq\":2}"),
              std::string::npos);
    std::remove(path.c_str());
}

TEST(SternheimerSIABWriter, RejectsPartialTask4ProvenanceBeforeTouchingTmp)
{
    const std::string path = test_path("partial_task4_provenance");
    const std::string tmp = path + ".tmp";
    write_text(tmp, "sentinel");

    siab::Provenance provenance = canonical_provenance();
    provenance.executable_sha256 = "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc";
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                 std::invalid_argument);
    EXPECT_EQ(read_text(tmp), "sentinel");

    provenance.exx_pca_thr = 1.0e-6;
    provenance.sternheimer_nfreq = 1;
    provenance.frequency_ha = {0.2};
    provenance.frequency_weights_ha = {};
    provenance.mpi_ranks = 1;
    provenance.omp_threads = 1;
    EXPECT_THROW(write_canonical(path, canonical_rows_reversed(), canonical_overlap_s(), provenance),
                 std::invalid_argument);
    EXPECT_EQ(read_text(tmp), "sentinel");
    std::remove(tmp.c_str());
}

TEST(SternheimerSIABInput, RegisteredCheckEnforcesLcaoDeltaCombination)
{
    ModuleIO::ReadInput read_input(0);
    const std::vector<std::pair<std::string, ModuleIO::Input_Item>>& items = read_input.get_input_lists();
    const auto found = std::find_if(items.begin(), items.end(), [](const auto& item) {
        return item.first == "out_sternheimer_siab";
    });
    ASSERT_NE(found, items.end());
    const ModuleIO::Input_Item& item = found->second;
    ASSERT_TRUE(static_cast<bool>(item.check_value));
    EXPECT_NE(item.description.find("basis_type=lcao"), std::string::npos);
    EXPECT_NE(item.availability.find("basis_type=lcao"), std::string::npos);

    Parameter parameter;
    Input_para& input = const_cast<Input_para&>(parameter.inp);
    EXPECT_FALSE(input.out_sternheimer_siab);
    item.check_value(item, parameter);

    input.out_sternheimer_siab = true;
    input.basis_type = "lcao";
    input.out_sternheimer_librpa = true;
    input.sternheimer_delta = true;
    input.bessel_nao_rcuts = {8.0};
    item.check_value(item, parameter);

    input.bessel_nao_rcuts.clear();
    EXPECT_EXIT(item.check_value(item, parameter), ::testing::ExitedWithCode(1), "");
    input.bessel_nao_rcuts = {8.0, 10.0};
    EXPECT_EXIT(item.check_value(item, parameter), ::testing::ExitedWithCode(1), "");
    input.bessel_nao_rcuts = {8.0};

    input.basis_type = "pw";
    EXPECT_EXIT(item.check_value(item, parameter), ::testing::ExitedWithCode(1), "");

    input.basis_type = "lcao";
    input.out_sternheimer_librpa = false;
    EXPECT_EXIT(item.check_value(item, parameter), ::testing::ExitedWithCode(1), "");

    input.out_sternheimer_librpa = true;
    input.sternheimer_delta = false;
    EXPECT_EXIT(item.check_value(item, parameter), ::testing::ExitedWithCode(1), "");
}
