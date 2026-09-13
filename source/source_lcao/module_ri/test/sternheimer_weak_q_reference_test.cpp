#include <gtest/gtest.h>

#include "source_lcao/module_ri/sternheimer_weak_q_reference.h"
#include "source_lcao/module_ri/sternheimer_siab_provenance.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace
{
using namespace ModuleRI;
using Bytes = std::vector<unsigned char>;

class SternheimerWeakQReferenceTest : public ::testing::Test
{
  protected:
    std::string directory;
    std::string path;
    const std::string contract = std::string(64, 'a');
    SternheimerWeakQReferenceDimensions dimensions{3, 2, 2};
    std::vector<double> potential{1.25, -0.0, -2.5};
    std::vector<SternheimerLCAOOccupiedKPoint> records;

    void SetUp() override
    {
        char pattern[] = "/tmp/st-weak-q-reference-XXXXXX";
        const char* made = ::mkdtemp(pattern);
        ASSERT_NE(made, nullptr);
        directory = made;
        path = directory + "/reference.bin";
        SternheimerLCAOOccupiedKPoint r;
        r.local_k_index = -1;
        r.global_k_index = 4;
        r.zero_order_k_index = 2;
        r.symmetry_spatial_isym = 7;
        r.symmetry_time_reversal = true;
        r.spin_index = 1;
        r.kpoint = {0.25, -0.5, -0.0};
        r.has_grid_kpoint_override = true;
        r.grid_kpoint = {-0.75, 0.5, 1.0};
        r.kweight = 0.125;
        r.eigenvalues = {-1.25, -0.5};
        r.occupations = {2.0, 0.25};
        r.coefficients = {{{0.25, -0.0}, {1.5, -2}}, {{-3, 4}, {5, 6}}};
        r.unoccupied_eigenvalues = {0.25};
        r.unoccupied_coefficients = {{{7, 8}, {-9, -10}}};
        records.push_back(r);
        r.local_k_index = 3;
        r.global_k_index = 5;
        r.symmetry_time_reversal = false;
        r.has_grid_kpoint_override = false;
        r.unoccupied_eigenvalues.clear();
        r.unoccupied_coefficients.clear();
        records.push_back(r);
    }

    void TearDown() override
    {
        std::remove(path.c_str());
        std::remove((directory + "/copy.bin").c_str());
        if (!directory.empty()) EXPECT_EQ(::rmdir(directory.c_str()), 0);
    }

    void write()
    {
        write_sternheimer_weak_q_reference(path, contract, dimensions, potential, records);
    }

    SternheimerWeakQReference read()
    {
        return read_sternheimer_weak_q_reference(path, contract, dimensions);
    }

    Bytes bytes() const
    {
        std::ifstream in(path, std::ios::binary);
        return Bytes(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }

    void replace(const Bytes& data)
    {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out.write(reinterpret_cast<const char*>(data.data()), data.size());
        ASSERT_TRUE(out.good());
    }

    static void put64(Bytes& data, const std::size_t offset, const std::uint64_t value)
    {
        for (unsigned i = 0; i < 8; ++i) data.at(offset + i) = (value >> (8 * i)) & 255;
    }

    static void reseal(Bytes& data)
    {
        ASSERT_GE(data.size(), 64U);
        module_ri::sternheimer_siab::Sha256 sha;
        sha.update(data.data(), data.size() - 64);
        const auto digest = sha.finish();
        std::copy(digest.begin(), digest.end(), data.end() - 64);
    }

    std::size_t dimension_offset() const { return 8 + 4 + 8 + contract.size(); }
    std::size_t potential_offset() const { return dimension_offset() + 24; }
    std::size_t record_offset() const { return potential_offset() + 8 * potential.size(); }
};

void expect_bits(const double actual, const double expected)
{
    EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(double)), 0);
}

void expect_reals(const std::vector<double>& actual, const std::vector<double>& expected)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) expect_bits(actual[i], expected[i]);
}

void expect_coefficients(const std::vector<std::vector<std::complex<double>>>& actual,
                         const std::vector<std::vector<std::complex<double>>>& expected)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        ASSERT_EQ(actual[i].size(), expected[i].size());
        for (std::size_t j = 0; j < actual[i].size(); ++j)
        {
            expect_bits(actual[i][j].real(), expected[i][j].real());
            expect_bits(actual[i][j].imag(), expected[i][j].imag());
        }
    }
}

TEST_F(SternheimerWeakQReferenceTest, RoundTripsEveryActualFieldBitExactly)
{
    write();
    const auto restored = read();
    expect_reals(restored.fine_potential, potential);
    ASSERT_EQ(restored.full_records.size(), records.size());
    for (std::size_t i = 0; i < records.size(); ++i)
    {
        const auto& a = restored.full_records[i];
        const auto& b = records[i];
        EXPECT_EQ(a.local_k_index, b.local_k_index);
        EXPECT_EQ(a.global_k_index, b.global_k_index);
        EXPECT_EQ(a.zero_order_k_index, b.zero_order_k_index);
        EXPECT_EQ(a.symmetry_spatial_isym, b.symmetry_spatial_isym);
        EXPECT_EQ(a.symmetry_time_reversal, b.symmetry_time_reversal);
        EXPECT_EQ(a.spin_index, b.spin_index);
        EXPECT_EQ(a.has_grid_kpoint_override, b.has_grid_kpoint_override);
        for (std::size_t j = 0; j < 3; ++j)
        {
            expect_bits(a.kpoint[j], b.kpoint[j]);
            expect_bits(a.grid_kpoint[j], b.grid_kpoint[j]);
        }
        expect_bits(a.kweight, b.kweight);
        expect_reals(a.eigenvalues, b.eigenvalues);
        expect_reals(a.occupations, b.occupations);
        expect_coefficients(a.coefficients, b.coefficients);
        expect_reals(a.unoccupied_eigenvalues, b.unoccupied_eigenvalues);
        expect_coefficients(a.unoccupied_coefficients, b.unoccupied_coefficients);
    }
    const auto original = bytes();
    write_sternheimer_weak_q_reference(directory + "/copy.bin", contract, dimensions,
                                      restored.fine_potential, restored.full_records);
    std::ifstream copy(directory + "/copy.bin", std::ios::binary);
    EXPECT_EQ(original, Bytes(std::istreambuf_iterator<char>(copy), std::istreambuf_iterator<char>()));
}

TEST_F(SternheimerWeakQReferenceTest, RejectsWrongContractAndEachDimension)
{
    write();
    EXPECT_THROW(read_sternheimer_weak_q_reference(path, std::string(64, 'b'), dimensions), std::runtime_error);
    for (const auto d : {SternheimerWeakQReferenceDimensions{4, 2, 2},
                         SternheimerWeakQReferenceDimensions{3, 3, 2},
                         SternheimerWeakQReferenceDimensions{3, 2, 3}})
        EXPECT_THROW(read_sternheimer_weak_q_reference(path, contract, d), std::runtime_error);
}

TEST_F(SternheimerWeakQReferenceTest, CallerReplacesBothIndependentSCFCopiesBeforeHashing)
{
    write();
    auto independent_potential = potential;
    auto independent_records = records;
    independent_potential[0] = std::nextafter(potential[0], 2.0);
    independent_records[0].eigenvalues[0] = std::nextafter(records[0].eigenvalues[0], 0.0);
    independent_records[0].coefficients[0][0] = {0.25, 1e-16};
    auto reference = read();
    independent_potential.swap(reference.fine_potential);
    independent_records.swap(reference.full_records);
    expect_reals(independent_potential, potential);
    expect_reals(independent_records[0].eigenvalues, records[0].eigenvalues);
    expect_coefficients(independent_records[0].coefficients, records[0].coefficients);
    // A caller's previous reference-hash calculation belongs after both swaps.
    write_sternheimer_weak_q_reference(directory + "/copy.bin", contract, dimensions,
                                      independent_potential, independent_records);
    EXPECT_EQ(module_ri::sternheimer_siab::sha256_file(path),
              module_ri::sternheimer_siab::sha256_file(directory + "/copy.bin"));
}

TEST_F(SternheimerWeakQReferenceTest, RejectsEmptyContractsZeroDimensionsAndNonRegularFiles)
{
    EXPECT_THROW(write_sternheimer_weak_q_reference(path, "", dimensions, potential, records), std::runtime_error);
    EXPECT_THROW(write_sternheimer_weak_q_reference(path, std::string(65537, 'a'), dimensions,
                                                   potential, records), std::runtime_error);
    for (const auto d : {SternheimerWeakQReferenceDimensions{0, 2, 2},
                         SternheimerWeakQReferenceDimensions{3, 0, 2},
                         SternheimerWeakQReferenceDimensions{3, 2, 0}})
        EXPECT_THROW(write_sternheimer_weak_q_reference(path, contract, d, potential, records), std::runtime_error);
    EXPECT_THROW(read(), std::runtime_error);
    ASSERT_EQ(::mkfifo(path.c_str(), 0600), 0);
    EXPECT_THROW(read(), std::runtime_error);
    EXPECT_THROW(read_sternheimer_weak_q_reference(directory, contract, dimensions), std::runtime_error);
}

TEST_F(SternheimerWeakQReferenceTest, PreservesExtremeIndicesAndEmptyBandLists)
{
    records[0].local_k_index = std::numeric_limits<int>::min();
    records[0].global_k_index = std::numeric_limits<int>::max();
    records[1].eigenvalues.clear();
    records[1].occupations.clear();
    records[1].coefficients.clear();
    write();
    const auto restored = read();
    EXPECT_EQ(restored.full_records[0].local_k_index, std::numeric_limits<int>::min());
    EXPECT_EQ(restored.full_records[0].global_k_index, std::numeric_limits<int>::max());
    EXPECT_TRUE(restored.full_records[1].eigenvalues.empty());
    EXPECT_TRUE(restored.full_records[1].occupations.empty());
    EXPECT_TRUE(restored.full_records[1].coefficients.empty());
}

TEST_F(SternheimerWeakQReferenceTest, RejectsEveryTruncationCorruptionAndTrailingBytes)
{
    write();
    const auto original = bytes();
    for (std::size_t n = 0; n < original.size(); ++n)
    {
        SCOPED_TRACE(n);
        replace(Bytes(original.begin(), original.begin() + n));
        EXPECT_THROW(read(), std::runtime_error);
    }
    for (std::size_t i = 0; i < original.size(); ++i)
    {
        SCOPED_TRACE(i);
        auto changed = original;
        changed[i] ^= 1;
        replace(changed);
        EXPECT_THROW(read(), std::runtime_error);
    }
    auto trailing = original;
    trailing.push_back(0);
    replace(trailing);
    EXPECT_THROW(read(), std::runtime_error);
    trailing = original;
    trailing.insert(trailing.end() - 64, 0);
    reseal(trailing);
    replace(trailing);
    EXPECT_THROW(read(), std::runtime_error);
}

TEST_F(SternheimerWeakQReferenceTest, RejectsResealedUnknownSchemaBooleanCountsAndNonfiniteValues)
{
    write();
    const auto original = bytes();
    for (const auto offset : {std::size_t(0), std::size_t(8), record_offset() + 20, record_offset() + 21})
    {
        auto changed = original;
        changed[offset] = 42;
        reseal(changed);
        replace(changed);
        EXPECT_THROW(read(), std::runtime_error);
    }
    for (const auto offset : {std::size_t(12), dimension_offset(), dimension_offset() + 8,
                              dimension_offset() + 16, record_offset() + 78, record_offset() + 86})
    {
        auto changed = original;
        put64(changed, offset, std::numeric_limits<std::uint64_t>::max());
        reseal(changed);
        replace(changed);
        EXPECT_THROW(read(), std::runtime_error);
    }
    // Potential, both k vectors, weight, eigenvalue, occupation, complex parts, virtual eigenvalue.
    for (const auto offset : {potential_offset(), record_offset() + 22, record_offset() + 46,
                              record_offset() + 70, record_offset() + 94, record_offset() + 110,
                              record_offset() + 126, record_offset() + 134, record_offset() + 190,
                              record_offset() + 198, record_offset() + 206})
    {
        for (const auto bits : {UINT64_C(0x7ff0000000000000), UINT64_C(0x7ff8000000000001)})
        {
            auto changed = original;
            put64(changed, offset, bits);
            reseal(changed);
            replace(changed);
            EXPECT_THROW(read(), std::runtime_error);
        }
    }
}

TEST_F(SternheimerWeakQReferenceTest, RejectsHugeDimensionsEvenWhenCallerExpectsThem)
{
    write();
    const auto original = bytes();
    const auto huge = UINT64_C(1) << 40;
    for (int field = 0; field < 3; ++field)
    {
        auto changed = original;
        auto expected = dimensions;
        if (field == 0) expected.fine_grid_size = huge;
        if (field == 1) expected.full_record_count = huge;
        if (field == 2) expected.orbital_count = huge;
        put64(changed, dimension_offset() + 8 * field, huge);
        reseal(changed);
        replace(changed);
        EXPECT_THROW(read_sternheimer_weak_q_reference(path, contract, expected), std::runtime_error);
    }
    // Below the default band limit: rejection must also use available file bytes.
    auto changed = original;
    put64(changed, record_offset() + 78, 65536);
    reseal(changed);
    replace(changed);
    EXPECT_THROW(read(), std::runtime_error);
}

TEST_F(SternheimerWeakQReferenceTest, EnforcesReadAndWriteBudgetsWithoutPublishingPartialFiles)
{
    write();
    const auto original = bytes();
    for (int mode = 0; mode < 3; ++mode)
    {
        SternheimerWeakQReferenceLimits limits;
        if (mode == 0) limits.max_file_bytes = original.size() - 1;
        if (mode == 1) limits.max_allocation_bytes = sizeof(SternheimerLCAOOccupiedKPoint);
        if (mode == 2) limits.max_bands_per_record = 1;
        EXPECT_THROW(read_sternheimer_weak_q_reference(path, contract, dimensions, limits), std::runtime_error);
        EXPECT_THROW(write_sternheimer_weak_q_reference(directory + "/copy.bin", contract, dimensions,
                                                       potential, records, limits), std::runtime_error);
        EXPECT_NE(::access((directory + "/copy.bin").c_str(), F_OK), 0);
    }
    EXPECT_EQ(bytes(), original);
}

TEST_F(SternheimerWeakQReferenceTest, RejectsInvalidExporterData)
{
    records[0].occupations.pop_back();
    EXPECT_THROW(write(), std::runtime_error);
    records[0].occupations.push_back(0.25);
    records[0].coefficients[0].pop_back();
    EXPECT_THROW(write(), std::runtime_error);
    records[0].coefficients[0].push_back({1.5, -2});
    records[0].unoccupied_coefficients.clear();
    EXPECT_THROW(write(), std::runtime_error);
    records[0].unoccupied_coefficients = {{{7, 8}, {-9, -10}}};
    records[0].grid_kpoint[2] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(write(), std::runtime_error);
    records[0].grid_kpoint[2] = 1.0;
    potential[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(write(), std::runtime_error);
    potential[0] = 1.25;
    records[0].coefficients[0][0] = {0, std::numeric_limits<double>::infinity()};
    EXPECT_THROW(write(), std::runtime_error);
    EXPECT_NE(::access(path.c_str(), F_OK), 0);
}

TEST_F(SternheimerWeakQReferenceTest, NeverOverwritesExistingFileOrSymlink)
{
    write();
    const auto original = bytes();
    potential[0] = 99;
    EXPECT_THROW(write(), std::runtime_error);
    EXPECT_EQ(bytes(), original);
    ASSERT_EQ(::symlink("missing-target", (directory + "/copy.bin").c_str()), 0);
    EXPECT_THROW(write_sternheimer_weak_q_reference(directory + "/copy.bin", contract, dimensions,
                                                   potential, records), std::runtime_error);
    struct stat st;
    ASSERT_EQ(::lstat((directory + "/copy.bin").c_str(), &st), 0);
    EXPECT_TRUE(S_ISLNK(st.st_mode));
}

TEST_F(SternheimerWeakQReferenceTest, ConcurrentPublishHasExactlyOneWinner)
{
    const pid_t child = ::fork();
    ASSERT_GE(child, 0);
    if (child == 0)
    {
        try { write(); ::_exit(0); }
        catch (...) { ::_exit(1); }
    }
    int parent_success = 0;
    try { write(); parent_success = 1; }
    catch (const std::runtime_error&) {}
    int status = 0;
    ASSERT_EQ(::waitpid(child, &status, 0), child);
    ASSERT_TRUE(WIFEXITED(status));
    EXPECT_EQ(parent_success + (WEXITSTATUS(status) == 0 ? 1 : 0), 1);
    EXPECT_NO_THROW(read());
}
} // namespace
