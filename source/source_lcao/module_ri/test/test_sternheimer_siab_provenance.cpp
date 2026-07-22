#include "source_lcao/module_ri/sternheimer_siab_provenance.h"

#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

namespace siab = module_ri::sternheimer_siab;

class TemporaryFile
{
  public:
    TemporaryFile(const std::string& path, const std::string& contents) : path_(path)
    {
        std::ofstream out(path_.c_str(), std::ios::binary | std::ios::trunc);
        out.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    }

    ~TemporaryFile()
    {
        std::remove(path_.c_str());
    }

    const std::string& path() const
    {
        return path_;
    }

  private:
    std::string path_;
};

TEST(SternheimerSIABProvenance, HashesKnownBytesAndFile)
{
    const std::string expected = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
    const std::vector<unsigned char> bytes{'a', 'b', 'c'};
    EXPECT_EQ(siab::sha256_bytes(bytes), expected);

    const TemporaryFile file("sternheimer_siab_sha256_abc.tmp", "abc");
    EXPECT_EQ(siab::sha256_file(file.path()), expected);
    EXPECT_EQ(siab::sha256_file_manifest({file.path()}), expected);

    siab::Sha256 chunked;
    const unsigned char a = 'a';
    const unsigned char bc[2] = {'b', 'c'};
    chunked.update(&a, 1);
    chunked.update(bc, 2);
    EXPECT_EQ(chunked.finish(), expected);
}

TEST(SternheimerSIABProvenance, FileManifestIsOrderedAndRejectsMissingInput)
{
    const TemporaryFile first("sternheimer_siab_sha256_first.tmp", "first");
    const TemporaryFile second("sternheimer_siab_sha256_second.tmp", "second");
    EXPECT_NE(siab::sha256_file_manifest({first.path(), second.path()}),
              siab::sha256_file_manifest({second.path(), first.path()}));
    EXPECT_THROW(siab::sha256_file_manifest({}), std::invalid_argument);
    EXPECT_THROW(siab::sha256_file("sternheimer_siab_missing_file.tmp"), std::runtime_error);
}

TEST(SternheimerSIABProvenance, UniqueFileManifestIgnoresRepeatedAuxiliaryBasisContents)
{
    const TemporaryFile first("sternheimer_siab_sha256_unique_first.tmp", "same auxiliary basis");
    const TemporaryFile duplicate("sternheimer_siab_sha256_unique_duplicate.tmp", "same auxiliary basis");
    const TemporaryFile second("sternheimer_siab_sha256_unique_second.tmp", "different auxiliary basis");

    const std::string single_hash = siab::sha256_file(first.path());
    EXPECT_EQ(siab::sha256_unique_file_manifest({first.path(), first.path()}), single_hash);
    EXPECT_EQ(siab::sha256_unique_file_manifest({first.path(), duplicate.path()}), single_hash);
    EXPECT_EQ(siab::sha256_unique_file_manifest({first.path(), duplicate.path(), second.path()}),
              siab::sha256_file_manifest({first.path(), second.path()}));
    EXPECT_THROW(siab::sha256_unique_file_manifest({}), std::invalid_argument);
}

TEST(SternheimerSIABProvenance, ResolvesExplicitFilesWithoutPlaceholders)
{
    const TemporaryFile file("sternheimer_siab_resolve.tmp", "x");
    EXPECT_EQ(siab::resolve_required_input_files("", {file.path()}, "orbital"),
              (std::vector<std::string>{file.path()}));
    EXPECT_THROW(siab::resolve_required_input_files("", {}, "ABFS"), std::runtime_error);
    EXPECT_THROW(siab::resolve_required_input_files("", {"auto"}, "pseudopotential"), std::runtime_error);
}

TEST(SternheimerSIABProvenance, RequiresExactlyOneExplicitPrimitiveRcut)
{
    EXPECT_DOUBLE_EQ(siab::require_single_primitive_rcut({8.0}), 8.0);
    EXPECT_THROW(siab::require_single_primitive_rcut({}), std::invalid_argument);
    EXPECT_THROW(siab::require_single_primitive_rcut({8.0, 10.0}), std::invalid_argument);
    EXPECT_THROW(siab::require_single_primitive_rcut({0.0}), std::invalid_argument);
}

TEST(SternheimerSIABProvenance, RequiresFullCompiledCommitMetadata)
{
    const std::string commit = "c407b55f9dc4ae83377b4c37b4f3656fa91f447a";
    EXPECT_EQ(siab::require_source_commit(commit + " (2026-07-18)"), commit);
    EXPECT_THROW(siab::require_source_commit("c407b55f9"), std::runtime_error);
    EXPECT_THROW(siab::require_source_commit("unknown"), std::runtime_error);
}

} // namespace
