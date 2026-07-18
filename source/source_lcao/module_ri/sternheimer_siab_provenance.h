#ifndef STERNHEIMER_SIAB_PROVENANCE_H
#define STERNHEIMER_SIAB_PROVENANCE_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace module_ri
{
namespace sternheimer_siab
{

class Sha256
{
  public:
    Sha256();

    void update(const unsigned char* data, std::size_t size);
    std::string finish();

  private:
    void transform(const unsigned char* block);

    std::array<std::uint32_t, 8> state_;
    std::array<unsigned char, 64> buffer_;
    std::size_t buffer_size_ = 0;
    std::uint64_t total_size_ = 0;
    bool finished_ = false;
};

std::string sha256_bytes(const std::vector<unsigned char>& bytes);
std::string sha256_file(const std::string& path);
std::string sha256_file_manifest(const std::vector<std::string>& paths);

std::vector<std::string> resolve_required_input_files(const std::string& directory,
                                                      const std::vector<std::string>& filenames,
                                                      const std::string& label);

double require_single_primitive_rcut(const std::vector<double>& rcuts);

std::string resolve_executable_path();
std::string require_source_commit(const std::string& compiled_commit);

} // namespace sternheimer_siab
} // namespace module_ri

#endif
