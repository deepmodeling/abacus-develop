#ifndef STERNHEIMER_WEAK_MATRIX_AUDIT_H
#define STERNHEIMER_WEAK_MATRIX_AUDIT_H

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ModuleRI
{

// Pure ownership policy. No environment access, sampling, or MPI communication.
class SternheimerWeakAuditShard
{
  public:
    const int index;
    const int count;

    static SternheimerWeakAuditShard parse(const char* raw_index, const char* raw_count, const int channels)
    {
        const int index = parse_integer(raw_index, 0, "WEAK_AUDIT_SHARD_INDEX");
        const int count = parse_integer(raw_count, 1, "WEAK_AUDIT_SHARD_COUNT");
        if (channels <= 0 || count <= 0 || count > channels || index >= count)
            throw std::invalid_argument("Weak matrix audit requires 0 <= shard index < count <= channels.");
        return SternheimerWeakAuditShard(index, count, channels);
    }

    bool owns(const int column) const
    {
        return column >= 0 && column < channels_ && column % count == index;
    }

    std::vector<int> columns() const
    {
        std::vector<int> result;
        result.reserve((channels_ - 1 - index) / count + 1);
        // Widen the increment so even an INT_MAX-sized space cannot wrap.
        for (long long column = index; column < channels_; column += count)
            result.push_back(static_cast<int>(column));
        return result;
    }

  private:
    const int channels_;

    SternheimerWeakAuditShard(const int index, const int count, const int channels)
        : index(index), count(count), channels_(channels)
    {
    }

    static int parse_integer(const char* raw, const int fallback, const char* name)
    {
        if (raw == nullptr) return fallback;
        const auto invalid = [name]() {
            return std::invalid_argument(std::string(name) + " must be a nonempty unsigned decimal integer.");
        };
        if (*raw == '\0') throw invalid();
        int value = 0;
        for (const char* p = raw; *p != '\0'; ++p)
        {
            if (*p < '0' || *p > '9') throw invalid();
            const int digit = *p - '0';
            if (value > (std::numeric_limits<int>::max() - digit) / 10) throw invalid();
            value = value * 10 + digit;
        }
        return value;
    }
};

} // namespace ModuleRI

#endif
