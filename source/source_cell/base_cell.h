#ifndef BASE_CELL_H
#define BASE_CELL_H

#include "source_base/matrix3.h"

#include <cstdint>

class BaseCell
{
public:
    enum class Kind
    {
        unit_cell,
        md_cell
    };

    virtual ~BaseCell() = default;

    Kind kind() const
    {
        return get_kind();
    }

    std::int64_t nat() const
    {
        return get_nat();
    }

    double lat0() const
    {
        return get_lat0();
    }

    double omega() const
    {
        return get_omega();
    }

    const ModuleBase::Matrix3& latvec() const
    {
        return get_latvec();
    }

    const ModuleBase::Matrix3& GT() const
    {
        return get_GT();
    }

    void require_kind(const Kind& expected, const char* caller) const;

private:
    virtual Kind get_kind() const = 0;
    virtual std::int64_t get_nat() const = 0;
    virtual double get_lat0() const = 0;
    virtual double get_omega() const = 0;
    virtual const ModuleBase::Matrix3& get_latvec() const = 0;
    virtual const ModuleBase::Matrix3& get_GT() const = 0;
};

#endif
