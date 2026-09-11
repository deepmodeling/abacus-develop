#ifndef BASECELL_H
#define BASECELL_H

#include <cstdint>

class BaseCell
{
public:
    enum class Kind
    {
        unitcell,
        mdcell
    };

    virtual ~BaseCell() = default;

    virtual Kind kind() const = 0;

    void require_kind(const Kind& expected, const char* caller) const;
};

#endif
