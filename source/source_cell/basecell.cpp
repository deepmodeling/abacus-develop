#include "source_cell/basecell.h"

#include "source_base/tool_quit.h"

void BaseCell::require_kind(const Kind& expected, const char* caller) const
{
    if (this->kind() != expected)
    {
        const char* required_cell = expected == Kind::unitcell ? "UnitCell" : "MDCell";
        ModuleBase::WARNING_QUIT(caller, std::string("This operation only supports ") + required_cell + ".");
    }
}
