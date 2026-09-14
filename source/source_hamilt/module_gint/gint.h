#pragma once
#include <cassert>
#include <memory>
#include "gint_info.h"
#include "gint_type.h"

namespace ModuleGint
{

class Gint
{
    public:
    Gint() = default;
    virtual ~Gint() = default;

    // note that gint_info_ is a static member variable
    // it is shared by all instances of Gint
    static void set_gint_info(GintInfo* gint_info)
    {
        gint_info_ = gint_info;
    }

    static const GintInfo& get_gint_info()
    {
        // set_gint_info() must have been called by the owning ESolver before any
        // grid integration runs; dereferencing a null gint_info_ here would be UB.
        assert(gint_info_ != nullptr && "Gint::set_gint_info() has not been called");
        return *gint_info_;
    }

    protected:
    static GintInfo* gint_info_;
};

}