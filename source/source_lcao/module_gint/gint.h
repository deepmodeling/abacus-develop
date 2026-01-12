#pragma once
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

    protected:
    static GintInfo* gint_info_;
};

}