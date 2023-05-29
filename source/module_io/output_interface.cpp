#include "module_io/output_interface.h"

namespace ModuleIO
{
void Output_Interface::setFilename(const std::string& fn)
{
    _fn = fn;
}
void Output_Interface::setPrecision(const int precision)
{
    _precision = precision;
}
} // namespace ModuleIO