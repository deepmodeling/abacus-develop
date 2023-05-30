#ifndef OUTPUT_INTERFACE_H
#define OUTPUT_INTERFACE_H

namespace ModuleIO
{
/// @brief base class for interface output classes
class Output_Interface
{
  public:
    virtual ~Output_Interface()
    {
    }
    virtual void write() = 0;
};
} // namespace ModuleIO

#endif