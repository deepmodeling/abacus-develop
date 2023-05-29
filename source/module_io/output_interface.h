#ifndef OUTPUT_INTERFACE_H
#define OUTPUT_INTERFACE_H

#include <string>

namespace ModuleIO
{
class Output_Interface
{
  public:
    virtual ~Output_Interface()
    {
    }
    virtual void setFilename(const std::string& fn);
    virtual void setPrecision(const int precision);
    virtual void write() = 0;

  public:
    std::string _fn;
    int _precision;
};
} // namespace ModuleIO

#endif