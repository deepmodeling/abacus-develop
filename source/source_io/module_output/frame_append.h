#ifndef FRAME_APPEND_H
#define FRAME_APPEND_H

#include <stdexcept>
#include <string>

namespace ModuleIO
{

/// @brief Tell whether an output written for ionic step `istep` should append.
///
/// Mode 1 appends every frame after the first ionic step, keeping one trajectory
/// file. Mode 0 truncates to the latest snapshot and mode 2 writes a separate
/// file per step, so neither of them appends.
///
/// @param out_app_flag Output append mode, must be 0, 1 or 2.
/// @param istep Zero-based ionic step index, must be non-negative.
/// @return true when the frame must be appended to an existing file.
/// @throws std::invalid_argument if out_app_flag or istep is out of range.
inline bool frame_appends(const int out_app_flag, const int istep)
{
    if (out_app_flag < 0 || out_app_flag > 2)
    {
        throw std::invalid_argument("frame_appends: out_app_flag must be 0, 1 or 2, got "
                                    + std::to_string(out_app_flag));
    }
    if (istep < 0)
    {
        throw std::invalid_argument("frame_appends: istep must be non-negative, got "
                                    + std::to_string(istep));
    }

    return out_app_flag == 1 && istep > 0;
}

} // namespace ModuleIO

#endif // FRAME_APPEND_H
