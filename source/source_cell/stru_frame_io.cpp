/**
 * @file stru_frame_io.cpp
 * @brief Implementation of the cumulative STRU frame extraction helpers.
 */
#include "stru_frame_io.h"

#include <fstream>
#include <string>

namespace stru_frame_io
{
bool extract_last_frame(const std::string& cumulative_path, const std::string& frame_path)
{
    std::ifstream ifs(cumulative_path.c_str());
    if (!ifs)
    {
        return false;
    }
    const std::string frame_marker = "# MD step ";
    std::string line;
    std::string last_frame;
    std::string current_frame;
    bool in_frame = false;
    while (std::getline(ifs, line))
    {
        const bool is_marker = line.compare(0, frame_marker.size(), frame_marker) == 0;
        if (is_marker)
        {
            if (in_frame)
            {
                last_frame = current_frame;
            }
            current_frame.clear();
            in_frame = true;
            continue;
        }
        if (in_frame)
        {
            current_frame += line;
            current_frame += '\n';
        }
    }
    if (in_frame)
    {
        last_frame = current_frame;
    }
    if (last_frame.empty())
    {
        return false;
    }
    std::ofstream ofs(frame_path.c_str(), std::ios::out | std::ios::trunc);
    if (!ofs)
    {
        return false;
    }
    ofs << last_frame;
    return ofs.good();
}
} // namespace stru_frame_io
