/**
 * @file stru_frame_io.cpp
 * @brief Implementation of the cumulative STRU frame extraction helpers.
 */
#include "stru_frame_io.h"

#include <cstdlib>
#include <fstream>
#include <string>

namespace
{
const char* const frame_marker = "# MD step ";

/// Parse the step number following a "# MD step " marker; 0 when malformed.
int parse_marker_step(const std::string& line)
{
    const std::string number = line.substr(std::string(frame_marker).size());
    if (number.empty())
    {
        return 0;
    }
    char* end = nullptr;
    const long value = std::strtol(number.c_str(), &end, 10);
    if (end == number.c_str() || value <= 0)
    {
        return 0;
    }
    return static_cast<int>(value);
}
} // namespace

namespace stru_frame_io
{
bool extract_last_frame(const std::string& cumulative_path, const std::string& frame_path)
{
    std::ifstream ifs(cumulative_path.c_str());
    if (!ifs)
    {
        return false;
    }
    std::string line;
    std::string last_frame;
    std::string current_frame;
    bool in_frame = false;
    while (std::getline(ifs, line))
    {
        const bool is_marker = line.compare(0, std::string(frame_marker).size(), frame_marker) == 0;
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

bool extract_frame_by_step(const std::string& cumulative_path, int step, const std::string& frame_path)
{
    if (step <= 0)
    {
        return false;
    }
    std::ifstream ifs(cumulative_path.c_str());
    if (!ifs)
    {
        return false;
    }
    std::string line;
    std::string current_frame;
    bool in_frame = false;
    int current_step = 0;
    while (std::getline(ifs, line))
    {
        const bool is_marker = line.compare(0, std::string(frame_marker).size(), frame_marker) == 0;
        if (is_marker)
        {
            if (in_frame && current_step == step)
            {
                // we just passed the end of the frame we were looking for
                break;
            }
            current_step = parse_marker_step(line);
            current_frame.clear();
            in_frame = true;
            continue;
        }
        if (in_frame && current_step == step)
        {
            current_frame += line;
            current_frame += '\n';
        }
    }
    if (!in_frame || current_step != step || current_frame.empty())
    {
        return false;
    }
    std::ofstream ofs(frame_path.c_str(), std::ios::out | std::ios::trunc);
    if (!ofs)
    {
        return false;
    }
    ofs << current_frame;
    return ofs.good();
}

int last_frame_step(const std::string& cumulative_path)
{
    std::ifstream ifs(cumulative_path.c_str());
    if (!ifs)
    {
        return 0;
    }
    int last_step = 0;
    std::string line;
    while (std::getline(ifs, line))
    {
        const bool is_marker = line.compare(0, std::string(frame_marker).size(), frame_marker) == 0;
        if (is_marker)
        {
            const int step = parse_marker_step(line);
            if (step > 0)
            {
                last_step = step;
            }
        }
    }
    return last_step;
}
} // namespace stru_frame_io
