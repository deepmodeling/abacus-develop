#ifndef STRU_FRAME_IO_H
#define STRU_FRAME_IO_H

#include <string>

/**
 * @file stru_frame_io.h
 * @brief Helpers to read frames from fixed-name cumulative STRU files.
 *
 * In out_app_flag == 1 mode run_md appends one structure frame per output
 * step to a single "STRU" file, each frame introduced by a one-based
 * "# MD step N" comment line. Restart readers such as MDCellReader::read_stru
 * expect a single-frame STRU file, so the frames written this way must be
 * split before they can be consumed. This header owns that frame-splitting
 * knowledge so callers outside source_cell do not parse STRU text themselves.
 */
namespace stru_frame_io
{
/// @brief Extract the last MD frame from a fixed-name cumulative STRU file.
/// @param cumulative_path path of the mode-1 "STRU" file holding one frame per
///        output step, each introduced by a "# MD step N" comment line
/// @param frame_path path of the single-frame file to (over)write with the
///        last frame; the caller passes this path to the structure readers
/// @return true when at least one frame was found and written
/// @details The last frame is copied verbatim, without its "# MD step N"
/// comment line, into frame_path. The side file is rewritten on every call so
/// repeated restarts always see the most recent frame.
bool extract_last_frame(const std::string& cumulative_path, const std::string& frame_path);
} // namespace stru_frame_io

#endif // STRU_FRAME_IO_H
