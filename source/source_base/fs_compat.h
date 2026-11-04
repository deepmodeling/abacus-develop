#ifndef MODULEBASE_FS_COMPAT_H
#define MODULEBASE_FS_COMPAT_H

//==========================================================
// Small filesystem-portability helpers.
//
// The POSIX `mkdir(path, mode)` takes a permission-mode argument that
// does not exist in the Windows CRT (`_mkdir`/MinGW `mkdir` take only a
// path). This header provides a single cross-platform directory-creation
// helper so call sites stay identical on every platform.
//==========================================================

#include <cerrno>
#include <string>

#ifdef _WIN32
#include <direct.h> // _mkdir
#else
#include <sys/stat.h>  // mkdir
#include <sys/types.h>
#endif

namespace ModuleBase
{

/**
 * @brief Create a single directory, portably.
 *
 * @param path directory path to create
 * @return 0 on success; -1 on failure with `errno` set (e.g. EEXIST when
 *         the directory already exists), matching POSIX `mkdir` semantics.
 *
 * On Windows the permission mode is not applicable and is ignored; on
 * POSIX systems the directory is created with mode 0755 (subject to umask),
 * preserving the previous behaviour of the call sites.
 */
inline int make_directory(const std::string& path)
{
#ifdef _WIN32
    return _mkdir(path.c_str());
#else
    return mkdir(path.c_str(), 0755);
#endif
}

/**
 * @brief Ensure a directory exists, creating it lazily if needed.
 *
 * Intended for write-time directory creation: output writers call this just
 * before writing a file into a directory that may not exist yet, instead of
 * relying on a central startup routine to predict which directories will be
 * needed. An already-existing directory (EEXIST) is not an error, so callers
 * on different MPI ranks may invoke it concurrently for the same path. The
 * parent directory must already exist (only one level is created).
 *
 * @param path directory path to ensure; an empty path is accepted as a no-op
 * @return true if the directory now exists (created or already present),
 *         false on a real creation failure (errno is set)
 */
inline bool ensure_directory(const std::string& path)
{
    if (path.empty())
    {
        return true;
    }
    const int ret = make_directory(path);
    return ret == 0 || errno == EEXIST;
}

} // namespace ModuleBase

#endif // MODULEBASE_FS_COMPAT_H
