#include "module_base/element_name.h"
#include "module_base/global_variable.h"
#include "module_io/cube_io.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

#include <algorithm>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
int choose_cube_rows_per_chunk(const int nz, const int precision)
{
    constexpr std::size_t target_chunk_bytes = 8 * 1024 * 1024;
    const std::size_t approx_bytes_per_value = static_cast<std::size_t>(std::max(precision + 9, 16));
    const std::size_t approx_line_breaks = static_cast<std::size_t>(std::max(nz / 6, 1));
    const std::size_t approx_bytes_per_row = static_cast<std::size_t>(std::max(nz, 1)) * approx_bytes_per_value
                                             + approx_line_breaks + 1;
    const std::size_t rows = target_chunk_bytes / std::max<std::size_t>(approx_bytes_per_row, 1);
    return static_cast<int>(std::max<std::size_t>(rows, 1));
}

// Estimates buffer size needed for a chunk of cube data
std::size_t estimate_cube_chunk_size(const int row_count,
                                      const int nz,
                                      const int ndata_line,
                                      const int precision)
{
    const std::size_t bytes_per_value = static_cast<std::size_t>(std::max(precision + 16, 32));
    const int inner_line_breaks = (ndata_line > 0 && nz > 0) ? (nz - 1) / ndata_line : 0;
    return static_cast<std::size_t>(row_count)
           * (static_cast<std::size_t>(nz) * bytes_per_value
              + static_cast<std::size_t>(inner_line_breaks + 1));
}

void format_cube_data_rows_to_buffer(const std::vector<double>& data,
                                     const int row_begin,
                                     const int row_end,
                                     const int nz,
                                     const int precision,
                                     const int ndata_line,
                                     std::string& buffer)
{
    const int row_count = row_end - row_begin;
    buffer.resize(estimate_cube_chunk_size(row_count, nz, ndata_line, precision));
    if (buffer.empty())
    {
        return;
    }

    char* ptr = &buffer[0];
    char* end = ptr + buffer.size();

    const auto ensure_space = [&buffer, &ptr, &end](const std::size_t required) {
        const std::size_t used = static_cast<std::size_t>(ptr - &buffer[0]);
        if (static_cast<std::size_t>(end - ptr) >= required)
        {
            return;
        }

        const std::size_t grown_size = std::max(buffer.size() + required, buffer.size() * 2 + 1);
        buffer.resize(grown_size);
        ptr = &buffer[0] + used;
        end = &buffer[0] + buffer.size();
    };

    for (int ixy = row_begin; ixy < row_end; ++ixy)
    {
        for (int iz = 0; iz < nz; ++iz)
        {
            ensure_space(static_cast<std::size_t>(std::max(precision + 16, 32)));
            int written = std::snprintf(ptr,
                                        static_cast<std::size_t>(end - ptr),
                                        " %.*e",
                                        precision,
                                        data[ixy * nz + iz]);
            if (written < 0)
            {
                written = 0;
            }
            if (written >= end - ptr)
            {
                ensure_space(static_cast<std::size_t>(written) + 1);
                written = std::snprintf(ptr,
                                        static_cast<std::size_t>(end - ptr),
                                        " %.*e",
                                        precision,
                                        data[ixy * nz + iz]);
                if (written < 0)
                {
                    written = 0;
                }
            }
            ptr += written;

            if ((iz + 1) % ndata_line == 0 && iz != nz - 1)
            {
                ensure_space(1);
                *ptr++ = '\n';
            }
        }
        ensure_space(1);
        *ptr++ = '\n';
    }

    buffer.resize(static_cast<std::size_t>(ptr - &buffer[0]));
}

std::string format_cube_data_chunk(const std::vector<double>& data,
                                   const int row_begin,
                                   const int row_end,
                                   const int nz,
                                   const int precision,
                                   const int ndata_line)
{
    std::string data_buffer;
    format_cube_data_rows_to_buffer(data, row_begin, row_end, nz, precision, ndata_line, data_buffer);
    return data_buffer;
}

std::vector<std::string> format_cube_data_chunk_parallel(const std::vector<double>& data,
                                                         const int row_begin,
                                                         const int row_end,
                                                         const int nz,
                                                         const int precision,
                                                         const int ndata_line)
{
#ifdef _OPENMP
    const int row_count = row_end - row_begin;
    const int max_threads = std::max(omp_get_max_threads(), 1);
    const int subchunk_count = std::min(max_threads, row_count);
    if (subchunk_count <= 1)
    {
        return std::vector<std::string>{format_cube_data_chunk(data, row_begin, row_end, nz, precision, ndata_line)};
    }

    std::vector<std::string> subchunk_buffers(subchunk_count);
    int actual_threads = 1;

#pragma omp parallel num_threads(subchunk_count)
    {
        const int thread_id = omp_get_thread_num();
        const int thread_count = omp_get_num_threads();
#pragma omp single
        {
            actual_threads = thread_count;
        }

        const int sub_begin = row_begin + row_count * thread_id / thread_count;
        const int sub_end = row_begin + row_count * (thread_id + 1) / thread_count;
        format_cube_data_rows_to_buffer(data,
                                        sub_begin,
                                        sub_end,
                                        nz,
                                        precision,
                                        ndata_line,
                                        subchunk_buffers[thread_id]);
    }
    subchunk_buffers.resize(static_cast<std::size_t>(actual_threads));
    return subchunk_buffers;
#else
    return std::vector<std::string>{format_cube_data_chunk(data, row_begin, row_end, nz, precision, ndata_line)};
#endif
}

class AsyncCubeWriter
{
  public:
    explicit AsyncCubeWriter(std::ofstream& ofs, const std::size_t max_queue_size)
        : ofs_(ofs), max_queue_size_(std::max<std::size_t>(max_queue_size, 1)), worker_(&AsyncCubeWriter::run, this)
    {
    }

    ~AsyncCubeWriter()
    {
        finish();
    }

    void push(std::vector<std::string>&& buffers)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_not_full_.wait(lock, [this] { return queue_.size() < max_queue_size_ || finished_; });
        if (finished_)
        {
            return;
        }
        queue_.push_back(std::move(buffers));
        cv_not_empty_.notify_one();
    }

    void finish()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (finished_)
            {
                return;
            }
            finished_ = true;
        }
        cv_not_empty_.notify_one();
        cv_not_full_.notify_all();
        if (worker_.joinable())
        {
            worker_.join();
        }
    }

  private:
    void run()
    {
        while (true)
        {
            std::vector<std::string> buffers;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_not_empty_.wait(lock, [this] { return finished_ || !queue_.empty(); });
                if (queue_.empty())
                {
                    if (finished_)
                    {
                        break;
                    }
                    continue;
                }

                buffers = std::move(queue_.front());
                queue_.pop_front();
                cv_not_full_.notify_one();
            }

            for (const std::string& buffer : buffers)
            {
                ofs_.write(buffer.data(), buffer.size());
            }
        }
    }

    std::ofstream& ofs_;
    const std::size_t max_queue_size_;
    std::deque<std::vector<std::string>> queue_;
    std::mutex mutex_;
    std::condition_variable cv_not_empty_;
    std::condition_variable cv_not_full_;
    bool finished_ = false;
    std::thread worker_;
};
} // namespace

void ModuleIO::write_vdata_palgrid(
    const Parallel_Grid& pgrid,
    const double* const data,
    const int is,
    const int nspin,
    const int iter,
    const std::string& fn,
    const double ef,
    const UnitCell*const ucell,
    const int precision,
    const int out_fermi)
{
    ModuleBase::TITLE("ModuleIO", "write_vdata_palgrid");

    const int my_rank = GlobalV::MY_RANK;
    const int my_pool = GlobalV::MY_POOL;

    time_t start;
    time_t end;
    std::stringstream ss;

    const int& nx = pgrid.nx;
    const int& ny = pgrid.ny;
    const int& nz = pgrid.nz;
    const int& nxyz = nx * ny * nz;

    start = time(nullptr);

    // reduce
    std::vector<double> data_xyz_full(nxyz);    // data to be written
#ifdef __MPI    // reduce to rank 0
    if (my_pool == 0 && GlobalV::MY_STOGROUP == 0)
    {
        pgrid.reduce(data_xyz_full.data(), data);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    std::memcpy(data_xyz_full.data(), data, nxyz * sizeof(double));
#endif

    // build the info structure
    if (my_rank == 0)
    {
        /// output header for cube file
        ss << "STEP: " << iter << "  Cubefile created from ABACUS. Inner loop is z, followed by y and x" << std::endl;

        ss << nspin << " (nspin) ";
        ss << std::fixed;
        ss << std::setprecision(6);
        if (out_fermi == 1)
        {
            if (PARAM.globalv.two_fermi)
            {
                if (is == 0)
                {
                    ss << ef << " (fermi energy for spin=1, in Ry)" << std::endl;
                }
                else if (is == 1)
                {
                    ss << ef << " (fermi energy for spin=2, in Ry)" << std::endl;
                }
            }
            else
            {
                ss << ef << " (fermi energy, in Ry)" << std::endl;
            }
        }
        else
        {
            ss << std::endl;
        }

        std::vector<std::string> comment(2);
        for (int i = 0;i < 2;++i) { std::getline(ss, comment[i]); }

        double fac = ucell->lat0;
        std::vector<double> dx = { fac * ucell->latvec.e11 / double(nx), fac * ucell->latvec.e12 / double(nx), fac * ucell->latvec.e13 / double(nx) };
        std::vector<double> dy = { fac * ucell->latvec.e21 / double(ny), fac * ucell->latvec.e22 / double(ny), fac * ucell->latvec.e23 / double(ny) };
        std::vector<double> dz = { fac * ucell->latvec.e31 / double(nz), fac * ucell->latvec.e32 / double(nz), fac * ucell->latvec.e33 / double(nz) };

        std::string element = "";
        std::vector<int> atom_type;
        std::vector<double> atom_charge;
        std::vector<std::vector<double>> atom_pos;
        for (int it = 0; it < ucell->ntype; it++)
        {
            // erase the number in label, such as Fe1.
            element = ucell->atoms[it].label;
            std::string::iterator temp = element.begin();
            while (temp != element.end())
            {
                if ((*temp >= '1') && (*temp <= '9'))
                {
                    temp = element.erase(temp);
                }
                else
                {
                    temp++;
                }
            }

            for (int ia = 0; ia < ucell->atoms[it].na; ia++)
            {
                // convert from label to atomic number
                int z = 0;
                for (int j = 0; j != ModuleBase::element_name.size(); j++)
                {
                    if (element == ModuleBase::element_name[j])
                    {
                        z = j + 1;
                        break;
                    }
                }
                atom_type.push_back(z);
                atom_charge.push_back(ucell->atoms[it].ncpp.zv);
                atom_pos.push_back({ fac * ucell->atoms[it].tau[ia].x, fac * ucell->atoms[it].tau[ia].y, fac * ucell->atoms[it].tau[ia].z });
            }
        }
        write_cube(fn, comment, ucell->nat, { 0.0, 0.0, 0.0 }, nx, ny, nz, dx, dy, dz, atom_type, atom_charge, atom_pos, data_xyz_full, precision);

        end = time(nullptr);
        ModuleBase::GlobalFunc::OUT_TIME("write_vdata_palgrid", start, end);
    }

    return;
}

void ModuleIO::write_cube(const std::string& file,
    const std::vector<std::string>& comment,
    const int& natom,
    const std::vector<double>& origin,
    const int& nx,
    const int& ny,
    const int& nz,
    const std::vector<double>& dx,
    const std::vector<double>& dy,
    const std::vector<double>& dz,
    const std::vector<int>& atom_type,
    const std::vector<double>& atom_charge,
    const std::vector<std::vector<double>>& atom_pos,
    const std::vector<double>& data,
    const int precision,
    const int ndata_line)
{
    assert(comment.size() >= 2);
    for (int i = 0;i < 2;++i) { assert(comment[i].find("\n") == std::string::npos); }
    assert(origin.size() >= 3);
    assert(dx.size() >= 3);
    assert(dy.size() >= 3);
    assert(dz.size() >= 3);
    assert(atom_type.size() >= natom);
    assert(atom_charge.size() >= natom);
    assert(atom_pos.size() >= natom);
    for (int i = 0;i < natom;++i) { assert(atom_pos[i].size() >= 3); }
    assert(data.size() >= nx * ny * nz);
    assert(precision >= 0);
    assert(ndata_line > 0);

    std::ofstream ofs(file);

    std::ostringstream header_stream;
    for (int i = 0;i < 2;++i) { header_stream << comment[i] << "\n"; }

    header_stream << std::fixed;
    header_stream << std::setprecision(1);    // as before

    header_stream << natom << " " << origin[0] << " " << origin[1] << " " << origin[2] << " \n";

    header_stream << std::setprecision(6);    //as before
    header_stream << nx << " " << dx[0] << " " << dx[1] << " " << dx[2] << "\n";
    header_stream << ny << " " << dy[0] << " " << dy[1] << " " << dy[2] << "\n";
    header_stream << nz << " " << dz[0] << " " << dz[1] << " " << dz[2] << "\n";

    for (int i = 0;i < natom;++i)
    {
        header_stream << " " << atom_type[i] << " " << atom_charge[i] << " " << atom_pos[i][0] << " " << atom_pos[i][1] << " " << atom_pos[i][2] << "\n";
    }
    const std::string header_buffer = header_stream.str();

    ofs.write(header_buffer.data(), header_buffer.size());

    const int nxy = nx * ny;
    const int rows_per_chunk = choose_cube_rows_per_chunk(nz, precision);
    AsyncCubeWriter async_writer(ofs, 2);
    for (int row_begin = 0; row_begin < nxy; row_begin += rows_per_chunk)
    {
        const int row_end = std::min(row_begin + rows_per_chunk, nxy);

        std::vector<std::string> data_buffers = format_cube_data_chunk_parallel(data,
                                                                                row_begin,
                                                                                row_end,
                                                                                nz,
                                                                                precision,
                                                                                ndata_line);
        async_writer.push(std::move(data_buffers));
    }

    async_writer.finish();

    ofs.close();
}
