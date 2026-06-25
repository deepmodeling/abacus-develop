#ifndef DIAGO_TRACE_H_
#define DIAGO_TRACE_H_

#include "source_base/global_variable.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>

namespace hsolver
{

class DiagoTrace
{
  public:
    explicit DiagoTrace(const std::string& solver_name)
        : enabled_(is_enabled()), solver_name_(solver_name)
    {
        if (!this->enabled_)
        {
            return;
        }

        const bool all_ranks = env_enabled("ABACUS_DIAGO_TRACE_ALL_RANKS");
        if (!all_ranks && GlobalV::MY_RANK != 0)
        {
            this->enabled_ = false;
            return;
        }

        std::string path = "diago_trace.csv";
        const char* filename = std::getenv("ABACUS_DIAGO_TRACE_FILE");
        if (filename != nullptr && filename[0] != '\0')
        {
            path = filename;
        }
        if (all_ranks)
        {
            path = rank_path(path, GlobalV::MY_RANK);
        }

        this->file_.open(path, std::ios::app);
        if (!this->file_)
        {
            this->enabled_ = false;
            return;
        }

        if (this->file_.tellp() == 0)
        {
            this->file_ << "solver,rank,iter,nband,max_residual,avg_residual,n_converged,orth_error,note\n";
        }
    }

    bool enabled() const
    {
        return this->enabled_;
    }

    template <typename Real>
    void record_iteration(const int iter,
                          const int nband,
                          const Real max_residual,
                          const Real avg_residual,
                          const int n_converged,
                          const Real orth_error,
                          const std::string& note = "")
    {
        if (!this->enabled_)
        {
            return;
        }
        this->file_ << this->solver_name_ << ','
                    << GlobalV::MY_RANK << ','
                    << iter << ','
                    << nband << ','
                    << std::setprecision(16) << max_residual << ','
                    << std::setprecision(16) << avg_residual << ','
                    << n_converged << ','
                    << std::setprecision(16) << orth_error << ','
                    << csv_note(note) << '\n';
        this->file_.flush();
    }

  private:
    static bool is_enabled()
    {
        return env_enabled("ABACUS_DIAGO_TRACE");
    }

    static bool env_enabled(const char* name)
    {
        const char* value = std::getenv(name);
        return value != nullptr && value[0] != '\0' && value[0] != '0';
    }

    static std::string csv_note(const std::string& note)
    {
        std::string out = note;
        for (char& ch : out)
        {
            if (ch == ',' || ch == '\n' || ch == '\r')
            {
                ch = ' ';
            }
        }
        return out;
    }

    static std::string rank_path(const std::string& path, const int rank)
    {
        const std::string suffix = ".rank" + std::to_string(rank);
        const std::string::size_type dot = path.find_last_of('.');
        const std::string::size_type slash = path.find_last_of("/\\");
        if (dot != std::string::npos && (slash == std::string::npos || dot > slash))
        {
            return path.substr(0, dot) + suffix + path.substr(dot);
        }
        return path + suffix;
    }

    bool enabled_ = false;
    std::ofstream file_;
    std::string solver_name_;
};

} // namespace hsolver

#endif // DIAGO_TRACE_H_
