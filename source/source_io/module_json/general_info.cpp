#include "general_info.h"

#include "para_json.h"
#include "abacusjson.h"
#include "source_io/module_parameter/parameter.h"

#ifdef __JSON
#include <nlohmann/json.hpp>
#endif
#include "source_base/parallel_global.h"
#include "source_main/version.h"

namespace Json
{

#ifdef __JSON
void gen_general_info(const Parameter& param)
{

#ifdef VERSION
    const std::string version = VERSION;
#else
    const std::string version = "unknown";
#endif
#ifdef COMMIT_INFO
#include "commit.h"
    const std::string commit = COMMIT;
#else
    const std::string commit = "unknown";
#endif

    // start_time
    std::time_t start_time = param.globalv.start_time;
    std::string start_time_str;
    convert_time(start_time, start_time_str);

    // end_time
    std::time_t time_now = std::time(nullptr);
    std::string end_time_str;
    convert_time(time_now, end_time_str);

#ifdef __MPI
    int mpi_num = Parallel_Global::mpi_number;
    int omp_num = Parallel_Global::omp_number;
#else
    int mpi_num = 1;
    int omp_num = 1;
#endif

    AbacusJson::document()["general_info"] = {
        {"version", version},
        {"commit", commit},
        {"device", param.inp.device},
        {"mpi_num", mpi_num},
        {"omp_num", omp_num},
        {"pseudo_dir", param.inp.pseudo_dir},
        {"orbital_dir", param.inp.orbital_dir},
        {"stru_file", param.globalv.global_in_stru},
        {"kpt_file", param.inp.kpoint_file},
        {"start_time", start_time_str},
        {"end_time", end_time_str}};
}
#endif
} // namespace Json
