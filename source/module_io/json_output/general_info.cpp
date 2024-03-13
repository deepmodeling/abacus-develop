#include "general_info.h"

#include "../para_json.h"
#include "abacusjson.h"
#include "module_base/parallel_global.h"
#include "module_io/input.h"
#include "version.h"

//Add json objects to gener_info
namespace Json
{

#ifdef __RAPIDJSON
void gen_general_info(Input *input)
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
    std::time_t start_time = input->get_start_time();
    std::string start_time_str;
    convert_time(start_time, start_time_str);

    // end_time
    std::time_t time_now = std::time(NULL);
    std::string end_time_str;
    convert_time(time_now, end_time_str);

#ifdef __MPI
    int mpi_num = Parallel_Global::mpi_number;
    int omp_num = Parallel_Global::omp_number;
#elif
    int mpi_num = 1;
    int omp_num = 1;
#endif

    AbacusJson::add_json({"general_info", "version"}, version,false);
    AbacusJson::add_json({"general_info", "commit"}, commit,false);
    AbacusJson::add_json({"general_info", "device"}, input->device,false);
    AbacusJson::add_json({"general_info", "mpi_num"}, mpi_num,false);
    AbacusJson::add_json({"general_info", "omp_num"}, omp_num,false);
    AbacusJson::add_json({"general_info", "pseudo_dir"}, input->pseudo_dir,false);
    AbacusJson::add_json({"general_info", "orbital_dir"}, input->orbital_dir,false);
    AbacusJson::add_json({"general_info", "stru_file"}, input->stru_file,false);
    AbacusJson::add_json({"general_info", "kpt_file"}, input->kpoint_file,false);
    AbacusJson::add_json({"general_info", "start_time"}, start_time_str,false);
    AbacusJson::add_json({"general_info", "end_time"}, end_time_str,false);
}
#endif
} // namespace Json