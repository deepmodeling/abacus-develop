#include "module_io/to_qo.h"
#ifdef __MPI
#include "../module_base/parallel_common.h"
#endif

// free function, to broadcast supercells_ and kvec_d_
void toQO::bcast_stdvector_ofvector3int(std::vector<ModuleBase::Vector3<int>>& vec)
{
    #ifdef __MPI
    int dim;
    std::vector<int> vec_1d;
    if(iproc_ == 0)
    {
        dim = vec.size();
        for(int i = 0; i < dim; i++)
        {
            vec_1d.push_back(vec[i].x);
            vec_1d.push_back(vec[i].y);
            vec_1d.push_back(vec[i].z);
        }
    }
    Parallel_Common::bcast_int(dim);
    if(iproc_ != 0) vec_1d.resize(dim * 3);
    Parallel_Common::bcast_int(vec_1d.data(), dim * 3);
    if(iproc_ != 0)
    {
        vec.clear(); vec.resize(dim);
        for(int i = 0; i < dim; i++)
        {
            vec.push_back(ModuleBase::Vector3<int>(vec_1d[i*3], vec_1d[i*3+1], vec_1d[i*3+2]));
        }
    }
    #endif
}

void toQO::bcast_stdvector_ofvector3double(std::vector<ModuleBase::Vector3<double>>& vec)
{
    #ifdef __MPI
    int dim;
    std::vector<double> vec_1d;
    if(iproc_ == 0)
    {
        dim = vec.size();
        for(int i = 0; i < dim; i++)
        {
            vec_1d.push_back(vec[i].x);
            vec_1d.push_back(vec[i].y);
            vec_1d.push_back(vec[i].z);
        }
    }
    Parallel_Common::bcast_int(dim);
    if(iproc_ != 0) vec_1d.resize(dim * 3);
    Parallel_Common::bcast_double(vec_1d.data(), dim * 3);
    if(iproc_ != 0)
    {
        vec.clear(); vec.resize(dim);
        for(int i = 0; i < dim; i++)
        {
            vec.push_back(ModuleBase::Vector3<double>(vec_1d[i*3], vec_1d[i*3+1], vec_1d[i*3+2]));
        }
    }
    #endif
}
