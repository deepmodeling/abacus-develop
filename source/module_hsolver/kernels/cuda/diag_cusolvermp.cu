#ifdef __CUSOLVERMP
#include "diag_cusolvermp.cuh"
#include "helper_cuda.h"
#include "module_base/global_variable.h"

#include <assert.h>

extern "C"
{
#include "module_hsolver/genelpa/Cblacs.h"
}
#include <iostream>
#include "helper_cusolver.h"
#include "module_base/global_function.h"
#include "module_base/module_device/device.h"
static calError_t allgather(void* src_buf, void* recv_buf, size_t size, void* data, void** request)
{
    MPI_Request req;
    int err = MPI_Iallgather(src_buf, size, MPI_BYTE, recv_buf, size, MPI_BYTE, (MPI_Comm)(data), &req);
    if (err != MPI_SUCCESS)
    {
        return CAL_ERROR;
    }
    *request = (void*)(req);
    return CAL_OK;
}

static calError_t request_test(void* request)
{
    MPI_Request req = (MPI_Request)(request);
    int completed;
    int err = MPI_Test(&req, &completed, MPI_STATUS_IGNORE);
    if (err != MPI_SUCCESS)
    {
        return CAL_ERROR;
    }
    return completed ? CAL_OK : CAL_ERROR_INPROGRESS;
}

static calError_t request_free(void* request)
{
    return CAL_OK;
}

template <typename inputT>
Diag_CusolverMP_gvd<inputT>::Diag_CusolverMP_gvd(const MPI_Comm mpi_comm,
                                                 const int narows,
                                                 const int nacols,
                                                 const int* desc)
{
    // 构造函数的实现
    this->cblacs_ctxt = desc[1];
    this->nFull = desc[2];

    // 20240529 zhanghaochong
    // set mb and nb is not nessary, but I keep it here for the code consistency.
    // Because in ABACUS, mb always equals to nb
    const int mb = desc[4];
    const int nb = desc[5];

    // 20240529 zhanghaochong
    // so far, cusolverMpSygvd only support rsrc == 0 and csrc == 0
    const int rsrc = desc[6];
    const int csrc = desc[7];

    this->lda = desc[8];

    MPI_Comm_size(mpi_comm, &this->globalMpiSize);
    MPI_Comm_rank(mpi_comm, &(this->globalMpiRank));
    int local_device_id = base_device::information::set_device_by_rank(mpi_comm);
    Cblacs_gridinfo(this->cblacs_ctxt, &this->nprows, &this->npcols, &this->myprow, &this->mypcol);

    this->cusolverCalComm = NULL;
    cal_comm_create_params_t params;
    params.allgather = allgather;
    params.req_test = request_test;
    params.req_free = request_free;
    params.data = (void*)(mpi_comm);
    params.rank = this->globalMpiRank;
    params.nranks = this->globalMpiSize;
    params.local_device = local_device_id;

    CAL_CHECK(cal_comm_create(params, &this->cusolverCalComm));

    checkCudaErrors(cudaStreamCreate(&this->localStream));
    CUSOLVER_CHECK(cusolverMpCreate(&cusolverMpHandle, local_device_id, this->localStream));

    // 20240529 zhanghaochong
    // so far, cusolvermp only support = 1
    this->matrix_i = 1;
    this->matrix_j = 1;
    this->m_local = narows;
    this->n_local = nacols;
    if (std::is_same<inputT, std::complex<double>>::value)
    {
        this->datatype = CUDA_C_64F;
    }
    else if (std::is_same<inputT, double>::value)
    {
        this->datatype = CUDA_R_64F;
    }
    else
    {
        ModuleBase::WARNING_QUIT("cusolvermp", "unsupported cusolvermp datatype");
    }

    // 20240529 zhanghaochong
    // the cpu mpi process blacs grid and multi gpu process blacs grid is the SAME
    // Setting them the same is not a natural result, but a result that I forced and artificially specified.
    // This is because the current implementation of the cusolvermp library is ONE process ONE GPU.
    // So, when we use cusolvermp, we must ensure that the number of processes is equal to the number of GPUs.
    // In a sense, the MPI usage strategy of ABACUS must be subject to the cusolvermp.
    CUSOLVER_CHECK(cusolverMpCreateDeviceGrid(cusolverMpHandle,
                                                   &this->grid,
                                                   this->cusolverCalComm,
                                                   this->nprows,
                                                   this->npcols,
                                                   CUSOLVERMP_GRID_MAPPING_COL_MAJOR));

    // 20240529 zhanghaochong
    // Actually, there should be three matrix descriptors, A matrix, B matrix, and output eigenvector matrix.
    // But in ABACUS the three matrices descriptors are the same.
    // So, I only create one matrix descriptor and use it for the three matrices.
    CUSOLVER_CHECK(cusolverMpCreateMatrixDesc(&this->desc_for_cusolvermp,
                               this->grid,
                               this->datatype,
                               nFull,
                               nFull,
                               mb,
                               nb,
                               rsrc,
                               csrc,
                               this->lda));
}

template <typename inputT>
Diag_CusolverMP_gvd<inputT>::~Diag_CusolverMP_gvd()
{
    CAL_CHECK(cal_comm_barrier(this->cusolverCalComm, this->localStream));
    CUSOLVER_CHECK(cusolverMpDestroyMatrixDesc(this->desc_for_cusolvermp));
    CUSOLVER_CHECK(cusolverMpDestroyGrid(this->grid));
    CUSOLVER_CHECK(cusolverMpDestroy(this->cusolverMpHandle));
    CAL_CHECK(cal_comm_destroy(this->cusolverCalComm));
    checkCudaErrors(cudaStreamDestroy(this->localStream));
}


template <typename inputT>
int Diag_CusolverMP_gvd<inputT>::generalized_eigenvector(inputT* A, inputT* B, outputT* EigenValue, inputT* EigenVector)
{
    void *d_A = NULL;
    void *d_B = NULL;
    void *d_D = NULL;
    void *d_Z = NULL;

    checkCudaErrors(cudaMalloc((void**)&d_A, this->n_local * this->m_local * sizeof(inputT)));
    checkCudaErrors(cudaMalloc((void**)&d_B, this->n_local * this->m_local * sizeof(inputT)));
    checkCudaErrors(cudaMalloc((void**)&d_D, this->nFull * sizeof(outputT)));
    checkCudaErrors(cudaMalloc((void**)&d_Z, this->n_local * this->m_local * sizeof(inputT)));

    checkCudaErrors(
        cudaMemcpy(d_A, (void*)A, this->n_local * this->m_local * sizeof(inputT), cudaMemcpyHostToDevice));
    checkCudaErrors(
        cudaMemcpy(d_B, (void*)B, this->n_local * this->m_local * sizeof(inputT), cudaMemcpyHostToDevice));
    CAL_CHECK(cal_stream_sync(this->cusolverCalComm, this->localStream));

    size_t sygvdWorkspaceInBytesOnDevice = 0;
    size_t sygvdWorkspaceInBytesOnHost = 0;

    CUSOLVER_CHECK(cusolverMpSygvd_bufferSize(cusolverMpHandle,
                               CUSOLVER_EIG_TYPE_1,
                               CUSOLVER_EIG_MODE_VECTOR,
                               CUBLAS_FILL_MODE_LOWER,
                               this->nFull,
                               this->matrix_i,
                               this->matrix_j,
                               this->desc_for_cusolvermp,
                               this->matrix_i,
                               this->matrix_j,
                               this->desc_for_cusolvermp,
                               this->matrix_i,
                               this->matrix_j,
                               this->desc_for_cusolvermp,
                               this->datatype,
                               &sygvdWorkspaceInBytesOnDevice,
                               &sygvdWorkspaceInBytesOnHost));

    /* Distributed host workspace */
    void* h_sygvdWork = NULL;
    h_sygvdWork = (void*)malloc(sygvdWorkspaceInBytesOnHost);

    /* Distributed device workspace */
    void* d_sygvdWork = NULL;
    checkCudaErrors(cudaMalloc((void**)&d_sygvdWork, sygvdWorkspaceInBytesOnDevice));
    checkCudaErrors(cudaMemset(d_sygvdWork, 0, sygvdWorkspaceInBytesOnDevice));

    int* d_sygvdInfo = NULL;
    checkCudaErrors(cudaMalloc((void**)&d_sygvdInfo, sizeof(int)));
    checkCudaErrors(cudaMemset(d_sygvdInfo, 0, sizeof(int)));

    /* sync wait for data to arrive to device */
    CAL_CHECK(cal_stream_sync(cusolverCalComm, localStream));

    CUSOLVER_CHECK(cusolverMpSygvd(cusolverMpHandle,
                    CUSOLVER_EIG_TYPE_1,
                    CUSOLVER_EIG_MODE_VECTOR,
                    CUBLAS_FILL_MODE_LOWER,
                    this->nFull,
                    d_A,
                    this->matrix_i,
                    this->matrix_j,
                    this->desc_for_cusolvermp,
                    d_B,
                    this->matrix_i,
                    this->matrix_j,
                    this->desc_for_cusolvermp,
                    d_D,
                    d_Z,
                    this->matrix_i,
                    this->matrix_j,
                    this->desc_for_cusolvermp,
                    this->datatype,
                    d_sygvdWork,
                    sygvdWorkspaceInBytesOnDevice,
                    h_sygvdWork,
                    sygvdWorkspaceInBytesOnHost,
                    d_sygvdInfo));

    int h_sygvdInfo = 0;
    checkCudaErrors(cudaMemcpyAsync(&h_sygvdInfo, d_sygvdInfo, sizeof(int), cudaMemcpyDeviceToHost, this->localStream));
    /* wait for d_sygvdInfo copy */
    checkCudaErrors(cudaStreamSynchronize(this->localStream));
    if (h_sygvdInfo != 0)
    {
        ModuleBase::WARNING_QUIT("cusolvermp", "cusolverMpSygvd failed with error");
    }
    CAL_CHECK(cal_stream_sync(this->cusolverCalComm, this->localStream));

    checkCudaErrors(cudaFree(d_sygvdWork));
    checkCudaErrors(cudaFree(d_sygvdInfo));

    free(h_sygvdWork);

    checkCudaErrors(cudaMemcpy((void*)EigenValue, d_D, this->nFull * sizeof(outputT), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy((void*)EigenVector,
                               d_Z,
                               this->n_local * this->m_local * sizeof(inputT),
                               cudaMemcpyDeviceToHost));
    // 20240529 zhanghaochong
    // I move the free operations from destructor to here.
    // Because I think it is more reasonable to free the memory in the function where it is allocated.
    // Destructor is used to release resources that allocated in the constructor.
    // And currently, we construct and destruct the object in every SCF iteration. Maybe one day we 
    // will construct the object only once during the whole program life cycle.
    // In that case, allocate and free memory in compute function is more reasonable.
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_D));
    checkCudaErrors(cudaFree(d_Z));

    return 0;
}
template <typename inputT>
void Diag_CusolverMP_gvd<inputT>::outputParameters()
{
    GlobalV::ofs_running << "nFull: " << this->nFull << std::endl
                         << "m_local: " << this->m_local << std::endl
                         << "n_local: " << this->n_local << std::endl
                         << "lda: " << this->lda << std::endl
                         << "nprows: " << this->nprows << std::endl
                         << "npcols: " << this->npcols << std::endl
                         << "myprow: " << this->myprow << std::endl
                         << "mypcol: " << this->mypcol << std::endl
                         << "globalMpiRank: " << this->globalMpiRank << std::endl
                         << "globalMpiSize: " << this->globalMpiSize << std::endl
                         << "matrix_i: " << this->matrix_i << std::endl
                         << "matrix_j: " << this->matrix_j << std::endl;
}

template class Diag_CusolverMP_gvd<double>;
template class Diag_CusolverMP_gvd<std::complex<double>>;
#endif