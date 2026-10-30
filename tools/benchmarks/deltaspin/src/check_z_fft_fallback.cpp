#define main spinor_reference_main
#include "check_spinor_potential.cpp"
#undef main

#include <unordered_set>

namespace
{
const ModulePW::PW_Basis_K* active_basis = nullptr;
bool reject_map = false;
int rejected = 0;
int allocated = 0;
int released = 0;
std::unordered_set<void*> live_maps;
}

extern "C" cudaError_t __real_cudaMalloc(void** pointer, size_t bytes);
extern "C" cudaError_t __real_cudaFree(void* pointer);

extern "C" cudaError_t __wrap_cudaMalloc(void** pointer, size_t bytes)
{
    const bool map_request = active_basis != nullptr && active_basis->nst > 0
        && bytes == size_t(active_basis->nst) * active_basis->nz * active_basis->nks * sizeof(int)
        && ModulePW::fused_fft(active_basis->nz) != nullptr;
    if (map_request && reject_map)
    {
        ++rejected;
        *pointer = nullptr;
        return cudaErrorMemoryAllocation;
    }
    const auto status = __real_cudaMalloc(pointer, bytes);
    if (map_request && status == cudaSuccess)
    {
        live_maps.insert(*pointer);
        ++allocated;
    }
    return status;
}

extern "C" cudaError_t __wrap_cudaFree(void* pointer)
{
    const auto status = __real_cudaFree(pointer);
    if (status == cudaSuccess && live_maps.erase(pointer)) { ++released; }
    return status;
}

void verify_map(const ModulePW::PW_Basis_K& basis)
{
    const size_t size = size_t(basis.nst) * basis.nz;
    std::vector<int> actual(size * basis.nks), expected(actual.size(), -1);
    require(live_maps.size() == 1, "Expected one tracked inverse map");
    check_cuda(cudaMemcpy(actual.data(), *live_maps.begin(), actual.size() * sizeof(int), cudaMemcpyDeviceToHost));
    for (int ik = 0; ik < basis.nks; ++ik)
    for (int ig = 0; ig < basis.npwk[ik]; ++ig)
    {
        const int index = basis.igl2isz_k[ik * basis.npwk_max + ig];
        require(index >= 0 && size_t(index) < size, "Plane-wave map out of bounds");
        auto& value = expected[size_t(ik) * size + index];
        require(value == -1, "Duplicate plane-wave grid position");
        value = ig;
    }
    require(actual == expected, "Inverse map differs from host plane-wave basis");
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    double error = 0, hermitian = 0;
    {
        ModulePW::PW_Basis_K basis("gpu", "double");
        active_basis = &basis;
        run(basis, 36, 36, 36, true, error, hermitian, 50.0, 3);
        require(live_maps.size() == 1, "Eligible map was not allocated");
        verify_map(basis);
        reject_map = true;
        run(basis, 72, 72, 144, true, error, hermitian, 100.0, 19);
        require(rejected == 1 && live_maps.empty(), "Allocation fallback did not release the preceding map");
        reject_map = false;
        run(basis, 72, 72, 31, true, error, hermitian, 50.0, 3);
        require(live_maps.empty(), "Unsupported z length allocated a fused map");
        run(basis, 72, 72, 72, true, error, hermitian, 100.0, 19);
        verify_map(basis);
        run(basis, 36, 36, 36, true, error, hermitian, 50.0, 3);
        verify_map(basis);
        require(live_maps.size() == 1, "Basis rebuild retained multiple maps");
    }
    active_basis = nullptr;
    require(live_maps.empty() && allocated == 3 && released == allocated, "Basis destruction retained a map");
    require(error < 1e-12 && hermitian < 1e-12, "Fallback potential does not match dense reference");
    std::cout << "maps_allocated=" << allocated << " maps_released=" << released << " maps_rejected=" << rejected
              << " error=" << error << " hermitian=" << hermitian << " checks=PASS" << std::endl;
    MPI_Finalize();
}
