#define main fixed_overlap_reference_main
#include "check_fixed_overlap.cpp"
#undef main
#include <map>
#include <cusolver_common.h>

namespace
{
int fail_after = -1;
bool track_allocations = false;
std::map<void*, std::size_t> tracked_allocations;
}

extern "C" cudaError_t __real_cudaMalloc(void** pointer, std::size_t bytes);
extern "C" cudaError_t __real_cudaFree(void* pointer);

extern "C" cudaError_t __wrap_cudaMalloc(void** pointer, std::size_t bytes)
{
    if (fail_after == 0)
    {
        fail_after = -1;
        *pointer = nullptr;
        return cudaErrorMemoryAllocation;
    }
    if (fail_after > 0) { --fail_after; }
    const auto status = __real_cudaMalloc(pointer, bytes);
    if (track_allocations && status == cudaSuccess) { tracked_allocations[*pointer] = bytes; }
    return status;
}

extern "C" cudaError_t __wrap_cudaFree(void* pointer)
{
    const auto status = __real_cudaFree(pointer);
    if (status == cudaSuccess) { tracked_allocations.erase(pointer); }
    return status;
}

void check_allocation_failure(int fail_at)
{
    const int bands = 7, points = 3;
    std::vector<Complex> h(bands * bands), s(bands * bands);
    for (int band = 0; band < bands; ++band)
    {
        h[band * bands + band] = 0.1 * band;
        s[band * bands + band] = 1;
    }
    DeviceBuffer input_h(h.size()), input_s(s.size()), output(s.size());
    input_h.set(h); input_s.set(s);
    hsolver::FixedOverlapSubspace cache(points, bands, bands);
    require(cache.available(), "Fallback fixture could not allocate");
    for (int point = 0; point < points; ++point)
    { cache.prepare(point, input_h.data(), input_s.data(), input_s.data()); }
    require(tracked_allocations.empty(), "Prior batch allocation remains live");
    track_allocations = true;
    fail_after = fail_at;
    const bool ready = cache.prepare_batch();
    track_allocations = false;
    require(fail_after == -1 && !ready, "Optional allocation failure did not select serial fallback");
    require(tracked_allocations.empty(), "Failed batch preparation leaked device storage");
    require(!cache.prepare_batch(), "Unavailable batch unexpectedly retried allocation");
    std::vector<double> values(bands);
    for (int point = points - 1; point >= 0; --point)
    {
        cache.diagonalize(input_h.data(), cache.projectors(point), output.data(), values.data());
        for (int band = 0; band < bands; ++band)
        { require(std::abs(values[band] - 0.1 * band) < 1e-13, "Serial fallback spectrum changed"); }
    }
    require(input_h.get() == h && input_s.get() == s, "Fallback changed an input");
    std::cout << "optional_allocation_failure=" << fail_at << " fallback=pass live_bytes=0" << std::endl;
}

void run_batch(int n, int points)
{
    const int rows = n + 2, value_leading = n + 5;
    const std::size_t matrix = std::size_t(n) * n;
    const long double pi = std::acos(-1.0L);
    std::vector<std::vector<Complex>> h(points, std::vector<Complex>(matrix));
    auto s = h;
    std::vector<Complex> projectors(std::size_t(rows) * n);
    for (int column = 0; column < n; ++column)
    {
        projectors[std::size_t(column) * rows + column] = 1;
        for (int row = n; row < rows; ++row)
        {
            projectors[std::size_t(column) * rows + row]
                = std::polar(0.1 / std::sqrt(n), 0.13 * (column + 1) * (row - n + 1));
        }
    }
    for (int point = 0; point < points; ++point)
    for (int column = 0; column < n; ++column)
    for (int row = 0; row < n; ++row)
    {
        Wide hs = 0, ss = 0;
        for (int state = 0; state < n; ++state)
        {
            const long double overlap = std::exp(-std::log(1000.0L) * (state + 0.03L * point) / n);
            const long double energy = 0.3L * point + 0.01L * state - 0.4L;
            const Wide phase = std::polar(1.0L / n,
                2 * pi * (row - column) * state / n + 0.17L * point * (row - column));
            hs += phase * overlap * energy;
            ss += phase * overlap;
        }
        h[point][std::size_t(column) * n + row] = Complex(hs);
        s[point][std::size_t(column) * n + row] = Complex(ss);
    }
    DeviceBuffer input_h(matrix), input_s(matrix), input_r(projectors.size());
    DeviceBuffer shifted(matrix), rotated(projectors.size());
    input_r.set(projectors);
    hsolver::FixedOverlapSubspace cache(points, n, rows);
    require(cache.available(), "Fixed response fixture could not allocate");
    std::vector<std::vector<Complex>> cached_h(points, std::vector<Complex>(matrix));
    std::vector<std::vector<Complex>> cached_r(points, std::vector<Complex>(projectors.size()));
    for (int point = 0; point < points; ++point)
    {
        input_h.set(h[point]); input_s.set(s[point]);
        cache.prepare(point, input_h.data(), input_s.data(), input_r.data());
        require(input_h.get() == h[point] && input_s.get() == s[point] && input_r.get() == projectors,
                "Batch preparation changed an input");
        check(cudaMemcpy(cached_h[point].data(), cache.hamiltonian(point), matrix * sizeof(Complex), cudaMemcpyDeviceToHost));
        check(cudaMemcpy(cached_r[point].data(), cache.projectors(point), projectors.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
    }
    track_allocations = true;
    const bool batched = cache.prepare_batch();
    track_allocations = false;
    require(batched == (points > 1 && CUSOLVER_VERSION >= 11703), "Unexpected batch eligibility");
    double maximum_residual = 0, maximum_orthogonality = 0, maximum_spectrum = 0, maximum_projection = 0;
    for (double field : {-0.017, 0.0, 0.023})
    {
        const double guard = 91.25;
        std::vector<double> values(std::size_t(points) * value_leading + 4, guard);
        for (int point = points - 1; point >= 0; --point)
        {
            auto candidate_h = cached_h[point];
            for (int row = 0; row < n; ++row) { candidate_h[std::size_t(row) * n + row] += field; }
            shifted.set(candidate_h);
            if (batched) { cache.capture_batch(point, shifted.data()); }
            require(shifted.get() == candidate_h, "Response capture changed an input");
        }
        if (batched) { cache.diagonalize_batch(values.data() + 2, value_leading); }
        require(values[0] == guard && values[1] == guard && values[values.size() - 1] == guard
                && values[values.size() - 2] == guard, "Host spectrum guard changed");
        for (int point = 0; point < points; ++point)
        {
            if (batched) { cache.rotate_batch(point, rotated.data()); }
            else
            {
                auto candidate_h = cached_h[point];
                for (int row = 0; row < n; ++row) { candidate_h[std::size_t(row) * n + row] += field; }
                shifted.set(candidate_h);
                cache.diagonalize(shifted.data(), cache.projectors(point), rotated.data(),
                                  values.data() + std::size_t(point) * value_leading + 2);
                require(shifted.get() == candidate_h, "Serial response changed an input");
            }
            const auto projections = rotated.get();
            std::vector<Complex> vectors(matrix), physical_h(h[point]);
            for (std::size_t index = 0; index < matrix; ++index) { physical_h[index] += field * s[point][index]; }
            for (int column = 0; column < n; ++column)
            {
                std::copy_n(projections.data() + std::size_t(column) * rows, n,
                            vectors.data() + std::size_t(column) * n);
                for (int row = n; row < rows; ++row)
                {
                    Wide expected = 0;
                    for (int inner = 0; inner < n; ++inner)
                    {
                        expected += Wide(projectors[std::size_t(inner) * rows + row])
                                    * Wide(vectors[std::size_t(column) * n + inner]);
                    }
                    const double error = double(std::abs(expected - Wide(projections[std::size_t(column) * rows + row])));
                    require(std::isfinite(error), "Nonfinite rotated projector");
                    maximum_projection = std::max(maximum_projection, error);
                }
            }
            const auto hv = multiply(physical_h, vectors, n), sv = multiply(s[point], vectors, n);
            long double norm_h = 0, norm_s = 0;
            for (std::size_t index = 0; index < matrix; ++index)
            {
                norm_h += std::norm(Wide(physical_h[index]));
                norm_s += std::norm(Wide(s[point][index]));
            }
            norm_h = std::sqrt(norm_h); norm_s = std::sqrt(norm_s);
            for (int column = 0; column < n; ++column)
            {
                const double value = values[std::size_t(point) * value_leading + column + 2];
                maximum_spectrum = std::max(maximum_spectrum,
                    std::abs(value - (0.3 * point + 0.01 * column - 0.4 + field)));
                long double residual = 0, vector_norm = 0;
                for (int row = 0; row < n; ++row)
                {
                    const std::size_t index = std::size_t(column) * n + row;
                    residual += std::norm(Wide(hv[index]) - static_cast<long double>(value) * Wide(sv[index]));
                    vector_norm += std::norm(Wide(vectors[index]));
                    Wide overlap = 0;
                    for (int inner = 0; inner < n; ++inner)
                    {
                        overlap += std::conj(Wide(vectors[std::size_t(row) * n + inner]))
                                   * Wide(sv[std::size_t(column) * n + inner]);
                    }
                    maximum_orthogonality = std::max(maximum_orthogonality,
                        double(std::abs(overlap - Wide(row == column))));
                }
                const double error = double(std::sqrt(residual / vector_norm) / (norm_h + std::abs(value) * norm_s));
                require(std::isfinite(value) && std::isfinite(error), "Nonfinite batch eigenpair");
                maximum_residual = std::max(maximum_residual, error);
            }
            for (int column = n; column < value_leading; ++column)
            { require(values[std::size_t(point) * value_leading + column + 2] == guard, "Spectrum padding changed"); }
            std::vector<Complex> preserved(matrix), preserved_r(projectors.size());
            check(cudaMemcpy(preserved.data(), cache.hamiltonian(point), matrix * sizeof(Complex), cudaMemcpyDeviceToHost));
            check(cudaMemcpy(preserved_r.data(), cache.projectors(point), projectors.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
            require(preserved == cached_h[point] && preserved_r == cached_r[point], "Cached basis products changed");
        }
        require(values[0] == guard && values[1] == guard && values[values.size() - 1] == guard
                && values[values.size() - 2] == guard, "Completed spectrum guard changed");
    }
    std::cout << "bands=" << n << " points=" << points << " batched=" << batched
              << " backward_error=" << maximum_residual << " orthogonality=" << maximum_orthogonality
              << " spectrum=" << maximum_spectrum << " projection=" << maximum_projection << std::endl;
    require(maximum_residual < 1e-12 && maximum_orthogonality < 1e-10 && maximum_spectrum < 1e-10
            && maximum_projection < 1e-10, "Batch response checks failed");
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    hsolver::createGpuBlasHandle();
    hsolver::createGpuSolverHandle();
    for (int bands : {1, 7, 31, 65})
    for (int points : {1, 3, 8})
    {
        run_batch(bands, points);
        require(tracked_allocations.empty(), "Batch destruction leaked device storage");
    }
    run_batch(257, 3);
    require(tracked_allocations.empty(), "Large batch destruction leaked device storage");
    if (CUSOLVER_VERSION >= 11703)
    {
        check_allocation_failure(0);
        check_allocation_failure(1);
    }
    hsolver::destroyGpuSolverHandle();
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
