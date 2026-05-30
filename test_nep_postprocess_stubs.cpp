// Stub for BlasConnector linking (not needed for our unit tests)
namespace base_device { enum class AbacusDevice_t {CpuDevice, GpuDevice}; }

struct BlasConnector {
    static void gemm(char, char, int, int, int, double, const double*, int,
                     const double*, int, double, double*, int,
                     base_device::AbacusDevice_t) {}
    static double nrm2(int, const double*, int, base_device::AbacusDevice_t) { return 0.0; }
};
