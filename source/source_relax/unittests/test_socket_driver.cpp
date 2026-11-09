#include "for_test.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "source_base/constants.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/parameter.h"
#include "source_relax/socket_constant_potential.h"
#include "source_relax/socket_driver.h"
#include "source_relax/socket_handlers.h"
#include "source_relax/socket_utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <unistd.h>

int XC_Functional::func_type = 2;

class TestParameters
{
  public:
    static Input_para& input(Parameter& parameters) { return parameters.input; }
};

namespace
{
constexpr std::size_t kIpiHeaderLen = 12;

std::string errno_message(const std::string& prefix)
{
    return prefix + ": " + std::strerror(errno);
}

void send_all(const int fd, const void* data, const std::size_t nbytes)
{
    const char* cursor = static_cast<const char*>(data);
    std::size_t done = 0;
    while (done < nbytes)
    {
#ifdef MSG_NOSIGNAL
        const int flags = MSG_NOSIGNAL;
#else
        const int flags = 0;
#endif
        const ssize_t sent = ::send(fd, cursor + done, nbytes - done, flags);
        if (sent < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("send failed"));
        }
        if (sent == 0)
        {
            throw std::runtime_error("send returned zero");
        }
        done += static_cast<std::size_t>(sent);
    }
}

template <typename T>
void send_value(const int fd, const T& value)
{
    send_all(fd, &value, sizeof(value));
}

void send_header(const int fd, const std::string& header)
{
    std::string padded = header;
    padded.resize(kIpiHeaderLen, ' ');
    send_all(fd, padded.data(), padded.size());
}

bool try_send_status(const int fd)
{
    try
    {
        send_header(fd, "STATUS");
        return true;
    }
    catch (const std::runtime_error&)
    {
        if (errno == EPIPE || errno == ECONNRESET)
        {
            return false;
        }
        throw;
    }
}

std::string read_header_or_close(const int fd)
{
    char header[kIpiHeaderLen];
    std::size_t done = 0;
    while (done < sizeof(header))
    {
        const ssize_t received = ::recv(fd, header + done, sizeof(header) - done, 0);
        if (received == 0 || (received < 0 && errno == ECONNRESET))
        {
            if (done == 0)
            {
                return "";
            }
            throw std::runtime_error("socket closed during response header");
        }
        if (received < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("receive failed"));
        }
        done += static_cast<std::size_t>(received);
    }

    std::string value(header, sizeof(header));
    while (!value.empty() && value.back() == ' ')
    {
        value.pop_back();
    }
    return value;
}

class UnixSocketServer
{
  public:
    UnixSocketServer()
    {
        char dir_template[] = "/tmp/abacus_socket_driver_test_XXXXXX";
        char* made_dir = ::mkdtemp(dir_template);
        if (made_dir == nullptr)
        {
            throw std::runtime_error(errno_message("mkdtemp failed"));
        }
        dir_ = made_dir;
        path_ = dir_ + "/ipi.sock";

        listen_fd_ = ::socket(AF_UNIX, SOCK_STREAM, 0);
        if (listen_fd_ < 0)
        {
            throw std::runtime_error(errno_message("socket failed"));
        }

        sockaddr_un address;
        std::memset(&address, 0, sizeof(address));
        address.sun_family = AF_UNIX;
        std::strncpy(address.sun_path, path_.c_str(), sizeof(address.sun_path) - 1);
        if (::bind(listen_fd_, reinterpret_cast<sockaddr*>(&address), sizeof(address)) != 0)
        {
            throw std::runtime_error(errno_message("bind failed"));
        }
        if (::listen(listen_fd_, 1) != 0)
        {
            throw std::runtime_error(errno_message("listen failed"));
        }
    }

    ~UnixSocketServer()
    {
        if (listen_fd_ >= 0)
        {
            ::close(listen_fd_);
        }
        if (!path_.empty())
        {
            ::unlink(path_.c_str());
        }
        if (!dir_.empty())
        {
            ::rmdir(dir_.c_str());
        }
    }

    UnixSocketServer(const UnixSocketServer&) = delete;
    UnixSocketServer& operator=(const UnixSocketServer&) = delete;

    std::string address() const
    {
        return path_ + ":UNIX";
    }

    int accept_once() const
    {
        const int fd = ::accept(listen_fd_, nullptr, nullptr);
        if (fd < 0)
        {
            throw std::runtime_error(errno_message("accept failed"));
        }
        return fd;
    }

  private:
    int listen_fd_ = -1;
    std::string dir_;
    std::string path_;
};

class FakeESolver : public ModuleESolver::ESolver
{
  public:
    explicit FakeESolver(const bool converged) : converged_(converged)
    {
    }

    void before_all_runners(BaseCell&, const Input_para&) override
    {
    }

    void runner(BaseCell& cell, const int step) override
    {
        UnitCell& ucell = dynamic_cast<UnitCell&>(cell);
        electronic_context_seen = ucell.electronic_step != nullptr;
        position_ = ucell.atoms[0].tau[0].x;
        this->conv_esolver = converged_ && step == 0;
    }

    void after_all_runners(BaseCell&) override
    {
    }

    double cal_energy() override
    {
        return 4.0 + position_;
    }

    void cal_force(BaseCell& cell, ModuleBase::matrix& force) override
    {
        force.create(dynamic_cast<UnitCell&>(cell).nat, 3);
        force(0, 0) = 4.0 + 2.0 * position_;
    }

    void cal_stress(BaseCell&, ModuleBase::matrix& stress) override
    {
        stress.create(3, 3);
        stress(0, 0) = 2.0 + 3.0 * position_;
        stress(1, 1) = 2.0;
        stress(2, 2) = 2.0;
    }

    bool electronic_context_seen = false;

  private:
    double position_ = 0.0;
    bool converged_;
};

enum class ElectronicFailure
{
    None,
    SecondFrame,
    IgnoreRequest,
    MissingResult
};

class ControlledESolver : public FakeESolver
{
  public:
    explicit ControlledESolver(const ElectronicFailure failure) : FakeESolver(true), failure_(failure) {}

    void runner(BaseCell& cell, const int step) override
    {
        UnitCell& ucell = dynamic_cast<UnitCell&>(cell);
        if (ucell.electronic_step == nullptr || ucell.electronic_step->request_applied
            || ucell.electronic_step->result_ready)
        {
            throw std::runtime_error("each runner needs one fresh electronic request");
        }
        if (failure_ != ElectronicFailure::IgnoreRequest)
        {
            nelect = ucell.electronic_step->control.nelect;
            ucell.electronic_step->request_applied = true;
            ++control_count;
        }
        FakeESolver::runner(cell, step);
        ++runner_count;
        conv_esolver = failure_ != ElectronicFailure::SecondFrame || step == 0;
    }

    double cal_energy() override { return FakeESolver::cal_energy() + nelect; }

    void cal_force(BaseCell& cell, ModuleBase::matrix& force) override
    {
        FakeESolver::cal_force(cell, force);
        elecstate::ElectronicStep* electronic = dynamic_cast<UnitCell&>(cell).electronic_step;
        if (!electronic || !electronic->request_applied || electronic->result_ready || !conv_esolver)
        {
            throw std::runtime_error("electronic result requires one converged request");
        }
        if (failure_ == ElectronicFailure::MissingResult) { return; }
        electronic->result.nelect = nelect;
        electronic->result.fermi_ry = -0.125;
        electronic->result_ready = true;
        ++collection_count;
    }

    double nelect = 0.0;
    int control_count = 0;
    int runner_count = 0;
    int collection_count = 0;

  private:
    ElectronicFailure failure_;
};

struct DriverResult
{
    int exit_code = -1;
    std::string response_header;
    std::string diagnostic;
};

struct ForceResponse
{
    std::string header;
    double energy_hartree = 0.0;
    std::int32_t nat = 0;
    std::vector<double> forces_hartree_per_bohr;
    std::vector<double> virial_wire_hartree;
    std::string extra;
};

void initialize_one_atom_cell(UnitCell& ucell)
{
    ucell.lat0 = 1.0;
    ucell.latvec.Identity();
    ucell.omega = 1.0;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms[0].na = 1;
    ucell.atoms[0].tau.resize(1);
    ucell.atoms[0].taud.resize(1);
    ucell.atoms[0].dis.resize(1);
}

void send_positions(const int fd, const double x)
{
    const double identity[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    const std::int32_t nat = 1;
    const double position[3] = {x, 0.0, 0.0};
    send_header(fd, "POSDATA");
    send_all(fd, identity, sizeof(identity));
    send_all(fd, identity, sizeof(identity));
    send_value(fd, nat);
    send_all(fd, position, sizeof(position));
}

void send_fixed_cell_frame(const int fd)
{
    const std::int32_t replica = 0;
    const std::int32_t parameter_bytes = 0;
    send_header(fd, "INIT");
    send_value(fd, replica);
    send_value(fd, parameter_bytes);

    send_positions(fd, 0.0);
}

void read_all(const int fd, void* data, const std::size_t nbytes)
{
    char* cursor = static_cast<char*>(data);
    std::size_t done = 0;
    while (done < nbytes)
    {
        const ssize_t received = ::recv(fd, cursor + done, nbytes - done, 0);
        if (received <= 0)
        {
            throw std::runtime_error("socket closed while reading response");
        }
        done += static_cast<std::size_t>(received);
    }
}

template <typename T>
T read_value(const int fd)
{
    T value;
    read_all(fd, &value, sizeof(value));
    return value;
}

std::vector<double> read_doubles(const int fd, const std::size_t count)
{
    std::vector<double> values(count);
    if (!values.empty())
    {
        read_all(fd, values.data(), values.size() * sizeof(double));
    }
    return values;
}

ForceResponse read_force_response(const int fd)
{
    ForceResponse response;
    response.header = read_header_or_close(fd);
    if (response.header.empty())
    {
        return response;
    }
    response.energy_hartree = read_value<double>(fd);
    response.nat = read_value<std::int32_t>(fd);
    response.forces_hartree_per_bohr = read_doubles(fd, static_cast<std::size_t>(3 * response.nat));
    response.virial_wire_hartree = read_doubles(fd, 9);
    const std::int32_t extra_bytes = read_value<std::int32_t>(fd);
    if (extra_bytes < 0)
    {
        throw std::runtime_error("negative extras length");
    }
    response.extra.resize(static_cast<std::size_t>(extra_bytes));
    if (!response.extra.empty())
    {
        read_all(fd, &response.extra[0], response.extra.size());
    }
    return response;
}

std::string read_pipe(const int fd)
{
    std::string output;
    char buffer[512];
    while (true)
    {
        const ssize_t nread = ::read(fd, buffer, sizeof(buffer));
        if (nread == 0)
        {
            break;
        }
        if (nread < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("pipe read failed"));
        }
        output.append(buffer, static_cast<std::size_t>(nread));
    }
    return output;
}

DriverResult run_driver_session(const std::function<void()>& driver_action,
                                const std::function<void(int)>& peer_action)
{
    UnixSocketServer server;
    int output_pipe[2];
    if (::pipe(output_pipe) != 0)
    {
        throw std::runtime_error(errno_message("pipe failed"));
    }

    const pid_t child = ::fork();
    if (child < 0)
    {
        ::close(output_pipe[0]);
        ::close(output_pipe[1]);
        throw std::runtime_error(errno_message("fork failed"));
    }
    if (child == 0)
    {
        ::close(output_pipe[0]);
        ::dup2(output_pipe[1], STDOUT_FILENO);
        ::dup2(output_pipe[1], STDERR_FILENO);
        ::close(output_pipe[1]);
        ::setenv("ABACUS_SOCKET_ADDRESS", server.address().c_str(), 1);

        driver_action();
        std::cout.flush();
        std::cerr.flush();
        ::_exit(0);
    }

    ::close(output_pipe[1]);
    DriverResult result;
    std::exception_ptr peer_error;
    int peer_fd = -1;
    try
    {
        peer_fd = server.accept_once();
        timeval timeout;
        timeout.tv_sec = 5;
        timeout.tv_usec = 0;
        if (::setsockopt(peer_fd, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) != 0)
        {
            throw std::runtime_error(errno_message("setsockopt failed"));
        }
        peer_action(peer_fd);
    }
    catch (...)
    {
        peer_error = std::current_exception();
    }
    if (peer_fd >= 0)
    {
        ::close(peer_fd);
    }

    result.diagnostic = read_pipe(output_pipe[0]);
    ::close(output_pipe[0]);
    int status = 0;
    while (::waitpid(child, &status, 0) < 0)
    {
        if (errno != EINTR)
        {
            throw std::runtime_error(errno_message("waitpid failed"));
        }
    }
    if (WIFEXITED(status))
    {
        result.exit_code = WEXITSTATUS(status);
    }

    if (peer_error)
    {
        std::rethrow_exception(peer_error);
    }
    return result;
}

DriverResult run_driver_frame(const bool converged,
                              const bool cal_force,
                              const bool cal_stress,
                              const std::function<void(int)>& peer_action)
{
    return run_driver_session([&]() {
        UnitCell cell;
        initialize_one_atom_cell(cell);
        Input_para input;
        input.cal_force = cal_force;
        input.cal_stress = cal_stress;
        FakeESolver solver(converged);
        std::ofstream running("/dev/null");
        Socket_Driver driver;
        driver.socket_driver(&solver, cell, input, running);
        if (solver.electronic_context_seen || cell.electronic_step != nullptr)
        {
            throw std::runtime_error("ordinary evaluation acquired electronic request state");
        }
    }, [&](const int fd) {
        send_fixed_cell_frame(fd);
        peer_action(fd);
    });
}

DriverResult run_constant_potential_session(const std::function<void(int)>& peer_action,
                                            const ElectronicFailure failure)
{
    return run_driver_session([&]() {
        UnitCell cell;
        initialize_one_atom_cell(cell);
        Parameter parameters;
        Input_para& input = TestParameters::input(parameters);
        input.cal_force = true;
        input.cal_stress = false;
        input.symmetry = "0";
        input.esolver_type = "ksdft";
        input.basis_type = "pw";
        input.device = "cpu";
        input.precision = "double";
        input.nspin = 1;
        input.nbands = 8;
        input.smearing_method = "gaussian";
        ControlledESolver solver(failure);
        std::ofstream running("/dev/null");
        Socket_Driver driver;
        driver.socket_driver(&solver, cell, parameters, parameters.inp, parameters.globalv, running);
        if (cell.electronic_step != nullptr)
        {
            throw std::runtime_error("completed electronic request remained bound to the cell");
        }
        std::cout << "CP calls: " << solver.control_count << " controls, "
                  << solver.runner_count << " runners, " << solver.collection_count
                  << " collections\n";
    }, [&](const int fd) {
        send_header(fd, "STATUS");
        EXPECT_EQ("NEEDINIT", read_header_or_close(fd));
        const std::string metadata = "NELECT:8,neutral_electrons:8";
        const std::int32_t replica = 0;
        const std::int32_t nbytes = static_cast<std::int32_t>(metadata.size());
        send_header(fd, "INIT");
        send_value(fd, replica);
        send_value(fd, nbytes);
        send_all(fd, metadata.data(), metadata.size());
        peer_action(fd);
    });
}
} // namespace

TEST(SocketElectronicFrameTest, BindingIsScopedAndRejectsNestedRequests)
{
    UnitCell cell;
    elecstate::ElectronicStep first;
    elecstate::ElectronicStep second;
    EXPECT_EQ(nullptr, cell.electronic_step);
    {
        SocketCP::ScopedElectronicStep binding(cell, &first);
        EXPECT_EQ(&first, cell.electronic_step);
        EXPECT_THROW(SocketCP::ScopedElectronicStep nested(cell, &second), std::runtime_error);
        EXPECT_EQ(&first, cell.electronic_step);
    }
    EXPECT_EQ(nullptr, cell.electronic_step);
    try
    {
        SocketCP::ScopedElectronicStep binding(cell, &second);
        throw std::runtime_error("evaluation failed");
    }
    catch (const std::runtime_error&)
    {
        EXPECT_EQ(nullptr, cell.electronic_step);
    }
}

TEST(SocketElectronicFrameTest, SeparateCellsDoNotShareRequestState)
{
    UnitCell first_cell;
    UnitCell second_cell;
    elecstate::ElectronicStep first;
    elecstate::ElectronicStep second;
    first.control.nelect = 8.25;
    second.control.nelect = 7.75;
    SocketCP::ScopedElectronicStep first_binding(first_cell, &first);
    {
        SocketCP::ScopedElectronicStep second_binding(second_cell, &second);
        first_cell.electronic_step->request_applied = true;
        EXPECT_FALSE(second_cell.electronic_step->request_applied);
        EXPECT_DOUBLE_EQ(8.25, first_cell.electronic_step->control.nelect);
        EXPECT_DOUBLE_EQ(7.75, second_cell.electronic_step->control.nelect);
    }
    EXPECT_EQ(nullptr, second_cell.electronic_step);
    EXPECT_EQ(&first, first_cell.electronic_step);
}

TEST(SocketDriverTest, RejectedElectronicStepLeavesSolverUntouched)
{
    UnitCell cell;
    initialize_one_atom_cell(cell);
    Parameter parameters;
    // Configuration comes from the explicit references, not the opaque handle.
    Input_para input;
    System_para system;
    input.esolver_type = "ksdft";
    input.basis_type = "pw";
    input.device = "cpu";
    input.precision = "double";
    input.nspin = 1;
    input.nbands = 8;
    input.smearing_method = "gaussian";
    ControlledESolver solver(ElectronicFailure::None);
    elecstate::ElectronicStep request;
    request.parameters = &parameters;
    request.control.nelect = 8.25;
    SocketCP::ScopedElectronicStep binding(cell, &request);

    input.assume_isolated = "makov-payne";
    try
    {
        SocketCP::run_electronic_step(solver, cell, input, system, 2, 0);
        FAIL() << "Makov-Payne request reached the solver";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_THAT(error.what(), testing::HasSubstr("energy-only Makov-Payne"));
    }
    input.assume_isolated = "none";
    input.nbands = 4;
    EXPECT_THROW(SocketCP::run_electronic_step(solver, cell, input, system, 2, 0), std::runtime_error);
    input.nbands = 8;
    input.smearing_method = "fixed";
    EXPECT_THROW(SocketCP::run_electronic_step(solver, cell, input, system, 2, 0), std::runtime_error);
    input.smearing_method = "gaussian";
    system.use_uspp = true;
    EXPECT_THROW(SocketCP::run_electronic_step(solver, cell, input, system, 2, 0), std::runtime_error);
    system.use_uspp = false;
    system.two_fermi = true;
    EXPECT_THROW(SocketCP::run_electronic_step(solver, cell, input, system, 2, 0), std::runtime_error);
    system.two_fermi = false;
    EXPECT_THROW(SocketCP::run_electronic_step(solver, cell, input, system, 0, 0), std::runtime_error);
    EXPECT_EQ(0, solver.control_count);
    EXPECT_EQ(0, solver.runner_count);
    EXPECT_DOUBLE_EQ(0.0, solver.nelect);
    EXPECT_FALSE(request.request_applied);
    EXPECT_FALSE(request.result_ready);
    EXPECT_DOUBLE_EQ(0.0, parameters.inp.nelec);

    input.smearing_method = "gaussian";
    SocketCP::run_electronic_step(solver, cell, input, system, 2, 0);
    EXPECT_EQ(1, solver.control_count);
    EXPECT_EQ(1, solver.runner_count);
    EXPECT_DOUBLE_EQ(8.25, solver.nelect);

    EXPECT_TRUE(request.request_applied);
    EXPECT_FALSE(request.result_ready);
    EXPECT_THROW(SocketCP::run_electronic_step(solver, cell, input, system, 2, 1), std::runtime_error);
    EXPECT_EQ(1, solver.runner_count);
}

TEST(SocketDriverTest, ConstantPotentialFramesApplyChargeBeforeEvaluation)
{
    const DriverResult result = run_constant_potential_session([](const int fd) {
        const double counts[] = {8.25, 7.75};
        const double positions[] = {0.0, 0.125};
        const std::string extras[] = {"\"nelect\":8.25", "\"nelect\":7.75"};
        const std::int32_t flags = 0;
        for (int frame = 0; frame < 2; ++frame)
        {
            send_header(fd, "STATUS");
            EXPECT_EQ("READY", read_header_or_close(fd));
            send_header(fd, "CHGDATA");
            send_value(fd, counts[frame]);
            send_value(fd, flags);
            send_header(fd, "STATUS");
            EXPECT_EQ("READY", read_header_or_close(fd));
            send_positions(fd, positions[frame]);
            for (int query = 0; query < 2; ++query)
            {
                send_header(fd, "STATUS");
                EXPECT_EQ("HAVEDATA", read_header_or_close(fd));
            }
            send_header(fd, "GETFORCE");
            const ForceResponse response = read_force_response(fd);
            EXPECT_EQ("FORCEREADY", response.header);
            EXPECT_DOUBLE_EQ((4.0 + positions[frame] + counts[frame]) * 0.5, response.energy_hartree);
            ASSERT_EQ(3u, response.forces_hartree_per_bohr.size());
            EXPECT_DOUBLE_EQ(2.0 + positions[frame], response.forces_hartree_per_bohr[0]);
            EXPECT_THAT(response.extra, testing::HasSubstr(extras[frame]));
            EXPECT_THAT(response.extra, testing::HasSubstr("\"scf_converged\":true"));
        }
        send_header(fd, "STATUS");
        EXPECT_EQ("READY", read_header_or_close(fd));
        send_header(fd, "EXIT");
        EXPECT_EQ("", read_header_or_close(fd));
    }, ElectronicFailure::None);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(result.diagnostic, testing::HasSubstr("CP calls: 2 controls, 2 runners, 2 collections"));
}

TEST(SocketDriverTest, ConstantPotentialPositionRequiresFreshCharge)
{
    const DriverResult result = run_constant_potential_session([](const int fd) {
        send_header(fd, "POSDATA");
        EXPECT_EQ("", read_header_or_close(fd));
    }, ElectronicFailure::None);
    EXPECT_NE(0, result.exit_code);
    EXPECT_THAT(result.diagnostic, testing::HasSubstr("constant-potential POSDATA requires one fresh CHGDATA"));
}

TEST(SocketDriverTest, FailedElectronicFrameCannotReturnPreviousFeedback)
{
    const DriverResult result = run_constant_potential_session([](const int fd) {
        const std::int32_t flags = 0;
        const double first_nelect = 8.25;
        const double second_nelect = 7.75;
        send_header(fd, "CHGDATA");
        send_value(fd, first_nelect);
        send_value(fd, flags);
        send_positions(fd, 0.0);
        send_header(fd, "GETFORCE");
        const ForceResponse first = read_force_response(fd);
        EXPECT_EQ("FORCEREADY", first.header);
        EXPECT_THAT(first.extra, testing::HasSubstr("\"nelect\":8.25"));
        send_header(fd, "CHGDATA");
        send_value(fd, second_nelect);
        send_value(fd, flags);
        send_positions(fd, 0.125);
        EXPECT_EQ("", read_header_or_close(fd));
    }, ElectronicFailure::SecondFrame);
    EXPECT_NE(0, result.exit_code);
    EXPECT_THAT(result.diagnostic, testing::HasSubstr("SCF did not converge"));
}

TEST(SocketDriverTest, ConvergedSolverMustApplyRequestAndCompleteFeedback)
{
    const ElectronicFailure failures[] = {ElectronicFailure::IgnoreRequest, ElectronicFailure::MissingResult};
    for (const ElectronicFailure failure : failures)
    {
        const DriverResult result = run_constant_potential_session([](const int fd) {
            const double nelect = 8.25;
            const std::int32_t flags = 0;
            send_header(fd, "CHGDATA");
            send_value(fd, nelect);
            send_value(fd, flags);
            send_positions(fd, 0.0);
            EXPECT_EQ("", read_header_or_close(fd));
        }, failure);
        EXPECT_NE(0, result.exit_code);
        const char* diagnostic = "did not complete this frame's electronic observables";
        if (failure == ElectronicFailure::IgnoreRequest)
        {
            diagnostic = "did not apply the constant-potential electron request";
        }
        EXPECT_THAT(result.diagnostic, testing::HasSubstr(diagnostic));
    }
}

TEST(SocketElectronicFrameTest, CollectionConvertsPhysicalUnitsAndKeepsZeroFermi)
{
    elecstate::ElectronicStep step;
    step.request_applied = true;
    step.result_ready = true;
    step.result.nelect = 8.125;
    step.result.fermi_ry = 0.0;
    step.result.reference_ry = 0.25;
    const auto config = SocketCP::parse_init("NELECT:8,neutral_electrons:8,"
                                            "potential_average_axis:y,potential_average_min_A:2,potential_average_max_A:3");
    step.window = config.window;
    const auto values = SocketCP::collect_observables(config, step);
    EXPECT_DOUBLE_EQ(0.0, values.fermi_eV);
    EXPECT_DOUBLE_EQ(8.125, values.nelect);
    EXPECT_TRUE(values.has_workfunction);
    EXPECT_DOUBLE_EQ(0.25 * ModuleBase::Ry_to_eV, values.reference_eV);
    EXPECT_DOUBLE_EQ(values.reference_eV, values.workfunction_eV);
}

TEST(SocketElectronicFrameTest, CollectionRejectsMissingOrInvalidPhysicalObservables)
{
    elecstate::ElectronicStep step;
    const auto config = SocketCP::parse_init("NELECT:8,neutral_electrons:8");
    const SocketCP::Configuration ordinary;
    EXPECT_THROW(SocketCP::collect_observables(ordinary, step), std::runtime_error);
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
    step.result_ready = true;
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
    step.request_applied = true;
    step.result.nelect = 0.0;
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
    step.result.nelect = 8.0;
    step.result.fermi_ry = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(SocketCP::collect_observables(config, step), std::runtime_error);
}

TEST(SocketElectronicFrameTest, SerializationUsesPublishedSnapshotWithoutCollectingAgain)
{
    elecstate::ElectronicStep step;
    step.request_applied = true;
    step.result_ready = true;
    step.result.nelect = 8.125;
    step.result.fermi_ry = -0.125;
    const auto config = SocketCP::parse_init("NELECT:8,neutral_electrons:8");
    SocketUtils::ComputedFrame frame;
    frame.valid = true;
    frame.forces_present = true;
    frame.electronic = SocketCP::collect_observables(config, step);
    frame.electronic_present = true;
    step.result.nelect = 9.0;
    step.result_ready = false;
    const std::string first = SocketUtils::properties_extra(frame);
    const std::string second = SocketUtils::properties_extra(frame);
    EXPECT_EQ(first, second);
    EXPECT_THAT(first, testing::HasSubstr("\"nelect\":8.125"));
    EXPECT_THAT(first, testing::HasSubstr("\"present\":[\"energy\",\"forces\"]"));
    EXPECT_THAT(first, testing::Not(testing::HasSubstr("workfunction_eV")));
}

TEST(SocketDriverTest, NonconvergedFrameIsPublishedWithMetadata)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(false, true, false, [&](const int fd) {
        send_header(fd, "GETFORCE");
        response = read_force_response(fd);
    });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"scf_converged\":false"));
}

TEST(SocketDriverTest, EnergyOnlyFrameMarksForceAndStressAbsent)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(true, false, false, [&](const int fd) {
        send_header(fd, "GETFORCE");
        response = read_force_response(fd);
    });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"present\":[\"energy\"]"));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"forces\"")));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"stress\"")));
    EXPECT_THAT(response.forces_hartree_per_bohr, testing::ElementsAre(0.0, 0.0, 0.0));
    EXPECT_THAT(response.virial_wire_hartree, testing::ElementsAre(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
}

TEST(SocketDriverTest, EnergyAndStressFrameDoesNotAdvertiseForce)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(true, false, true, [&](const int fd) {
        send_header(fd, "GETFORCE");
        response = read_force_response(fd);
    });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"present\":[\"energy\",\"stress\"]"));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"forces\"")));
    EXPECT_NE(0.0, response.virial_wire_hartree[0]);
}

TEST(SocketDriverTest, EnergyAndForceFrameAdvertisesOnlyForce)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(true, true, false, [&](const int fd) {
        send_header(fd, "GETFORCE");
        response = read_force_response(fd);
    });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"present\":[\"energy\",\"forces\"]"));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"stress\"")));
    EXPECT_THAT(response.forces_hartree_per_bohr, testing::ElementsAre(2.0, 0.0, 0.0));
    EXPECT_THAT(response.virial_wire_hartree, testing::ElementsAre(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
}

TEST(SocketDriverTest, EnergyForceAndStressFrameAdvertisesBothDerivatives)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(true, true, true, [&](const int fd) {
        send_header(fd, "GETFORCE");
        response = read_force_response(fd);
    });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"present\":[\"energy\",\"forces\",\"stress\"]"));
    EXPECT_EQ(3u, response.forces_hartree_per_bohr.size());
    EXPECT_THAT(response.forces_hartree_per_bohr, testing::ElementsAre(2.0, 0.0, 0.0));
    EXPECT_NE(0.0, response.virial_wire_hartree[0]);
}

TEST(SocketDriverTest, ConsecutiveFramesKeepGeometryResultsAndConvergenceTogether)
{
    ForceResponse first, second;
    const DriverResult result = run_driver_frame(true, true, true, [&](const int fd) {
        send_header(fd, "GETFORCE");
        first = read_force_response(fd);
        send_positions(fd, 0.25);
        send_header(fd, "GETFORCE");
        second = read_force_response(fd);
    });
    EXPECT_EQ(0, result.exit_code);
    EXPECT_DOUBLE_EQ(2.0, first.energy_hartree);
    EXPECT_DOUBLE_EQ(2.125, second.energy_hartree);
    EXPECT_DOUBLE_EQ(2.0, first.forces_hartree_per_bohr.at(0));
    EXPECT_DOUBLE_EQ(2.25, second.forces_hartree_per_bohr.at(0));
    EXPECT_DOUBLE_EQ(1.0, first.virial_wire_hartree.at(0));
    EXPECT_DOUBLE_EQ(1.375, second.virial_wire_hartree.at(0));
    EXPECT_THAT(first.extra, testing::HasSubstr("\"scf_converged\":true"));
    EXPECT_THAT(second.extra, testing::HasSubstr("\"scf_converged\":false"));
}

TEST(SocketDriverTest, ConsumedFrameCannotBeReturnedTwice)
{
    const DriverResult result = run_driver_frame(true, true, false, [&](const int fd) {
        send_header(fd, "GETFORCE");
        read_force_response(fd);
        send_header(fd, "GETFORCE");
        EXPECT_EQ("", read_header_or_close(fd));
    });
    EXPECT_NE(0, result.exit_code);
    EXPECT_THAT(result.diagnostic, testing::HasSubstr("GETFORCE requires HAVEDATA"));
}

TEST(SocketDriverTest, OrdinaryExitMayDiscardUnconsumedFrame)
{
    const DriverResult result = run_driver_frame(true, true, false, [&](const int fd) {
        send_header(fd, "STATUS");
        EXPECT_EQ("HAVEDATA", read_header_or_close(fd));
        send_header(fd, "EXIT");
        EXPECT_EQ("", read_header_or_close(fd));
    });
    EXPECT_EQ(0, result.exit_code);
}

TEST(SocketDriverTest, ConstantPotentialExitRequiresCompletedRequest)
{
    SocketHandlers::DriverContext context;
    std::ofstream running("/dev/null");
    context.cp.enabled = true;
    context.state = SocketUtils::DriverState::HasData;
    EXPECT_THROW(SocketHandlers::handle_exit(context, running), std::runtime_error);
    context.state = SocketUtils::DriverState::Ready;
    context.charge_pending = true;
    EXPECT_THROW(SocketHandlers::handle_exit(context, running), std::runtime_error);
    context.charge_pending = false;
    EXPECT_NO_THROW(SocketHandlers::handle_exit(context, running));
}

TEST(SocketDriverTest, InvalidNextFrameCannotReturnPreviousResults)
{
    const DriverResult result = run_driver_frame(true, true, false, [&](const int fd) {
        send_header(fd, "GETFORCE");
        read_force_response(fd);
        send_positions(fd, std::numeric_limits<double>::quiet_NaN());
        EXPECT_EQ("", read_header_or_close(fd));
    });
    EXPECT_NE(0, result.exit_code);
    EXPECT_THAT(result.diagnostic, testing::HasSubstr("finite"));
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
