#include "source_relax/socket_driver.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mpi.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"
#include "for_test.h"

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
        position_ = dynamic_cast<UnitCell&>(cell).atoms[0].tau[0].x;
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
        force.create(cell.nat(), 3);
        force(0, 0) = 4.0 + 2.0 * position_;
    }

    void cal_stress(BaseCell&, ModuleBase::matrix& stress) override
    {
        stress.create(3, 3);
        stress(0, 0) = 2.0 + 3.0 * position_;
        stress(1, 1) = 2.0;
        stress(2, 2) = 2.0;
    }

  private:
    double position_ = 0.0;
    bool converged_;
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
    const double identity[9] = {1.0, 0.0, 0.0,
                                0.0, 1.0, 0.0,
                                0.0, 0.0, 1.0};
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
    response.forces_hartree_per_bohr
        = read_doubles(fd, static_cast<std::size_t>(3 * response.nat));
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

DriverResult run_driver_frame(const bool converged,
                              const bool cal_force,
                              const bool cal_stress,
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

        UnitCell ucell;
        initialize_one_atom_cell(ucell);
        Input_para input;
        input.cal_force = cal_force;
        input.cal_stress = cal_stress;
        FakeESolver solver(converged);
        std::ofstream running("/dev/null");
        Socket_Driver driver;
        driver.socket_driver(&solver, ucell, input, running);
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
        send_fixed_cell_frame(peer_fd);
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
} // namespace

TEST(SocketDriverTest, NonconvergedFrameIsPublishedWithMetadata)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(
        false, true, false,
        [&](const int fd) {
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
    const DriverResult result = run_driver_frame(
        true, false, false,
        [&](const int fd) {
            send_header(fd, "GETFORCE");
            response = read_force_response(fd);
        });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"present\":[\"energy\"]"));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"forces\"")));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"stress\"")));
    EXPECT_THAT(response.forces_hartree_per_bohr,
                testing::ElementsAre(0.0, 0.0, 0.0));
    EXPECT_THAT(response.virial_wire_hartree,
                testing::ElementsAre(0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0));
}

TEST(SocketDriverTest, EnergyAndStressFrameDoesNotAdvertiseForce)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(
        true, false, true,
        [&](const int fd) {
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
    const DriverResult result = run_driver_frame(
        true, true, false,
        [&](const int fd) {
            send_header(fd, "GETFORCE");
            response = read_force_response(fd);
        });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra, testing::HasSubstr("\"present\":[\"energy\",\"forces\"]"));
    EXPECT_THAT(response.extra, testing::Not(testing::HasSubstr("\"stress\"")));
    EXPECT_THAT(response.forces_hartree_per_bohr,
                testing::ElementsAre(2.0, 0.0, 0.0));
    EXPECT_THAT(response.virial_wire_hartree,
                testing::ElementsAre(0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0));
}

TEST(SocketDriverTest, EnergyForceAndStressFrameAdvertisesBothDerivatives)
{
    ForceResponse response;
    const DriverResult result = run_driver_frame(
        true, true, true,
        [&](const int fd) {
            send_header(fd, "GETFORCE");
            response = read_force_response(fd);
        });

    EXPECT_EQ("FORCEREADY", response.header);
    EXPECT_EQ(0, result.exit_code);
    EXPECT_THAT(response.extra,
                testing::HasSubstr("\"present\":[\"energy\",\"forces\",\"stress\"]"));
    EXPECT_EQ(3u, response.forces_hartree_per_bohr.size());
    EXPECT_THAT(response.forces_hartree_per_bohr,
                testing::ElementsAre(2.0, 0.0, 0.0));
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
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
