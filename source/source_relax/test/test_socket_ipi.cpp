#include "../socket_ipi.h"

#include "gtest/gtest.h"

#include <cerrno>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/socket.h>
#include <sys/un.h>
#include <thread>
#include <unistd.h>
#include <vector>

namespace
{
constexpr std::size_t kIpiHeaderLen = 12;

std::string errno_message(const std::string& prefix)
{
    return prefix + ": " + std::strerror(errno);
}

void send_all(int fd, const void* data, std::size_t nbytes)
{
    const char* cursor = static_cast<const char*>(data);
    std::size_t done = 0;
    while (done < nbytes)
    {
        const ssize_t sent = ::send(fd, cursor + done, nbytes - done, 0);
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

void recv_all(int fd, void* data, std::size_t nbytes)
{
    char* cursor = static_cast<char*>(data);
    std::size_t done = 0;
    while (done < nbytes)
    {
        const ssize_t received = ::recv(fd, cursor + done, nbytes - done, 0);
        if (received < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("recv failed"));
        }
        if (received == 0)
        {
            throw std::runtime_error("socket closed while receiving test data");
        }
        done += static_cast<std::size_t>(received);
    }
}

std::string padded_header(const std::string& header)
{
    std::string padded = header;
    padded.resize(kIpiHeaderLen, ' ');
    return padded;
}

class UnixSocketServer
{
  public:
    UnixSocketServer()
    {
        char dir_template[] = "/tmp/abacus_socket_ipi_test_XXXXXX";
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

        sockaddr_un addr;
        std::memset(&addr, 0, sizeof(addr));
        addr.sun_family = AF_UNIX;
        std::strncpy(addr.sun_path, path_.c_str(), sizeof(addr.sun_path) - 1);
        if (::bind(listen_fd_, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) != 0)
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

    int accept_once()
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

void rethrow_thread_error(const std::exception_ptr& thread_error)
{
    if (thread_error)
    {
        std::rethrow_exception(thread_error);
    }
}
} // namespace

TEST(IpiSocketTest, WriteHeaderPadsToTwelveBytes)
{
    UnixSocketServer server;
    std::string received;
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            char buffer[kIpiHeaderLen];
            recv_all(fd, buffer, sizeof(buffer));
            received.assign(buffer, sizeof(buffer));
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    socket.write_header("READY");
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
    EXPECT_EQ(padded_header("READY"), received);
}

TEST(IpiSocketTest, CleanPeerCloseBeforeNextHeaderThrowsDedicatedSignal)
{
    UnixSocketServer server;
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            const std::string header = padded_header("STATUS");
            send_all(fd, header.data(), header.size());
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    EXPECT_EQ("STATUS", socket.read_header());
    EXPECT_THROW(socket.read_header(), IpiSocketClosed);
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
}

TEST(IpiSocketTest, PartialHeaderCloseStaysRuntimeError)
{
    UnixSocketServer server;
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            const std::string partial = "STAT";
            send_all(fd, partial.data(), partial.size());
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    try
    {
        static_cast<void>(socket.read_header());
        FAIL() << "partial header EOF should throw";
    }
    catch (const IpiSocketClosed&)
    {
        FAIL() << "partial header EOF must not be treated as clean peer close";
    }
    catch (const std::runtime_error& exc)
    {
        EXPECT_NE(std::string::npos, std::string(exc.what()).find("closed while reading header"));
    }
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
}

TEST(IpiSocketTest, Int32UsesExactlyFourNativeEndianBytes)
{
    UnixSocketServer server;
    const std::int32_t expected = INT32_C(0x12345678);
    std::vector<char> received(4);
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            recv_all(fd, received.data(), received.size());
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    socket.write_int32(expected);
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
    EXPECT_EQ(0, std::memcmp(received.data(), &expected, 4));
}

TEST(IpiSocketTest, DoubleUsesExactlyEightNativeEndianBytes)
{
    UnixSocketServer server;
    const double expected = -1234.5;
    std::vector<char> received(8);
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            recv_all(fd, received.data(), received.size());
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    socket.write_double(expected);
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
    EXPECT_EQ(0, std::memcmp(received.data(), &expected, 8));
}

TEST(IpiSocketTest, ReadInt32HandlesSplitPayload)
{
    UnixSocketServer server;
    const std::int32_t expected = INT32_C(0x12345678);
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            const char* bytes = reinterpret_cast<const char*>(&expected);
            send_all(fd, bytes, 2);
            send_all(fd, bytes + 2, 2);
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    EXPECT_EQ(expected, socket.read_int32());
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
}

TEST(IpiSocketTest, ReadInt32RejectsMidPayloadClose)
{
    UnixSocketServer server;
    const std::int32_t value = INT32_C(0x12345678);
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            send_all(fd, &value, 2);
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    EXPECT_THROW(socket.read_int32(), IpiSocketClosed);
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
}

TEST(IpiSocketTest, WriteDoublesCompletesLargePayloadWithSmallPeerReads)
{
    UnixSocketServer server;
    std::vector<double> expected(1 << 18);
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        expected[i] = -1234.5 + static_cast<double>(i) * 0.25;
    }
    std::vector<char> received(expected.size() * sizeof(double));
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            std::size_t done = 0;
            while (done < received.size())
            {
                const std::size_t remaining = received.size() - done;
                const std::size_t chunk = remaining < 37 ? remaining : 37;
                const ssize_t nread = ::recv(fd, received.data() + done, chunk, 0);
                if (nread < 0)
                {
                    if (errno == EINTR)
                    {
                        continue;
                    }
                    throw std::runtime_error(errno_message("recv failed"));
                }
                if (nread == 0)
                {
                    throw std::runtime_error("socket closed while receiving large payload");
                }
                done += static_cast<std::size_t>(nread);
            }
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    socket.write_doubles(expected);
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
    EXPECT_EQ(0, std::memcmp(received.data(), expected.data(), received.size()));
}

TEST(IpiSocketTest, ReadDoublesRejectsByteCountOverflow)
{
    IpiSocket socket;
    const std::size_t count = std::numeric_limits<std::size_t>::max() / sizeof(double) + 1;

    try
    {
        static_cast<void>(socket.read_doubles(count));
        FAIL() << "overflowing double payload size should throw";
    }
    catch (const std::overflow_error& exc)
    {
        EXPECT_NE(std::string::npos, std::string(exc.what()).find(std::to_string(count)));
    }
    catch (...)
    {
        FAIL() << "overflowing double payload size should throw std::overflow_error";
    }
}

TEST(IpiSocketTest, WriteStringSendsExactBytesWithoutTerminator)
{
    UnixSocketServer server;
    const std::string expected = "{\"scf_converged\":false}";
    std::vector<char> received(expected.size());
    std::exception_ptr thread_error;
    std::thread peer([&]() {
        try
        {
            const int fd = server.accept_once();
            recv_all(fd, received.data(), received.size());
            ::close(fd);
        }
        catch (...)
        {
            thread_error = std::current_exception();
        }
    });

    IpiSocket socket;
    socket.connect(server.address());
    socket.write_string(expected);
    socket.close();

    peer.join();
    rethrow_thread_error(thread_error);
    EXPECT_EQ(expected, std::string(received.begin(), received.end()));
}
