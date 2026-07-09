#include "source_relax/ipi_socket.h"

#include <arpa/inet.h>
#include <cerrno>
#include <cstring>
#include <netdb.h>
#include <stdexcept>
#include <string>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

namespace
{
constexpr std::size_t IPI_HEADER_LEN = 12;

std::string errno_message(const std::string& prefix)
{
    return prefix + ": " + std::strerror(errno);
}

std::string trim_header(const char* data)
{
    std::string value(data, IPI_HEADER_LEN);
    while (!value.empty() && value.back() == ' ')
    {
        value.pop_back();
    }
    return value;
}

std::string padded_header(const std::string& header)
{
    if (header.size() > IPI_HEADER_LEN)
    {
        throw std::runtime_error("i-PI header is longer than 12 bytes: " + header);
    }
    std::string out = header;
    out.resize(IPI_HEADER_LEN, ' ');
    return out;
}
} // namespace

IpiSocketClosed::IpiSocketClosed(const std::string& message) : std::runtime_error(message)
{
}

IpiSocket::~IpiSocket()
{
    this->close();
}

void IpiSocket::connect(const std::string& address)
{
    this->close();
    const std::size_t colon = address.rfind(':');
    if (colon == std::string::npos)
    {
        throw std::runtime_error("i-PI address must be host:port or path:UNIX, got " + address);
    }
    const std::string host = address.substr(0, colon);
    const std::string service = address.substr(colon + 1);

    if (service == "UNIX")
    {
        fd_ = ::socket(AF_UNIX, SOCK_STREAM, 0);
        if (fd_ < 0)
        {
            throw std::runtime_error(errno_message("failed to create UNIX socket"));
        }
        sockaddr_un addr;
        std::memset(&addr, 0, sizeof(addr));
        addr.sun_family = AF_UNIX;
        if (host.size() >= sizeof(addr.sun_path))
        {
            this->close();
            throw std::runtime_error("UNIX socket path too long: " + host);
        }
        std::strncpy(addr.sun_path, host.c_str(), sizeof(addr.sun_path) - 1);
        if (::connect(fd_, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) != 0)
        {
            const std::string msg = errno_message("failed to connect UNIX i-PI socket " + host);
            this->close();
            throw std::runtime_error(msg);
        }
        return;
    }

    addrinfo hints;
    std::memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;

    addrinfo* result = nullptr;
    const int gai = ::getaddrinfo(host.c_str(), service.c_str(), &hints, &result);
    if (gai != 0)
    {
        throw std::runtime_error("failed to resolve i-PI socket " + address + ": " + ::gai_strerror(gai));
    }

    std::string last_error;
    for (addrinfo* rp = result; rp != nullptr; rp = rp->ai_next)
    {
        fd_ = ::socket(rp->ai_family, rp->ai_socktype, rp->ai_protocol);
        if (fd_ < 0)
        {
            last_error = errno_message("failed to create INET socket");
            continue;
        }
        if (::connect(fd_, rp->ai_addr, rp->ai_addrlen) == 0)
        {
            ::freeaddrinfo(result);
            return;
        }
        last_error = errno_message("failed to connect INET i-PI socket " + address);
        this->close();
    }
    ::freeaddrinfo(result);
    throw std::runtime_error(last_error.empty() ? "failed to connect i-PI socket " + address : last_error);
}

void IpiSocket::close()
{
    if (fd_ >= 0)
    {
        ::close(fd_);
        fd_ = -1;
    }
}

std::string IpiSocket::read_header()
{
    char header[IPI_HEADER_LEN];
    std::size_t done = 0;
    while (done < sizeof(header))
    {
        const ssize_t nread = ::recv(fd_, header + done, sizeof(header) - done, 0);
        if (nread == 0)
        {
            if (done == 0)
            {
                throw IpiSocketClosed("i-PI socket closed before next header");
            }
            throw std::runtime_error("i-PI socket closed while reading header");
        }
        if (nread < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("i-PI socket header read failed"));
        }
        done += static_cast<std::size_t>(nread);
    }
    return trim_header(header);
}

void IpiSocket::write_header(const std::string& header)
{
    const std::string padded = padded_header(header);
    this->write_exact(padded.data(), padded.size());
}

int IpiSocket::read_int()
{
    int value = 0;
    this->read_exact(&value, sizeof(value));
    return value;
}

void IpiSocket::write_int(int value)
{
    this->write_exact(&value, sizeof(value));
}

double IpiSocket::read_double()
{
    double value = 0.0;
    this->read_exact(&value, sizeof(value));
    return value;
}

void IpiSocket::write_double(double value)
{
    this->write_exact(&value, sizeof(value));
}

std::vector<double> IpiSocket::read_doubles(std::size_t n)
{
    std::vector<double> values(n);
    if (!values.empty())
    {
        this->read_exact(values.data(), values.size() * sizeof(double));
    }
    return values;
}

void IpiSocket::write_doubles(const std::vector<double>& values)
{
    if (!values.empty())
    {
        this->write_exact(values.data(), values.size() * sizeof(double));
    }
}

std::string IpiSocket::read_string(std::size_t nbytes)
{
    std::string value(nbytes, '\0');
    if (nbytes > 0)
    {
        this->read_exact(&value[0], nbytes);
    }
    return value;
}

void IpiSocket::read_exact(void* data, std::size_t nbytes)
{
    char* cursor = static_cast<char*>(data);
    std::size_t done = 0;
    while (done < nbytes)
    {
        const ssize_t nread = ::recv(fd_, cursor + done, nbytes - done, 0);
        if (nread == 0)
        {
            throw IpiSocketClosed("i-PI socket closed while reading");
        }
        if (nread < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("i-PI socket read failed"));
        }
        done += static_cast<std::size_t>(nread);
    }
}

void IpiSocket::write_exact(const void* data, std::size_t nbytes)
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
        const ssize_t nwritten = ::send(fd_, cursor + done, nbytes - done, flags);
        if (nwritten == 0)
        {
            throw std::runtime_error("i-PI socket closed while writing");
        }
        if (nwritten < 0)
        {
            if (errno == EINTR)
            {
                continue;
            }
            throw std::runtime_error(errno_message("i-PI socket write failed"));
        }
        done += static_cast<std::size_t>(nwritten);
    }
}
