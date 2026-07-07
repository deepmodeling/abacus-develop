#ifndef ABACUS_IPI_SOCKET_H
#define ABACUS_IPI_SOCKET_H

#include <cstddef>
#include <string>
#include <vector>

class IpiSocket
{
  public:
    IpiSocket() = default;
    ~IpiSocket();

    IpiSocket(const IpiSocket&) = delete;
    IpiSocket& operator=(const IpiSocket&) = delete;

    void connect(const std::string& address);
    void close();

    std::string read_header();
    void write_header(const std::string& header);

    int read_int();
    void write_int(int value);

    double read_double();
    void write_double(double value);

    std::vector<double> read_doubles(std::size_t n);
    void write_doubles(const std::vector<double>& values);
    std::string read_string(std::size_t nbytes);

  private:
    int fd_ = -1;

    void read_exact(void* data, std::size_t nbytes);
    void write_exact(const void* data, std::size_t nbytes);
};

#endif
