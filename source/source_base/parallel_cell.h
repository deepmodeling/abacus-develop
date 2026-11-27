#ifndef PARALLEL_CELL_H
#define PARALLEL_CELL_H

#include <string>

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleBase
{
class CommunicationDomain
{
public:
    CommunicationDomain() = default;
#ifdef __MPI
    void initialize(MPI_Comm communicator);
    MPI_Comm communicator() const;
#endif
    int rank() const;
    int size() const;
    int max(int value) const;
    int sum(int value) const;
    void sum(int* values, int count) const;
    void sum(double* values, int count) const;
    double max(double value) const;
    /// Gather one bounded text batch on rank zero; other ranks receive an empty string.
    /// maximum_count bounds both local payloads and their total in MPI count units.
    std::string gather_text(const std::string& value, int maximum_count) const;

private:
#ifdef __MPI
    MPI_Comm communicator_ = MPI_COMM_NULL;
#endif
    int rank_ = 0;
};

CommunicationDomain world_comm_domain();
} // namespace ModuleBase

#endif
