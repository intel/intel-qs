#ifndef IQS_MPI_UTILS_HPP
#define IQS_MPI_UTILS_HPP

#ifdef INTELQS_HAS_MPI

/////////////////////////////////////////////////////////////////////////////////////////
// All methods involve MPI types in the arguments. Only available with MPI enabled.
// However there are two implementation depending on whether BIGMPI is used or not.
/////////////////////////////////////////////////////////////////////////////////////////

#include <complex>
#include <mpi.h>

#include "utils.hpp"

#ifdef BIGMPI
#include "bigmpi.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

namespace mpi {

/////////////////////////////////////////////////////////////////////////////////////////

#ifndef BIGMPI

/////////////////////////////////////////////////////////////////////////////////////////
// Definitions without BigMPI
/////////////////////////////////////////////////////////////////////////////////////////

static int MPI_Allreduce_x(float *sendbuf, float *recvbuf, int count,
                           MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce((void*)sendbuf, (void *)recvbuf, count, MPI_FLOAT, op, comm);
}

static int MPI_Allreduce_x(double *sendbuf, double *recvbuf, int count,
                           MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce((void*)sendbuf, (void *)recvbuf, count, MPI_DOUBLE, op, comm);
}

/////////////////////////////////////////////////////////////////////////////////////////

static
int MPI_Sendrecv_x(ComplexSP *sendbuf, size_t sendcount, size_t dest, size_t sendtag,
                   ComplexSP *recvbuf, size_t recvcount, size_t source, size_t recvtag,
                   MPI_Comm comm, MPI_Status *status)
{
  return MPI_Sendrecv((void *)sendbuf, sendcount, MPI_COMPLEX, dest, sendtag,
                      (void *)recvbuf, recvcount, MPI_COMPLEX, source, recvtag,
                       comm, status);
}

static
int MPI_Sendrecv_x(ComplexDP *sendbuf, size_t sendcount, size_t dest, size_t sendtag,
                   ComplexDP *recvbuf, size_t recvcount, size_t source, size_t recvtag,
                   MPI_Comm comm, MPI_Status *status)
{
  return MPI_Sendrecv((void *)sendbuf, sendcount, MPI_DOUBLE_COMPLEX, dest, sendtag,
                      (void *)recvbuf, recvcount, MPI_DOUBLE_COMPLEX, source, recvtag,
                      comm, status);
}

/////////////////////////////////////////////////////////////////////////////////////////

static int MPI_Bcast_x(ComplexSP *data, int root, MPI_Comm comm)
{
  return MPI_Bcast((void*)data, 1, MPI_COMPLEX, root, comm);
}

static int MPI_Bcast_x(ComplexDP *data, int root, MPI_Comm comm)
{
  return MPI_Bcast((void*)data, 1, MPI_DOUBLE_COMPLEX, root, comm);
}

#else

/////////////////////////////////////////////////////////////////////////////////////////
// Definitions with BigMPI
/////////////////////////////////////////////////////////////////////////////////////////

static int MPI_Allreduce_x(float *sendbuf, float *recvbuf, int count,
                           MPI_Op op, MPI_Comm comm)
{
  return MPIX_Allreduce_x((void*)sendbuf, (void *)recvbuf, count, MPI_FLOAT, op, comm);
}

static int MPI_Allreduce_x(double *sendbuf, double *recvbuf, int count,
                           MPI_Op op, MPI_Comm comm)
{
  return MPIX_Allreduce_x((void*)sendbuf, (void *)recvbuf, count, MPI_DOUBLE, op, comm);
}

/////////////////////////////////////////////////////////////////////////////////////////

static
 int MPI_Sendrecv_x(ComplexSP *sendbuf, size_t sendcount, size_t dest, size_t sendtag,
                   ComplexSP *recvbuf, size_t recvcount, size_t source, size_t recvtag,
                   MPI_Comm comm, MPI_Status *status)
{
  return MPIX_Sendrecv_x((void *)sendbuf, sendcount, MPI_COMPLEX, dest, sendtag,
                         (void *)recvbuf, recvcount, MPI_COMPLEX, source, recvtag,
                         comm, status);
}

static
int MPI_Sendrecv_x(ComplexDP *sendbuf, size_t sendcount, size_t dest, size_t sendtag,
                   ComplexDP *recvbuf, size_t recvcount, size_t source, size_t recvtag,
                   MPI_Comm comm, MPI_Status *status)
{
  return MPIX_Sendrecv_x((void *)sendbuf, sendcount, MPI_DOUBLE_COMPLEX, dest, sendtag,
                         (void *)recvbuf, recvcount, MPI_DOUBLE_COMPLEX, source, recvtag,
                         comm, status);
}

/////////////////////////////////////////////////////////////////////////////////////////

static int MPI_Bcast_x(ComplexSP *data, int root, MPI_Comm comm)
{
 // Not sure of how it is defined with BigMPI.
 assert(0);
}

//

static int MPI_Bcast_x(ComplexDP *data, int root, MPI_Comm comm)
{
 // Not sure of how it is defined with BigMPI.
 assert(0);
}

#endif //BIGMPI

/////////////////////////////////////////////////////////////////////////////////////////

}	// end namespace mpi
}	// end namespace qhipster

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// INTELQS_HAS_MPI

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard IQS_MPI_UTILS_HPP
