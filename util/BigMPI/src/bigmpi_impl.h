#ifndef BIGMPI_IMPL_H
#define BIGMPI_IMPL_H

#include "bigmpiconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>

#include <mpi.h>

#include "bigmpi.h"

#include "likely.h"

#ifdef BIGMPI_MAX_INT
static const MPI_Count bigmpi_int_max   = BIGMPI_MAX_INT;
static const MPI_Count bigmpi_count_max = (MPI_Count)BIGMPI_MAX_INT*BIGMPI_MAX_INT;
#else
#include <limits.h>
#include <stdint.h>
static const MPI_Count bigmpi_int_max   = INT_MAX;
/* SIZE_MAX corresponds to size_t, which should be what MPI_Aint is. */
static const MPI_Count bigmpi_count_max = SIZE_MAX;
#endif

void BigMPI_Error_impl(const char *file, const int line, const char *func, const char *msg, ...);

#define BigMPI_Error(...) BigMPI_Error_impl(__FILE__,__LINE__,__func__,__VA_ARGS__)

void BigMPI_Convert_vectors(int                num,
                            int                splat_old_count,
                            const MPI_Count    oldcount,
                            const MPI_Count    oldcounts[],
                            int                splat_old_type,
                            const MPI_Datatype oldtype,
                            const MPI_Datatype oldtypes[],
                            int                zero_new_displs,
                            const MPI_Aint     olddispls[],
                                  int          newcounts[],
                                  MPI_Datatype newtypes[],
                                  MPI_Aint     newdispls[]);

#endif // BIGMPI_IMPL_H
