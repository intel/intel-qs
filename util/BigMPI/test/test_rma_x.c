#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include <mpi.h>
#include "bigmpi.h"
#include "verify_buffer.h"

/* Yes, it is technically unsafe to cast MPI_Count to MPI_Aint or size_t without checking,
 * given that MPI_Count might be 128b and MPI_Aint and size_t might be 64b, but BigMPI
 * does not aspire to support communication of more than 8 EiB messages at a time. */

int main(int argc, char * argv[])
{
    const MPI_Count test_int_max = BigMPI_Get_max_int();

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int l = (argc > 1) ? atoi(argv[1]) : 2;
    int m = (argc > 2) ? atoi(argv[2]) : 17777;
    MPI_Count n = l * test_int_max + m;

    double * baseptr = NULL;
    MPI_Win win;
    MPI_Aint winsize = (rank==0 ? n : 0);
#if MPI_VERSION >= 3
    /* Allocate all the window memory on rank 0 */
    MPI_Win_allocate(winsize, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &baseptr, &win);
    MPI_Win_lock_all(0, win);
#else
    MPI_Alloc_mem(winsize*sizeof(double), MPI_INFO_NULL, &baseptr);
    MPI_Win_create(baseptr, size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
#endif

    if (rank==0) {
        for (size_t i=0; i<(n/sizeof(double)); i++) {
            baseptr[i] = 0.0;
        }
#if MPI_VERSION >= 3
        MPI_Win_sync(win);
#endif
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double * buf = NULL;
    MPI_Alloc_mem((MPI_Aint)n, MPI_INFO_NULL, &buf);
    for (size_t i=0; i<(n/sizeof(double)); i++) {
        buf[i] = 1.0;
    }

#if MPI_VERSION < 3
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0 /* target */, 0 /* assert */, win);
#endif
    MPIX_Accumulate_x(buf, n/sizeof(double), MPI_DOUBLE,
                      0 /* target */, 0 /* disp */, n/sizeof(double), MPI_DOUBLE, MPI_SUM, win);
#if MPI_VERSION >= 3
    MPI_Win_flush(0,win);
#else
    MPI_Win_unlock(0,win);
#endif

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0) {
#if MPI_VERSION >= 3
        MPI_Win_sync(win);
#endif
        double expected = size;
        size_t errors = verify_doubles(baseptr, n/sizeof(double), expected);
        if (errors > 0) {
            printf("There were %zu errors!", errors);
            for (size_t i=0; i<(n/(sizeof(double))); i++) {
                printf("baseptr[%zu] = %lf (expected %lf)\n", i, baseptr[i], expected);
            }
        }
        if (errors==0) {
            printf("SUCCESS\n");
        }
    }

#if MPI_VERSION >= 3
    MPI_Win_unlock_all(win);
#endif
    MPI_Win_free(&win);
#if MPI_VERSION < 3
    MPI_Free_mem(baseptr);
#endif

    MPI_Finalize();

    return 0;
}
