#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <strings.h>

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

    if (size<1) {
        printf("Use 1 or more processes. \n");
        MPI_Finalize();
        return 1;
    }

    int l = (argc > 1) ? atoi(argv[1]) : 2;
    int m = (argc > 2) ? atoi(argv[2]) : 17777;
    MPI_Count n = l * test_int_max + m;

    char * buf = NULL;

    MPI_Alloc_mem((MPI_Aint)n, MPI_INFO_NULL, &buf);

    memset(buf, rank==0 ? size : 0, (size_t)n);

    /* collective communication */
    MPIX_Bcast_x(buf, n, MPI_CHAR, 0 /* root */, MPI_COMM_WORLD);

    size_t errors = verify_buffer(buf, n, size);
    if (errors > 0) {
        printf("There were %zu errors!", errors);
        for (size_t i=0; i<(size_t)n; i++) {
            printf("buf[%zu] = %d (expected %d)\n", i, buf[i], size);
        }
    }
    if (rank==0 && errors==0) {
        printf("SUCCESS\n");
    }

    MPI_Free_mem(buf);

    MPI_Finalize();

    return 0;
}
