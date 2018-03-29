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
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const MPI_Count test_int_max = BigMPI_Get_max_int();

    if (size<2) {
        printf("Use 2 or more processes. \n");
        MPI_Finalize();
        return 1;
    }

    int l = (argc > 1) ? atoi(argv[1]) : 2;
    int m = (argc > 2) ? atoi(argv[2]) : 17777;
    MPI_Count n = l * test_int_max + m;

    char * buf = NULL;

    MPI_Alloc_mem((MPI_Aint)n, MPI_INFO_NULL, &buf);

    memset(buf, rank, (size_t)n);

    size_t errors = 0;
    for (int r = 1; r < size; r++) {

        /* pairwise communication */
        if (rank==r) {
            MPIX_Rsend_x(buf, n, MPI_CHAR, 0 /* dst */, r /* tag */, MPI_COMM_WORLD);
        }
        else if (rank==0) {
            MPIX_Recv_x(buf, n, MPI_CHAR, r /* src */, r /* tag */, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            errors += verify_buffer(buf, n, r);
            if (errors > 0) {
                printf("There were %zu errors!", errors);
            }
        }
    }

    MPI_Free_mem(buf);

    if (rank==0 && errors==0) {
        printf("SUCCESS\n");
    }

    MPI_Finalize();

    return 0;
}
