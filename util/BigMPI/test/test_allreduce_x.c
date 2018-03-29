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

    double * sbuf = NULL;
    double * rbuf = NULL;

    MPI_Aint bytes = n*sizeof(double);
    MPI_Alloc_mem(bytes, MPI_INFO_NULL, &sbuf);
    MPI_Alloc_mem(bytes, MPI_INFO_NULL, &rbuf);

    for (MPI_Count i=0; i<n; i++) {
        sbuf[i] = (double)rank+1.;
    }
    for (MPI_Count i=0; i<n; i++) {
        rbuf[i] = 0.0;
    }

    /* collective communication */
    MPIX_Allreduce_x(sbuf, rbuf, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    size_t errors = 0;
    double val = (double)size*(size+1.)/2.;
    errors = verify_doubles(rbuf, n, val);
    if (errors) {
        printf("There were %zu errors out of %zu elements!\n", errors, (size_t)n);
        for (MPI_Count i=0; i<n; i++) {
            printf("rbuf[%zu] = %lf (expected %lf - %s)\n",
                    (size_t)i, rbuf[i], val, rbuf[i]==val ? "RIGHT" : "WRONG");
        }
        fflush(stdout);
    }

    MPI_Free_mem(sbuf);
    MPI_Free_mem(rbuf);

    /* TODO: reduce errors across all ranks in case root result is correct
     *       but others are wrong. */
    if (rank==0 && errors==0) {
        printf("SUCCESS\n");
    }

    MPI_Finalize();

    return 0;
}
