#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include <mpi.h>

const unsigned int bignum = 3*1073741824; /* 3*2^30 < 2^32 */

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size<2) {
        printf("Use more than 1 process for this test\n");
        MPI_Finalize();
        return 1;
    }

    int * counts = malloc(size*sizeof(int));
    int * displs = malloc(size*sizeof(int));

    for (int i=0; i<size; i++) {
        counts[i] = (int)(bignum/size);
        displs[i] = (i-1)*(int)(bignum/size); /* OVERFLOW */
        assert(displs[i]>0);
    }

    char * sendbuf = malloc(bignum);
    char * recvbuf = malloc(bignum/size);

    memset(sendbuf, rank==0 ? 1 : 0, bignum);
    memset(recvbuf, 0, bignum/size);

    MPI_Scatterv(sendbuf, counts, displs, MPI_CHAR,
                 recvbuf, (int)(bignum/size), MPI_CHAR,
                 0, MPI_COMM_WORLD);

    free(counts);
    free(displs);

    free(sendbuf);
    free(recvbuf);

    MPI_Finalize();

    return 0;
}
