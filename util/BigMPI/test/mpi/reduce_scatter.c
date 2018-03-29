#include <stdio.h>
#include <unistd.h>

#include <mpi.h>

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size!=4) {
        if (rank==0) printf("Use 4 processes\n");
        MPI_Finalize();
        return size;
    }

    {
        if (rank==0) printf("MPI_Reduce_scatter(sendbuf, recvbuf...\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        int junk = rank+1;
        int sendbuf[4] = {junk, junk*2, junk*3, junk*4};
        int recvbuf[1] = {0};
        int recvcounts[4] = {1,1,1,1};
        MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        printf("%d: sendbuf = {%d,%d,%d,%d}, recvbuf = {%d} \n",
                rank, sendbuf[0], sendbuf[1], sendbuf[2], sendbuf[3], recvbuf[0]);
    }

    fflush(stdout);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("===================\n");

    {
        if (rank==0) printf("MPI_Reduce_scatter(MPI_IN_PLACE, recvbuf...\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        int junk = rank+1;
        int recvbuf[4] = {junk, junk*2, junk*3, junk*4};
        int recvcounts[4] = {1,1,1,1};
        MPI_Reduce_scatter(MPI_IN_PLACE, recvbuf, recvcounts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        printf("%d: recvbuf = {%d,%d,%d,%d} \n",
                rank, recvbuf[0], recvbuf[1], recvbuf[2], recvbuf[3]);
    }

    fflush(stdout);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("===================\n");

    {
        if (rank==0) printf("MPI_Reduce_scatter_block(sendbuf, recvbuf...\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        int junk = rank+1;
        int sendbuf[4] = {junk, junk*2, junk*3, junk*4};
        int recvbuf[1] = {0};
        int recvcount = 1;
        MPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        printf("%d: sendbuf = {%d,%d,%d,%d}, recvbuf = {%d} \n",
                rank, sendbuf[0], sendbuf[1], sendbuf[2], sendbuf[3], recvbuf[0]);
    }

    fflush(stdout);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("===================\n");

    {
        if (rank==0) printf("MPI_Reduce_scatter_block(MPI_IN_PLACE, recvbuf...\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        int junk = rank+1;
        int recvbuf[4] = {junk, junk*2, junk*3, junk*4};
        int recvcount = 1;
        MPI_Reduce_scatter_block(MPI_IN_PLACE, recvbuf, recvcount, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        printf("%d: recvbuf = {%d,%d,%d,%d} \n",
                rank, recvbuf[0], recvbuf[1], recvbuf[2], recvbuf[3]);
    }

    fflush(stdout);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("===================\n");

    {
        if (rank==0) printf("MPI_Reduce(sendbuf, tempbuf... + MPI_Scatter(tempbuf, recvcount...\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        int junk = rank+1;
        int sendbuf[4] = {junk, junk*2, junk*3, junk*4};
        int tempbuf[4] = {0,0,0,0};
        int recvbuf[1] = {0};
        int recvcount = 1;
        MPI_Reduce(sendbuf, tempbuf, 4*recvcount, MPI_INT, MPI_SUM, 0 /* root */, MPI_COMM_WORLD);
        MPI_Scatter(tempbuf, recvcount, MPI_INT, recvbuf, recvcount, MPI_INT, 0 /* root */, MPI_COMM_WORLD);
        printf("%d: sendbuf = {%d,%d,%d,%d}, recvbuf = {%d} \n",
                rank, sendbuf[0], sendbuf[1], sendbuf[2], sendbuf[3], recvbuf[0]);
    }

    fflush(stdout);
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("===================\n");

    {
        if (rank==0) printf("MPI_Reduce(MPI_IN_PLACE, recvbuf... + MPI_Scatter(MPI_IN_PLACE, recvcount...\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        int junk = rank+1;
        int recvbuf[4] = {junk, junk*2, junk*3, junk*4};
        int recvcount = 1;
        MPI_Reduce(rank==0 ? MPI_IN_PLACE : recvbuf, rank==0 ? recvbuf : NULL,
                   4*recvcount, MPI_INT, MPI_SUM, 0 /* root */, MPI_COMM_WORLD);
        MPI_Scatter(recvbuf, recvcount, MPI_INT, rank==0 ? MPI_IN_PLACE : recvbuf, recvcount, MPI_INT, 0 /* root */, MPI_COMM_WORLD);
        printf("%d: recvbuf = {%d,%d,%d,%d} \n",
                rank, recvbuf[0], recvbuf[1], recvbuf[2], recvbuf[3]);
    }

    MPI_Finalize();

    return 0;
}
