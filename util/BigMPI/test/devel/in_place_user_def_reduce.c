#include <stdio.h>
#include <mpi.h>

void printer(void * invec, void * inoutvec, int * len, MPI_Datatype * type)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("%d: invec = %p inoutvec = %p len = %d \n", rank, invec, inoutvec, *len);
    return;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int commute = 1;
    MPI_Op myop;
    MPI_Op_create(printer, commute, &myop);

    int in[8]  = {1,2,3,4,5,6,7,8};
    int out[8] = {-1,-2,-3,-4,-5,-6,-7,-8};

    printf("%d: MPI_IN_PLACE = %p in = %p out = %p \n", rank, MPI_IN_PLACE, in, out);

    MPI_Barrier(MPI_COMM_WORLD);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0) printf("out-of-place Allreduce \n");
    MPI_Allreduce(in, out, 8, MPI_INT, myop, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0) printf("in-place Allreduce \n");
    MPI_Allreduce(MPI_IN_PLACE, out, 8, MPI_INT, myop, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Op_free(&myop);
    MPI_Finalize();
    return 0;
}
