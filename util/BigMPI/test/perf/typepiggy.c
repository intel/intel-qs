#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int          piggysize = 8;
MPI_Datatype piggytype = MPI_CHAR;
char piggybuf[8] = "01234567";

int MPIX_Type_piggy_x(int count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    int size;
    MPI_Type_size(oldtype,&size);

    int blocklengths[2]       = {count,piggysize};
    MPI_Aint displacements[2] = {0,(MPI_Aint)count * (MPI_Aint)size};
    MPI_Datatype types[2]     = {oldtype,piggytype};
    MPI_Type_create_struct(2, blocklengths, displacements, types, newtype);

    return MPI_SUCCESS;
}

int main(int argc, char* argv[])
{
    int rank=0, size=1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = (argc>1) ? atoi(argv[1]) : 10000;
    //MPI_Datatype * dtout = malloc(n*sizeof(MPI_Datatype));
    MPI_Datatype dtout;
    double t0 = MPI_Wtime();
    for (int i=0; i<n; i++) {
        MPIX_Type_piggy_x(i, MPI_DOUBLE, &dtout);
        MPI_Type_commit(&dtout);
        MPI_Type_free(&dtout);
    }
    double t1 = MPI_Wtime();
    double dt = t1-t0;
    printf("create, commit (free?) %d Type_contig_x in %lf s (%lf us per call)\n",
            n, dt, 1.e6*dt/(double)n);

    //for (int i=0; i<n; i++) {
    //    MPI_Type_free(&(dtout[i]));
    //}
    //free(dtout);

    MPI_Finalize();
    return 0;
}
