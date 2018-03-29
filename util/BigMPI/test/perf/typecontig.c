#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

const int bigmpi_int_max = INT_MAX;

int MPIX_Type_contiguous_x(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    MPI_Count c = count/bigmpi_int_max;
    MPI_Count r = count%bigmpi_int_max;

    MPI_Datatype chunks;
    MPI_Type_vector(c, bigmpi_int_max, bigmpi_int_max, oldtype, &chunks);

    MPI_Datatype remainder;
    MPI_Type_contiguous(r, oldtype, &remainder);

    MPI_Aint lb /* unused */, extent;
    MPI_Type_get_extent(oldtype, &lb, &extent);

    MPI_Aint remdisp          = (MPI_Aint)c*bigmpi_int_max*extent;
    int blocklengths[2]       = {1,1};
    MPI_Aint displacements[2] = {0,remdisp};
    MPI_Datatype types[2]     = {chunks,remainder};
    MPI_Type_create_struct(2, blocklengths, displacements, types, newtype);

    MPI_Type_free(&chunks);
    MPI_Type_free(&remainder);

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
        //MPIX_Type_contiguous_x((MPI_Count)i, MPI_DOUBLE, &(dtout[i]));
        //MPI_Type_commit(&(dtout[i]));
        MPI_Count bigcount = (i%10)*(MPI_Count)INT_MAX+(i%100);
        MPIX_Type_contiguous_x(bigcount, MPI_DOUBLE, &dtout);
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
