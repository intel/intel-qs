#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <mpi.h>

int MPIX_Type_contiguous_x(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    const MPI_Count bigmpi_int_max = (MPI_Count)INT_MAX;

    MPI_Count c = count/bigmpi_int_max;
    MPI_Count r = count%bigmpi_int_max;

    assert(c<bigmpi_int_max);
    assert(r<bigmpi_int_max);

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

int main(int argc, char** argv)
{
    int requested=MPI_THREAD_SINGLE, provided;
    MPI_Init_thread(&argc,&argv,requested,&provided);

    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Aint n = (argc>1) ? atol(argv[1]) : 1L<<31;
    if (rank==0) {
        printf("Test: %zu bytes\n",n);
    }

    char * ptr = NULL;
    MPI_Alloc_mem(n,MPI_INFO_NULL,&ptr);

    if (rank==0) {
        memset(ptr,255,n);
    } else {
        memset(ptr,0,n);
    }

    MPI_Datatype bigtype;
    MPIX_Type_contiguous_x(n,MPI_CHAR,&bigtype);
    MPI_Type_commit(&bigtype);

    if (0 && n<(size_t)INT_MAX) {
        if (rank==0) printf("Not using a contiguous datatype\n");
        MPI_Bcast(ptr,(int)n,MPI_CHAR,0,MPI_COMM_WORLD);
    } else {
        if (rank==0) printf("Using a contiguous datatype\n");
        MPI_Bcast(ptr,1,bigtype,0,MPI_COMM_WORLD);
    }

    size_t errors = 0;
    for (size_t i=0; i<n; i++) {
        errors += (ptr[i]!=(char)255);
    }

    if (errors>0) {
        printf("%d: There were %zu errors!\n", rank, errors);
        for (size_t i=0; i<n; i++) {
            printf("%d: ptr[%zu] = %c\n", rank, i, ptr[i]);
        }
        fflush(stdout);
    }

    MPI_Type_free(&bigtype);

    MPI_Free_mem(ptr);

    MPI_Finalize();

    return 0;
}
