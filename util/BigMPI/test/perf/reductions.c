#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include <mpi.h>

static int verify_doubles(double *buf, int count, double expected_value)
{
    int errors = 0;
    for (int i = 0; i < count; i++) {
        double absdiff = fabs(buf[i] - expected_value);
        if (absdiff>1.e-4) errors++;
    }
    return errors;
}

#define PASTE_BIGMPI_REDUCE_OP(OP)                                                   \
void BigMPI_##OP##_x(void * invec, void * inoutvec, int * len, MPI_Datatype * type)  \
{                                                                                    \
    MPI_Reduce_local(invec,inoutvec,*len,*type,MPI_##OP);                            \
}

/* Create a BigMPI_<op>_x for all built-in ops. */
PASTE_BIGMPI_REDUCE_OP(MAX)
PASTE_BIGMPI_REDUCE_OP(MIN)
PASTE_BIGMPI_REDUCE_OP(SUM)
PASTE_BIGMPI_REDUCE_OP(PROD)
PASTE_BIGMPI_REDUCE_OP(LAND)
PASTE_BIGMPI_REDUCE_OP(BAND)
PASTE_BIGMPI_REDUCE_OP(LOR)
PASTE_BIGMPI_REDUCE_OP(BOR)
PASTE_BIGMPI_REDUCE_OP(LXOR)
PASTE_BIGMPI_REDUCE_OP(BXOR)
PASTE_BIGMPI_REDUCE_OP(MAXLOC)
PASTE_BIGMPI_REDUCE_OP(MINLOC)

#undef PASTE_BIGMPI_REDUCE_OP

int BigMPI_Op_create(MPI_Op op, MPI_Op * userop)
{
    int commute;
    MPI_Op_commutative(op, &commute);

    MPI_User_function * bigfn = NULL;

    if      (op==MPI_MAX)    bigfn = BigMPI_MAX_x;
    else if (op==MPI_MIN)    bigfn = BigMPI_MIN_x;
    else if (op==MPI_SUM)    bigfn = BigMPI_SUM_x;
    else if (op==MPI_PROD)   bigfn = BigMPI_PROD_x;
    else if (op==MPI_LAND)   bigfn = BigMPI_LAND_x;
    else if (op==MPI_BAND)   bigfn = BigMPI_BAND_x;
    else if (op==MPI_LOR)    bigfn = BigMPI_LOR_x;
    else if (op==MPI_BOR)    bigfn = BigMPI_BOR_x;
    else if (op==MPI_LXOR)   bigfn = BigMPI_LXOR_x;
    else if (op==MPI_BXOR)   bigfn = BigMPI_BXOR_x;
    else if (op==MPI_MAXLOC) bigfn = BigMPI_MAXLOC_x;
    else if (op==MPI_MINLOC) bigfn = BigMPI_MINLOC_x;
    else {
        fprintf(stderr,"BigMPI does not support this op.  Sorry. \n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    return MPI_Op_create(bigfn, commute, userop);
}

int BigMPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    MPI_Op userop;
    BigMPI_Op_create(op, &userop);
    MPI_Allreduce(sendbuf, recvbuf, count, datatype, userop, comm);
    MPI_Op_free(&userop);
    return MPI_SUCCESS;
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = (argc > 1) ? atoi(argv[1]) : 1000;

    double * sbuf1 = NULL;
    double * rbuf1 = NULL;
    double * sbuf2 = NULL;
    double * rbuf2 = NULL;

    MPI_Aint bytes = n*sizeof(double);
    MPI_Alloc_mem(bytes, MPI_INFO_NULL, &sbuf1);
    MPI_Alloc_mem(bytes, MPI_INFO_NULL, &rbuf1);
    MPI_Alloc_mem(bytes, MPI_INFO_NULL, &sbuf2);
    MPI_Alloc_mem(bytes, MPI_INFO_NULL, &rbuf2);

    for (int i=0; i<n; i++) {
        sbuf1[i] = (double)rank;
    }
    for (int i=0; i<n; i++) {
        rbuf1[i] = 0.0;
    }
    for (int i=0; i<n; i++) {
        sbuf2[i] = (double)rank;
    }
    for (int i=0; i<n; i++) {
        rbuf2[i] = 0.0;
    }

    /* Correctness checking first */

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(sbuf1, rbuf1, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    BigMPI_Allreduce(sbuf1, rbuf2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    const double val = (double)size*(size-1)/2.;
    int error1 = verify_doubles(rbuf1, n, val);
    if (error1>0) {
        printf("There were %d errors out of %d elements!\n", error1, n);
        for (int i=0; i<n; i++) {
            printf("rbuf1[%d] = %lf (expected %lf - %s)\n",
                    i, rbuf1[i], val, rbuf1[i]==val ? "RIGHT" : "WRONG");
        }
        fflush(stdout);
    }
    int error2 = verify_doubles(rbuf2, n, val);
    if (error2>0) {
        printf("There were %d errors out of %d elements!\n", error2, n);
        for (int i=0; i<n; i++) {
            printf("rbuf2[%d] = %lf (expected %lf - %s)\n",
                    i, rbuf2[i], val, rbuf2[i]==val ? "RIGHT" : "WRONG");
        }
        fflush(stdout);
    }

    /* collective communication */
    double t0, t1, dtmpi, dtusr;
    MPI_Barrier(MPI_COMM_WORLD);
    
    t0 = MPI_Wtime();
    MPI_Allreduce(sbuf1, rbuf1, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    dtmpi = t1-t0;
    MPI_Barrier(MPI_COMM_WORLD);
    
    t0 = MPI_Wtime();
    BigMPI_Allreduce(sbuf2, rbuf2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    dtusr = t1-t0;
    
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Free_mem(sbuf2);
    MPI_Free_mem(rbuf2);
    MPI_Free_mem(sbuf1);
    MPI_Free_mem(rbuf1);

    /* TODO: reduce errors across all ranks in case root result is correct
     *       but others are wrong. */
    if (rank==0 && error1==0 && error2==0) {
        printf("n = %d tmpi = %lf dtusr = %lf\n", n, dtmpi, dtusr);
    }

    MPI_Finalize();

    return 0;
}
