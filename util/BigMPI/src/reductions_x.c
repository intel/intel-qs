#include "bigmpi_impl.h"

/* There are different ways to implement large-count reductions.
 * The fully general and most correct way to do it is with user-defined
 * reductions, which are required to do reductions on user-defined types,
 * and to parse the datatype inside of a user-defined operation.
 * However, this appear is likely to lead to vastly reduced performance.
 *
 * A less general and potentially unsafe way to implement large-count
 * reductions is to chop them up into multiple messages.
 * This may formally violate the MPI standard, but since BigMPI is
 * not part of the standard, we are going to do it in the name
 * of performance and implementation simplicity. */

#define PASTE_BIGMPI_REDUCE_OP(OP)                                                      \
void BigMPI_##OP##_x(void * invec, void * inoutvec, int * len, MPI_Datatype * bigtype)  \
{                                                                                       \
    /* We are reducing a single element of bigtype... */                                \
    assert(*len==1);                                                                    \
                                                                                        \
    MPI_Count count;                                                                    \
    MPI_Datatype basetype;                                                              \
    BigMPI_Decode_contiguous_x(*bigtype, &count, &basetype);                            \
                                                                                        \
    int c = (int)(count/bigmpi_int_max);                                                \
    int r = (int)(count%bigmpi_int_max);                                                \
                                                                                        \
    /* Can use typesize rather than extent here because all built-ins lack holes. */    \
    int typesize;                                                                       \
    MPI_Type_size(basetype, &typesize);                                                 \
    for (int i=0; i<c; i++) {                                                           \
        MPI_Reduce_local(invec+(size_t)i*bigmpi_int_max*(size_t)typesize,               \
                         inoutvec+i*bigmpi_int_max*(size_t)typesize,                    \
                         bigmpi_int_max, basetype, MPI_##OP);                           \
    }                                                                                   \
    MPI_Reduce_local(invec+(size_t)c*bigmpi_int_max*(size_t)typesize,                   \
                     inoutvec+c*bigmpi_int_max*(size_t)typesize,                        \
                     r, basetype, MPI_##OP);                                            \
                                                                                        \
    return;                                                                             \
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

int BigMPI_Op_create(MPI_Op op, MPI_Op * bigop)
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
        BigMPI_Error("BigMPI does not support this op.  Sorry. \n");
    }
    return MPI_Op_create(bigfn, commute, bigop);
}

int MPIX_Reduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                  MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
    if (likely (count <= bigmpi_int_max )) {
        return MPI_Reduce(sendbuf, recvbuf, (int)count, datatype, op, root, comm);
    } else {
#ifdef BIGMPI_CLEAVER
        int c = (int)(count/bigmpi_int_max);
        int r = (int)(count%bigmpi_int_max);
        int typesize;
        MPI_Type_size(datatype, &typesize);
        if (sendbuf==MPI_IN_PLACE) {
            int commrank;
            MPI_Comm_rank(comm, &commrank);
            for (ptrdiff_t i=0; i<c; i++) {
                MPI_Reduce(commrank==root ? MPI_IN_PLACE : recvbuf+i*bigmpi_int_max*(size_t)typesize,
                           recvbuf+i*bigmpi_int_max*(size_t)typesize,
                           bigmpi_int_max, datatype, op, root, comm);
            }
            MPI_Reduce(commrank==root ? MPI_IN_PLACE : recvbuf+c*bigmpi_int_max*(size_t)typesize,
                       recvbuf+c*bigmpi_int_max*(size_t)typesize,
                       r, datatype, op, root, comm);
        } else {
            for (ptrdiff_t i=0; i<c; i++) {
                MPI_Reduce(sendbuf+i*bigmpi_int_max*(size_t)typesize,
                           recvbuf+i*bigmpi_int_max*(size_t)typesize,
                           bigmpi_int_max, datatype, op, root, comm);
            }
            MPI_Reduce(sendbuf+c*bigmpi_int_max*(size_t)typesize,
                       recvbuf+c*bigmpi_int_max*(size_t)typesize,
                       r, datatype, op, root, comm);
        }
        return MPI_SUCCESS;
#else /* BIGMPI_CLEAVER */

        MPI_Datatype bigtype;
        BigMPI_Type_contiguous(0,count, datatype, &bigtype);
        MPI_Type_commit(&bigtype);

        MPI_Op bigop;
        BigMPI_Op_create(op, &bigop);

        /* This is overkill.  Only the root needs to copy the buffer. */
        void * tempbuf = NULL;
        if (sendbuf==MPI_IN_PLACE) {
            MPI_Aint lb /* unused */, extent;
            MPI_Type_get_extent(datatype, &lb, &extent);
            MPI_Aint buf_size = (MPI_Aint)count * extent;

            MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &tempbuf);
            assert(tempbuf!=NULL);
            memcpy(tempbuf, recvbuf, (size_t)buf_size);
        }

        int rc = MPI_Reduce(sendbuf==MPI_IN_PLACE ? tempbuf : sendbuf,
                            recvbuf, 1, bigtype, bigop, root, comm);

        if (sendbuf==MPI_IN_PLACE) {
            MPI_Free_mem(&tempbuf);
        }

        MPI_Type_free(&bigtype);
        MPI_Op_free(&bigop);

        return rc;

#endif /* BIGMPI_CLEAVER */
    }
}

int MPIX_Allreduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    if (likely (count <= bigmpi_int_max )) {
        return MPI_Allreduce(sendbuf, recvbuf, (int)count, datatype, op, comm);
    } else {
#ifdef BIGMPI_CLEAVER
        int c = (int)(count/bigmpi_int_max);
        int r = (int)(count%bigmpi_int_max);
        int typesize;
        MPI_Type_size(datatype, &typesize);
        if (sendbuf==MPI_IN_PLACE) {
            for (ptrdiff_t i=0; i<c; i++) {
                MPI_Allreduce(MPI_IN_PLACE, recvbuf+i*bigmpi_int_max*(size_t)typesize,
                              bigmpi_int_max, datatype, op, comm);
            }
            MPI_Allreduce(MPI_IN_PLACE, recvbuf+c*bigmpi_int_max*(size_t)typesize,
                          r, datatype, op, comm);
        } else {
            for (ptrdiff_t i=0; i<c; i++) {
                MPI_Allreduce(sendbuf+i*bigmpi_int_max*(size_t)typesize,
                              recvbuf+i*bigmpi_int_max*(size_t)typesize,
                              bigmpi_int_max, datatype, op, comm);
            }
            MPI_Allreduce(sendbuf+c*bigmpi_int_max*(size_t)typesize,
                          recvbuf+c*bigmpi_int_max*(size_t)typesize,
                          r, datatype, op, comm);
        }
        return MPI_SUCCESS;
#else /* BIGMPI_CLEAVER */

        MPI_Datatype bigtype;
        BigMPI_Type_contiguous(0,count, datatype, &bigtype);
        MPI_Type_commit(&bigtype);

        MPI_Op bigop;
        BigMPI_Op_create(op, &bigop);

        void * tempbuf = NULL;
        if (sendbuf==MPI_IN_PLACE) {
            MPI_Aint lb /* unused */, extent;
            MPI_Type_get_extent(datatype, &lb, &extent);
            MPI_Aint buf_size = (MPI_Aint)count * extent;

            MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &tempbuf);
            assert(tempbuf!=NULL);
            memcpy(tempbuf, recvbuf, (size_t)buf_size);
        }
        int rc = MPI_Allreduce(sendbuf==MPI_IN_PLACE ? tempbuf : sendbuf,
                               recvbuf, 1, bigtype, bigop, comm);
        if (sendbuf==MPI_IN_PLACE) {
            MPI_Free_mem(&tempbuf);
        }

        MPI_Type_free(&bigtype);
        MPI_Op_free(&bigop);

        return rc;

#endif /* BIGMPI_CLEAVER */
    }
}

/* MPI-3 Section 5.10
 * Advice to implementers:
 * The MPI_REDUCE_SCATTER_BLOCK routine is functionally equivalent to:
 * an MPI_REDUCE collective operation with count equal to recvcount*n,
 * followed by an MPI_SCATTER with sendcount equal to recvcount. */

/* The previous statement is untrue when sendbuf=MPI_IN_PLACE so we
 * are forced to buffer even in the in-place case. */

int MPIX_Reduce_scatter_block_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    if (likely (recvcount <= bigmpi_int_max )) {
        return MPI_Reduce_scatter_block(sendbuf, recvbuf, (int)recvcount, datatype, op, comm);
    } else {
        int root = 0;

        int commsize;
        MPI_Comm_size(comm, &commsize);
        MPI_Count sendcount = recvcount * commsize;

        MPI_Aint lb /* unused */, extent;
        MPI_Type_get_extent(datatype, &lb, &extent);
        MPI_Aint buf_size = (MPI_Aint)sendcount * extent;

        void * tempbuf = NULL;
        MPI_Alloc_mem(buf_size, MPI_INFO_NULL, &tempbuf);
        assert(tempbuf!=NULL);

        MPIX_Reduce_x(sendbuf==MPI_IN_PLACE ? recvbuf : sendbuf,
                      tempbuf, sendcount, datatype, op, root, comm);
        MPIX_Scatter_x(tempbuf, recvcount, datatype, recvbuf, recvcount, datatype, root, comm);

        MPI_Free_mem(&tempbuf);
    }
    return MPI_SUCCESS;
}

#if MPI_VERSION >= 3

int MPIX_Ireduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request * req)
{
    if (likely (count <= bigmpi_int_max )) {
        return MPI_Ireduce(sendbuf, recvbuf, (int)count, datatype, op, root, comm, req);
    } else {

        MPI_Datatype bigtype;
        BigMPI_Type_contiguous(0,count, datatype, &bigtype);
        MPI_Type_commit(&bigtype);

        MPI_Op bigop;
        BigMPI_Op_create(op, &bigop);

        if (sendbuf==MPI_IN_PLACE)
            BigMPI_Error("BigMPI does not support in-place here yet.  Sorry. \n");

        int rc = MPI_Ireduce(sendbuf, recvbuf, 1, bigtype, bigop, root, comm, req);

        MPI_Type_free(&bigtype);
        MPI_Op_free(&bigop);

        return rc;

    }
}

int MPIX_Iallreduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request * req)
{
    if (likely (count <= bigmpi_int_max )) {
        return MPI_Iallreduce(sendbuf, recvbuf, (int)count, datatype, op, comm, req);
    } else {

        MPI_Datatype bigtype;
        BigMPI_Type_contiguous(0,count, datatype, &bigtype);
        MPI_Type_commit(&bigtype);

        MPI_Op bigop;
        BigMPI_Op_create(op, &bigop);

        if (sendbuf==MPI_IN_PLACE)
            BigMPI_Error("BigMPI does not support in-place here yet.  Sorry. \n");

        int rc = MPI_Iallreduce(sendbuf, recvbuf, 1, bigtype, bigop, comm, req);

        MPI_Type_free(&bigtype);
        MPI_Op_free(&bigop);

        return rc;

    }
}

#endif
