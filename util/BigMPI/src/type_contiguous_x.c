#include "bigmpi_impl.h"

/* This function does all the heavy lifting in BigMPI. */

#ifdef BIGMPI_AVOID_TYPE_CREATE_STRUCT
#include <math.h>
/*
 * Synopsis
 *
 * int BigMPI_Factorize_count(MPI_Count c, int * a, int *b)
 *
 *  Input Parameter
 *
 *   c                  large count
 *
 * Output Parameters
 *
 *   a, b               integers such that c=a*b and a,b<INT_MAX
 *   rc                 returns 0 if a,b found (success), else 1 (failure)
 *
 */
static int BigMPI_Factorize_count(MPI_Count in, int * a, int *b)
{
    /* THIS FUNCTION IS NOT OPTIMIZED AND MAY RUN VERY SLOWLY IN MANY CASES */
    /* TODO Implement something other than brute-force search for prime factors. */

    /* Is it better to do the division as MPI_Count (often long long) or double? */
    MPI_Count lo = in/bigmpi_int_max+1;
    MPI_Count hi = (MPI_Count)floor(sqrt((double)in));

    /* FIXME This is not safe.  Must test for overflow before casting to int. */
    for (MPI_Count g=hi; g>lo; g--) {
        MPI_Count rem = in%g;
        if (rem==0) {
            *a = (int)g;
            *b = (int)(in/g);
            return 0;
        }
    }
    *a = 1;
    *b = -1;
    return 1;
}
#endif

/*
 * Synopsis
 *
 * int BigMPI_Type_contiguous(MPI_Aint offset,
 *                            MPI_Count count,
 *                            MPI_Datatype   oldtype,
 *                            MPI_Datatype * newtype)
 *
 *  Input Parameters
 *
 *   offset            byte offset of the start of the contiguous chunk
 *   count             replication count (nonnegative integer)
 *   oldtype           old datatype (handle)
 *
 * Output Parameter
 *
 *   newtype           new datatype (handle)
 *
 * Notes
 *
 *   Following the addition of the offset argument, this function no longer
 *   matches the signature of MPI_Type_contiguous.  This may constitute
 *   breaking user experience for some people.  However, the value of
 *   adding it simplies the primary purpose of this function, which is to
 *   do the heavy lifting _inside_ of BigMPI.  In particular, it allows
 *   us to use MPI_Alltoallw instead of MPI_Neighborhood_alltoallw.
 *
 */
int BigMPI_Type_contiguous(MPI_Aint offset, MPI_Count count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    /* The count has to fit into MPI_Aint for BigMPI to work. */
    if ((uint64_t)count>(uint64_t)bigmpi_count_max) {
        printf("count (%llu) exceeds bigmpi_count_max (%llu)\n",
                (long long unsigned)count, (long long unsigned)bigmpi_count_max);
        fflush(stdout);
    }

#ifdef BIGMPI_AVOID_TYPE_CREATE_STRUCT
    if (offset==0) {
        /* There is no need for this code path in homogeneous execution,
         * but it is useful to exercise anyways. */
        int a, b;
        int prime = BigMPI_Factorize_count(count, &a, &b);
        if (!prime) {
            MPI_Type_vector(a, b, b, oldtype, newtype);
            return MPI_SUCCESS;
        }
    }
#endif
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
    MPI_Aint displacements[2] = {offset,offset+remdisp};
    MPI_Datatype types[2]     = {chunks,remainder};
    MPI_Type_create_struct(2, blocklengths, displacements, types, newtype);

    MPI_Type_free(&chunks);
    MPI_Type_free(&remainder);

    return MPI_SUCCESS;
}

/*
 * Synopsis
 *
 * This function inverts BigMPI_Type_contiguous, i.e. it provides
 * the original arguments for that call so that we know how many
 * built-in types are in the user-defined datatype.
 *
 * This function is primary used inside of BigMPI and does not
 * correspond to an MPI function, so we do avoid the use of the
 * MPIX namespace.
 *
 * int BigMPI_Decode_contiguous_x(MPI_Datatype   intype,
 *                                MPI_Count    * count,
 *                                MPI_Datatype * basetype)
 *
 *  Input Parameters
 *
 *   newtype           new datatype (handle)
 *
 * Output Parameter
 *
 *   count             replication count (nonnegative integer)
 *   oldtype           old datatype (handle)
 *
 */
int BigMPI_Decode_contiguous_x(MPI_Datatype intype, MPI_Count * count, MPI_Datatype * basetype)
{
    int nint, nadd, ndts, combiner;

    /* Step 1: Decode the type_create_struct call. */

    MPI_Type_get_envelope(intype, &nint, &nadd, &ndts, &combiner);
    assert(combiner==MPI_COMBINER_STRUCT || combiner==MPI_COMBINER_VECTOR);
#ifdef BIGMPI_AVOID_TYPE_CREATE_STRUCT
    if (combiner==MPI_COMBINER_VECTOR) {
        assert(nint==3);
        assert(nadd==0);
        assert(ndts==1);

        int cbs[3]; /* {count,blocklength,stride} */
        MPI_Datatype vbasetype[1];
        MPI_Type_get_contents(intype, 3, 0, 1, cbs, NULL, vbasetype);
        MPI_Count a = cbs[0];   /* count */
        MPI_Count b = cbs[1];   /* blocklength */
        assert(cbs[1]==cbs[2]); /* blocklength==stride */

        *count = a*b;
        *basetype = vbasetype[0];
        return MPI_SUCCESS;
    }
#else
    assert(combiner==MPI_COMBINER_STRUCT);
#endif
    assert(nint==3);
    assert(nadd==2);
    assert(ndts==2);

    int cnbls[3]; /* {count, blocklengths[]} */
    MPI_Aint displacements[2]; /* {0,remdisp} */
    MPI_Datatype types[2]; /* {chunks,remainder} */;
    MPI_Type_get_contents(intype, 3, 2, 2, cnbls, displacements, types);
    assert(cnbls[0]==2);
    assert(cnbls[1]==1);
    assert(cnbls[2]==1);
    assert(displacements[0]==0);

    /* Step 2: Decode the type_vector call. */

    MPI_Type_get_envelope(types[0], &nint, &nadd, &ndts, &combiner);
    assert(combiner==MPI_COMBINER_VECTOR);
    assert(nint==3);
    assert(nadd==0);
    assert(ndts==1);

    int cbs[3]; /* {count,blocklength,stride} */
    MPI_Datatype vbasetype[1];
    MPI_Type_get_contents(types[0], 3, 0, 1, cbs, NULL, vbasetype);
    assert(/* blocklength = */ cbs[1]==bigmpi_int_max);
    assert(/* stride = */ cbs[2]==bigmpi_int_max);

    /* chunk count - see above */
    MPI_Count c = cbs[0];

    /* Step 3: Decode the type_contiguous call. */

    MPI_Type_get_envelope(types[1], &nint, &nadd, &ndts, &combiner);
    assert(combiner==MPI_COMBINER_CONTIGUOUS);
    assert(nint==1);
    assert(nadd==0);
    assert(ndts==1);

    int ccc[1]; /* {count} */
    MPI_Datatype cbasetype[1];
    MPI_Type_get_contents(types[1], 1, 0, 1, ccc, NULL, cbasetype);

    /* remainder - see above */
    MPI_Count r = ccc[0];

    /* The underlying type of the vector and contig types must match. */
    assert(cbasetype[0]==vbasetype[0]);
    *basetype = cbasetype[0];

    /* This should not overflow because everything is already MPI_Count type. */
    *count = c*bigmpi_int_max+r;

    return MPI_SUCCESS;
}

/* MPIX_Type_contiguous_x is consistent with MPI_Type_contiguous... */
int MPIX_Type_contiguous_x(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    return BigMPI_Type_contiguous(0, count, oldtype, newtype);
}
