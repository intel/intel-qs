#include "bigmpi_impl.h"

/* This is a workaround for tests so that BIGMPI_MAX_INT is visible without header inclusion. */
MPI_Count BigMPI_Get_max_int(void)
{
#ifdef BIGMPI_MAX_INT
    return BIGMPI_MAX_INT;
#else
    return INT_MAX;
#endif
}

/* Raise an internal fatal BigMPI error.
 *
 * @param[in] file Current file name (__FILE__)
 * @param[in] line Current line numeber (__LINE__)
 * @param[in] func Current function name (__func__)
 * @param[in] msg  Message to be printed
 * @param[in] code Exit error code
 */
void BigMPI_Error_impl(const char *file, const int line, const char *func, const char *msg, ...)
{
    va_list ap;
    int  disp;
    char string[500];

    disp  = 0;
    va_start(ap, msg);
    disp += vsnprintf(string, 500, msg, ap);
    va_end(ap);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    fprintf(stderr, "[%d] BigMPI Internal error in %s (%s:%d)\n[%d] Message: %s\n",
            rank, func, file, line, rank, string);
    MPI_Abort(MPI_COMM_WORLD, 100);
}

/*
 * Synopsis
 *
 * void BigMPI_Convert_vectors(..)
 *
 *  Input Parameter
 *
 *  int          num                length of all vectors (unless splat true)
 *  int          splat_old_count    if non-zero, use oldcount instead of iterating over vector (v-to-w)
 *  MPI_Count    oldcount           single count (ignored if splat_old_count==0)
 *  MPI_Count    oldcounts          vector of counts
 *  int          splat_old_type     if non-zero, use oldtype instead of iterating over vector (v-to-w)
 *  MPI_Datatype oldtype            single type (MPI_DATATYPE_NULL if splat_old_type==0)
 *  MPI_Datatype oldtypes           vector of types (NULL if splat_old_type!=0)
 *  int          zero_new_displs    set the displacement to zero (scatter/gather)
 *  MPI_Aint     olddispls          vector of displacements (NULL if zero_new_displs!=0)
 *
 * Output Parameters
 *
 *  int          newcounts
 *  MPI_Datatype newtypes
 *  MPI_Aint     newdispls
 *
 */
void BigMPI_Convert_vectors(int                num,
                            int                splat_old_count,
                            const MPI_Count    oldcount,
                            const MPI_Count    oldcounts[],
                            int                splat_old_type,
                            const MPI_Datatype oldtype,
                            const MPI_Datatype oldtypes[],
                            int                zero_new_displs,
                            const MPI_Aint     olddispls[],
                                  int          newcounts[],
                                  MPI_Datatype newtypes[],
                                  MPI_Aint     newdispls[])
{
    assert(splat_old_count || (oldcounts!=NULL));
    assert(splat_old_type  || (oldtypes!=NULL));
    assert(zero_new_displs || (olddispls!=NULL));

    MPI_Aint lb /* unused */, oldextent;
    if (splat_old_type) {
        MPI_Type_get_extent(oldtype, &lb, &oldextent);
    } else {
        /* !splat_old_type implies ALLTOALLW, which implies no displacement zeroing. */
        assert(!zero_new_displs);
    }

    for (int i=0; i<num; i++) {
        /* counts */
        newcounts[i] = 1;

        /* types */
        BigMPI_Type_contiguous(0, splat_old_count ? oldcount : oldcounts[i], 
                                  splat_old_type  ? oldtype  : oldtypes[i], &newtypes[i]);
        MPI_Type_commit(&newtypes[i]);

        /* displacements */
        MPI_Aint newextent;
        /* If we are not splatting old type, it implies ALLTOALLW,
         * which does not scale the displacement by the type extent,
         * nor would we ever zero the displacements. */
        if (splat_old_type) {
            MPI_Type_get_extent(newtypes[i], &lb, &newextent);
            newdispls[i] = (zero_new_displs ? 0 : olddispls[i]*oldextent/newextent);
        } else {
            newdispls[i] = olddispls[i];
        }
    }
    return;
}

#if MPI_VERSION >= 3

/*
 * Synopsis
 *
 * int BigMPI_Create_graph_comm(MPI_Comm comm_old, int root, MPI_Comm * graph_comm)
 *
 *  Input Parameter
 *
 *   comm_old           MPI communicator from which to create a graph comm
 *   root               integer id of root.  if -1, create fully connected graph,
 *                      which is appropriate for the all___ collectives.
 *
 * Output Parameters
 *
 *   graph_comm         MPI topology communicator associated with input communicator
 *   rc                 returns the rc from the MPI graph comm create function.
 *
 */
int BigMPI_Create_graph_comm(MPI_Comm comm_old, int root, MPI_Comm * comm_dist_graph)
{
    int rank, size;
    MPI_Comm_rank(comm_old, &rank);
    MPI_Comm_size(comm_old, &size);

    /* in the all case (root == -1), every rank is a destination for every other rank;
     * otherwise, only the root is a destination. */
    int indegree  = (root == -1 || root==rank) ? size : 0;
    /* in the all case (root == -1), every rank is a source for every other rank;
     * otherwise, all non-root processes are the source for only one rank (the root). */
    int outdegree = (root == -1 || root==rank) ? size : 1;

    int * sources      = malloc(indegree*sizeof(int));  assert(sources!=NULL);
    int * destinations = malloc(outdegree*sizeof(int)); assert(destinations!=NULL);

    for (int i=0; i<indegree; i++) {
        sources[i]      = i;
    }
    for (int i=0; i<outdegree; i++) {
        destinations[i] = (root == -1 || root==rank) ? i : root;
    }

    int rc = MPI_Dist_graph_create_adjacent(comm_old,
                indegree,  sources,      indegree==0  ? MPI_WEIGHTS_EMPTY : MPI_UNWEIGHTED,
                outdegree, destinations, outdegree==0 ? MPI_WEIGHTS_EMPTY : MPI_UNWEIGHTED,
                MPI_INFO_NULL, 0 /* reorder */, comm_dist_graph);

    free(sources);
    free(destinations);

    return rc;
}

#endif
