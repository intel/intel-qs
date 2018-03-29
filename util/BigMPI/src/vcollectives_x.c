#include "bigmpi_impl.h"
#include "pthread.h"

/* The displacements vector cannot be represented in the existing set of MPI-3
   functions because it is an integer rather than an MPI_Aint. */

typedef enum { GATHERV,
               SCATTERV,
               ALLGATHERV,
               ALLTOALLV,
               ALLTOALLW } bigmpi_collective_t;

typedef enum { ALLTOALLW_OFFSET,
#if MPI_VERSION >= 3
               NEIGHBORHOOD_ALLTOALLW,
#endif
               P2P,
               RMA } bigmpi_method_t;

static pthread_once_t BigMPI_vcollectives_method_is_initialized = PTHREAD_ONCE_INIT;
bigmpi_method_t BigMPI_vcollectives_method;

/* Tries to deduce collective operation implementation strategy from
   environment */
void BigMPI_Detect_default_vcollectives_method()
{
    char *env_var = getenv("BIGMPI_DEFAULT_METHOD");

    if (env_var != NULL) {
#if MPI_VERSION >= 3
        if (strcmp(env_var, "NEIGHBORHOOD_ALLTOALLW")) {
            BigMPI_vcollectives_method = NEIGHBORHOOD_ALLTOALLW;
            return;
        }
#endif

        if (strcmp(env_var, "P2P")) {
            BigMPI_vcollectives_method = P2P;
            return;
        }

        if (strcmp(env_var, "RMA")) {
            BigMPI_vcollectives_method = RMA;
            return;
        }

        fprintf(stderr, "Unknown value \"%s\" for environment variable BIGMPI_DEFAULT_METHOD\n", env_var);
    }

    // fallback to default:
    BigMPI_vcollectives_method = P2P;
}

bigmpi_method_t BigMPI_Get_default_vcollectives_method()
{
    pthread_once(&BigMPI_vcollectives_method_is_initialized, BigMPI_Detect_default_vcollectives_method);
    return BigMPI_vcollectives_method;
}

int BigMPI_Collective(bigmpi_collective_t coll, bigmpi_method_t method,
                      BIGMPI_CONST void *sendbuf,
                      const MPI_Count sendcount, const MPI_Count sendcounts[],
                      const MPI_Aint senddispls[],
                      const MPI_Datatype sendtype, const MPI_Datatype sendtypes[],
                      void *recvbuf,
                      const MPI_Count recvcount, const MPI_Count recvcounts[],
                      const MPI_Aint recvdispls[],
                      const MPI_Datatype recvtype, const MPI_Datatype recvtypes[],
                      int root,
                      MPI_Comm comm)
{
    int rc = MPI_SUCCESS;

    int is_intercomm;
    MPI_Comm_test_inter(comm, &is_intercomm);
    if (is_intercomm)
        BigMPI_Error("BigMPI does not support intercommunicators yet.\n");

    if (sendbuf==MPI_IN_PLACE)
        BigMPI_Error("BigMPI does not support in-place in the v-collectives.  Sorry. \n");

    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (method==P2P) {

        switch(coll) {
            case ALLTOALLW: /* See page 173 of MPI-3 */
                {
                    MPI_Request * reqs = malloc(2*size*sizeof(MPI_Request)); assert(reqs!=NULL);
                    /* No extent calculation because alltoallw does not use that. */
                    /* Use tag=0 because there is perfect pair-wise matching. */
                    for (int i=0; i<size; i++) {
                        /* Pre-post all receives... */
                        MPIX_Irecv_x(recvbuf+recvdispls[i], recvcounts[i], recvtypes[i],
                                     i /* source */, 0 /* tag */, comm, &reqs[i]);
                    }
                    for (int j=rank; j<(size+rank); j++) {
                        /* Schedule communication in balanced way... */
                        int i = j%size;
                        MPIX_Isend_x(sendbuf+senddispls[i], sendcounts[i], sendtypes[i],
                                     i /* target */, 0 /* tag */, comm, &reqs[size+i]);
                    }
                    MPI_Waitall(2*size, reqs, MPI_STATUSES_IGNORE);
                    free(reqs);
                }
                break;
            case ALLTOALLV: /* See page 171 of MPI-3 */
                {
                    MPI_Request * reqs = malloc(2*size*sizeof(MPI_Request)); assert(reqs!=NULL);
                    MPI_Aint lb /* unused */, sendextent, recvextent;
                    MPI_Type_get_extent(sendtype, &lb, &sendextent);
                    MPI_Type_get_extent(recvtype, &lb, &recvextent);
                    /* Use tag=0 because there is perfect pair-wise matching without it. */
                    for (int i=0; i<size; i++) {
                        /* Pre-post all receives... */
                        MPIX_Irecv_x(recvbuf+recvdispls[i]*recvextent, recvcounts[i], recvtype,
                                     i /* source */, 0 /* tag */, comm, &reqs[i]);
                    }
                    for (int j=rank; j<(size+rank); j++) {
                        /* Schedule communication in balanced way... */
                        int i = j%size;
                        MPIX_Isend_x(sendbuf+senddispls[i]*sendextent, sendcounts[i], sendtype,
                                     i /* target */, 0 /* tag */, comm, &reqs[size+i]);
                    }
                    MPI_Waitall(2*size, reqs, MPI_STATUSES_IGNORE);
                    free(reqs);
                }
                break;
            case ALLGATHERV:
                {
                    MPI_Request * reqs = malloc(2*size*sizeof(MPI_Request)); assert(reqs!=NULL);
                    MPI_Aint lb /* unused */, recvextent;
                    MPI_Type_get_extent(recvtype, &lb, &recvextent);
                    /* Use tag=0 because there is perfect pair-wise matching without it. */
                    for (int i=0; i<size; i++) {
                        /* Pre-post all receives... */
                        MPIX_Irecv_x(recvbuf+recvdispls[i]*recvextent, recvcounts[i], recvtype,
                                     i /* source */, 0 /* tag */, comm, &reqs[i]);
                    }
                    for (int j=rank; j<(size+rank); j++) {
                        /* Schedule communication in balanced way... */
                        int i = j%size;
                        MPIX_Isend_x(sendbuf, sendcount, sendtype,
                                     i /* target */, 0 /* tag */, comm, &reqs[size+i]);
                    }
                    MPI_Waitall(2*size, reqs, MPI_STATUSES_IGNORE);
                    free(reqs);
                }
                break;
            case GATHERV:
                {
                    int nreqs = (rank==root ? size+1 : 1);
                    MPI_Request * reqs = malloc(nreqs*sizeof(MPI_Request)); assert(reqs!=NULL);
                    if (rank==root) {
                        MPI_Aint lb /* unused */, recvextent;
                        MPI_Type_get_extent(recvtype, &lb, &recvextent);
                        for (int i=0; i<size; i++) {
                            /* Use tag=0 because there is perfect pair-wise matching without it. */
                            MPIX_Irecv_x(recvbuf+recvdispls[i]*recvextent, recvcounts[i], recvtype,
                                         i /* source */, 0 /* tag */, comm, &reqs[i+1]);
                        }
                    }
                    MPIX_Isend_x(sendbuf, sendcount, sendtype,
                                 root /* target */, 0 /* tag */, comm, &reqs[0]);
                    MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
                    free(reqs);
                }
                break;
            case SCATTERV:
                {
                    int nreqs = (rank==root ? size+1 : 1);
                    MPI_Request * reqs = malloc(nreqs*sizeof(MPI_Request)); assert(reqs!=NULL);
                    if (rank==root) {
                        MPI_Aint lb /* unused */, sendextent;
                        MPI_Type_get_extent(sendtype, &lb, &sendextent);
                        for (int i=0; i<size; i++) {
                            /* Use tag=0 because there is perfect pair-wise matching without it. */
                            MPIX_Isend_x(sendbuf+senddispls[i]*sendextent, sendcounts[i], sendtype,
                                         i /* target */, 0 /* tag */, comm, &reqs[i+1]);
                        }
                    }
                    MPIX_Irecv_x(recvbuf, recvcount, recvtype,
                                 root /* source */, 0 /* tag */, comm, &reqs[0]);
                    MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
                    free(reqs);
                }
                break;
            default:
                BigMPI_Error("Invalid collective chosen. \n");
                break;
        }
#if 0
    } else if (method==ALLTOALLW_OFFSET) {

        BigMPI_Error("ALLTOALL_OFFSET implementation of v-collectives is incomplete!\n");

        int          * newsendcounts = malloc(size*sizeof(int));          assert(newsendcounts!=NULL);
        MPI_Datatype * newsendtypes  = malloc(size*sizeof(MPI_Datatype)); assert(newsendtypes!=NULL);
        int          * newsdispls    = malloc(size*sizeof(int));          assert(newsdispls!=NULL);

        int          * newrecvcounts = malloc(size*sizeof(int));          assert(newrecvcounts!=NULL);
        MPI_Datatype * newrecvtypes  = malloc(size*sizeof(MPI_Datatype)); assert(newrecvtypes!=NULL);
        int          * newrdispls    = malloc(size*sizeof(int));          assert(newrdispls!=NULL);

        switch(coll) {
            case ALLTOALLW:
                assert(root == -1);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, sendcounts,
                                       0 /* splat type */, 0, sendtypes,
                                       0 /* zero displs */, senddispls,
                                       newsendcounts, newsendtypes, newsdispls);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, recvcounts,
                                       0 /* splat type */, 0, recvtypes,
                                       0 /* zero displs */, recvdispls,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            case ALLTOALLV:
                assert(root == -1);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, sendcounts,
                                       1 /* splat type */, sendtype, NULL,
                                       0 /* zero displs */, senddispls,
                                       newsendcounts, newsendtypes, newsdispls);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, recvcounts,
                                       1 /* splat type */, recvtype, NULL,
                                       0 /* zero displs */, recvdispls,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            case ALLGATHERV:
                assert(root == -1);
                BigMPI_Convert_vectors(size,
                                       1 /* splat count */, sendcount, NULL,
                                       1 /* splat type */, sendtype, NULL,
                                       1 /* zero displs */, NULL,
                                       newsendcounts, newsendtypes, newsdispls);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, recvcounts,
                                       1 /* splat type */, recvtype, NULL,
                                       0 /* zero displs */, recvdispls,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            case GATHERV:
                assert(root != -1);
                BigMPI_Convert_vectors(size,
                                       1 /* splat count */, sendcount, NULL,
                                       1 /* splat type */, sendtype, NULL,
                                       1 /* zero displs */, NULL,
                                       newsendcounts, newsendtypes, newsdispls);
                /* Gatherv: Only the root receives data. */
                if (rank==root) {
                    BigMPI_Convert_vectors(size,
                                           0 /* splat count */, 0, recvcounts,
                                           1 /* splat type */, recvtype, NULL,
                                           0 /* zero displs */, recvdispls,
                                           newrecvcounts, newrecvtypes, newrdispls);
                } else {
                    BigMPI_Convert_vectors(size,
                                           1 /* splat count */, 0, NULL,
                                           1 /* splat type */, MPI_DATATYPE_NULL, NULL,
                                           1 /* zero displs */, NULL,
                                           newrecvcounts, newrecvtypes, newrdispls);
                }
                break;
            case SCATTERV:
                assert(root != -1);
                /* Scatterv: Only the root sends data. */
                if (rank==root) {
                    BigMPI_Convert_vectors(size,
                                           0 /* splat count */, 0, sendcounts,
                                           1 /* splat type */, sendtype, NULL,
                                           0 /* zero displs */, senddispls,
                                           newsendcounts, newsendtypes, newsdispls);
                } else {
                    BigMPI_Convert_vectors(size,
                                           1 /* splat count */, 0, NULL,
                                           1 /* splat type */, MPI_DATATYPE_NULL, NULL,
                                           1 /* zero displs */, NULL,
                                           newsendcounts, newsendtypes, newsdispls);
                }
                BigMPI_Convert_vectors(size,
                                       1 /* splat count */, recvcount, NULL,
                                       1 /* splat type */, recvtype, NULL,
                                       1 /* zero displs */, NULL,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            default:
                BigMPI_Error("Invalid collective chosen. \n");
                break;
        }

        rc = MPI_Alltoallw(sendbuf, newsendcounts, newsdispls, newsendtypes,
                           recvbuf, newrecvcounts, newrdispls, newrecvtypes, comm);

        for (int i=0; i<size; i++) {
            MPI_Type_free(&newsendtypes[i]);
            MPI_Type_free(&newrecvtypes[i]);
        }
        free(newsendcounts);
        free(newsdispls);
        free(newsendtypes);

        free(newrecvcounts);
        free(newrecvtypes);
        free(newrdispls);

#endif

#if MPI_VERSION >= 3
    } else if (method==NEIGHBORHOOD_ALLTOALLW) {

        int          * newsendcounts = malloc(size*sizeof(int));          assert(newsendcounts!=NULL);
        MPI_Datatype * newsendtypes  = malloc(size*sizeof(MPI_Datatype)); assert(newsendtypes!=NULL);
        MPI_Aint     * newsdispls    = malloc(size*sizeof(MPI_Aint));     assert(newsdispls!=NULL);

        int          * newrecvcounts = malloc(size*sizeof(int));          assert(newrecvcounts!=NULL);
        MPI_Datatype * newrecvtypes  = malloc(size*sizeof(MPI_Datatype)); assert(newrecvtypes!=NULL);
        MPI_Aint     * newrdispls    = malloc(size*sizeof(MPI_Aint));     assert(newrdispls!=NULL);

        switch(coll) {
            case ALLTOALLW:
                assert(root == -1);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, sendcounts,
                                       0 /* splat type */, 0, sendtypes,
                                       0 /* zero displs */, senddispls,
                                       newsendcounts, newsendtypes, newsdispls);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, recvcounts,
                                       0 /* splat type */, 0, recvtypes,
                                       0 /* zero displs */, recvdispls,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            case ALLTOALLV:
                assert(root == -1);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, sendcounts,
                                       1 /* splat type */, sendtype, NULL,
                                       0 /* zero displs */, senddispls,
                                       newsendcounts, newsendtypes, newsdispls);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, recvcounts,
                                       1 /* splat type */, recvtype, NULL,
                                       0 /* zero displs */, recvdispls,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            case ALLGATHERV:
                assert(root == -1);
                BigMPI_Convert_vectors(size,
                                       1 /* splat count */, sendcount, NULL,
                                       1 /* splat type */, sendtype, NULL,
                                       1 /* zero displs */, NULL,
                                       newsendcounts, newsendtypes, newsdispls);
                BigMPI_Convert_vectors(size,
                                       0 /* splat count */, 0, recvcounts,
                                       1 /* splat type */, recvtype, NULL,
                                       0 /* zero displs */, recvdispls,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            case GATHERV:
                assert(root != -1);
                BigMPI_Convert_vectors(size,
                                       1 /* splat count */, sendcount, NULL,
                                       1 /* splat type */, sendtype, NULL,
                                       1 /* zero displs */, NULL,
                                       newsendcounts, newsendtypes, newsdispls);
                /* Gatherv: Only the root receives data. */
                if (rank==root) {
                    BigMPI_Convert_vectors(size,
                                           0 /* splat count */, 0, recvcounts,
                                           1 /* splat type */, recvtype, NULL,
                                           0 /* zero displs */, recvdispls,
                                           newrecvcounts, newrecvtypes, newrdispls);
                } else {
                    BigMPI_Convert_vectors(size,
                                           1 /* splat count */, 0, NULL,
                                           1 /* splat type */, MPI_DATATYPE_NULL, NULL,
                                           1 /* zero displs */, NULL,
                                           newrecvcounts, newrecvtypes, newrdispls);
                }
                break;
            case SCATTERV:
                assert(root != -1);
                /* Scatterv: Only the root sends data. */
                if (rank==root) {
                    BigMPI_Convert_vectors(size,
                                           0 /* splat count */, 0, sendcounts,
                                           1 /* splat type */, sendtype, NULL,
                                           0 /* zero displs */, senddispls,
                                           newsendcounts, newsendtypes, newsdispls);
                } else {
                    BigMPI_Convert_vectors(size,
                                           1 /* splat count */, 0, NULL,
                                           1 /* splat type */, MPI_DATATYPE_NULL, NULL,
                                           1 /* zero displs */, NULL,
                                           newsendcounts, newsendtypes, newsdispls);
                }
                BigMPI_Convert_vectors(size,
                                       1 /* splat count */, recvcount, NULL,
                                       1 /* splat type */, recvtype, NULL,
                                       1 /* zero displs */, NULL,
                                       newrecvcounts, newrecvtypes, newrdispls);
                break;
            default:
                BigMPI_Error("Invalid collective chosen. \n");
                break;
        }

        MPI_Comm comm_dist_graph;
        BigMPI_Create_graph_comm(comm, root, &comm_dist_graph);
        rc = MPI_Neighbor_alltoallw(sendbuf, newsendcounts, newsdispls, newsendtypes,
                                    recvbuf, newrecvcounts, newrdispls, newrecvtypes, comm_dist_graph);
        MPI_Comm_free(&comm_dist_graph);

        for (int i=0; i<size; i++) {
            MPI_Type_free(&newsendtypes[i]);
            MPI_Type_free(&newrecvtypes[i]);
        }
        free(newsendcounts);
        free(newsdispls);
        free(newsendtypes);

        free(newrecvcounts);
        free(newrecvtypes);
        free(newrdispls);

#endif
    } else if (method==RMA) {

        /* TODO Add MPI_Win_create_dynamic version of this?
         *      This may entail an MPI_Allgather over base addresses,
         *      but MPI_Win_create has to do that, plus some, so it
         *      may work out to be faster. */

        printf("RMA implementation of v-collectives is incomplete!\n");
        MPI_Abort(comm, 1);

        /* In the RMA implementation, we will treat send as source (buf) and recv as target (win). */
        MPI_Win win;
        /* This is the most (?) conservative approach possible, and assumes that datatypes are
         * noncontiguous and potentially out-of-order. */
        MPI_Aint max_size = 0;
        for (int i=0; i<size; i++) {
            MPI_Aint lb /* unused */, extent;
            MPI_Type_get_extent(recvtypes[i], &lb, &extent);
            MPI_Aint offset = recvdispls[i]+recvcounts[i]*extent;
            max_size = ((offset > max_size) ? offset : max_size);
        }
        MPI_Win_create(recvbuf, max_size, 1, MPI_INFO_NULL, comm, &win);
        MPI_Win_fence(MPI_MODE_NOPRECEDE | MPI_MODE_NOSTORE, win);
        for (int i=0; i<size; i++) {
            MPI_Put(sendbuf+senddispls[i], sendcounts[i], sendtypes[i],
                    i, recvdispls[i], recvcounts[i], recvtypes[i], win);
        }
        MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, win);
        MPI_Win_free(&win);

    } else {
        /* This should be unreachable... */
        BigMPI_Error("Invalid method for v-collectives chosen. \n");
    }
    return rc;
}

int MPIX_Gatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                   void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
{
    bigmpi_method_t method = BigMPI_Get_default_vcollectives_method();
    return BigMPI_Collective(GATHERV, method,
                             sendbuf, sendcount, NULL, NULL, sendtype, NULL,
                             recvbuf, -1 /* recvcount */, recvcounts, rdispls, recvtype, NULL,
                             root, comm);
}

int MPIX_Allgatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], MPI_Datatype recvtype,
                      MPI_Comm comm)
{
    bigmpi_method_t method = BigMPI_Get_default_vcollectives_method();
    return BigMPI_Collective(ALLGATHERV, method,
                             sendbuf, sendcount, NULL, NULL, sendtype, NULL,
                             recvbuf, -1 /* recvcount */, recvcounts, rdispls, recvtype, NULL,
                             -1 /* root */, comm);
}

int MPIX_Scatterv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    bigmpi_method_t method = BigMPI_Get_default_vcollectives_method();
    return BigMPI_Collective(SCATTERV, method,
                             sendbuf, -1 /* sendcount */, sendcounts, sdispls, sendtype, NULL,
                             recvbuf, recvcount, NULL, NULL, recvtype, NULL,
                             root, comm);
}

int MPIX_Alltoallv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                     void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], MPI_Datatype recvtype,
                     MPI_Comm comm)
{
    bigmpi_method_t method = BigMPI_Get_default_vcollectives_method();
    return BigMPI_Collective(ALLTOALLV, method,
                             sendbuf, -1 /* sendcount */, sendcounts, sdispls, sendtype, NULL,
                             recvbuf, -1 /* recvcount */, recvcounts, rdispls, recvtype, NULL,
                             -1 /* root */, comm);
}

int MPIX_Alltoallw_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                     void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                     MPI_Comm comm)
{
    bigmpi_method_t method = BigMPI_Get_default_vcollectives_method();
    return BigMPI_Collective(ALLTOALLW, method,
                             sendbuf, -1 /* sendcount */, sendcounts, sdispls, MPI_DATATYPE_NULL, sendtypes,
                             recvbuf, -1 /* recvcount */, recvcounts, rdispls, MPI_DATATYPE_NULL, recvtypes,
                             -1 /* root */, comm);
}
