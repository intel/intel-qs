#ifndef BIGMPI_H_INCLUDED
#define BIGMPI_H_INCLUDED

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/* MPI_Count does not exist in MPI-2.  Our implementation
 * does not require it and in any case uses MPI_Aint in
 * place of MPI_Count in many places. */
#if MPI_VERSION < 3
typedef MPI_Aint MPI_Count;
#endif

/* MPI-3 added const to input arguments, which causes
 * incompatibilities if BigMPI passes in const arguments. */
#if MPI_VERSION >= 3
#define BIGMPI_CONST const
#else
#define BIGMPI_CONST
#endif

/* This function does the heavy lifting in BigMPI. */

int BigMPI_Type_contiguous(MPI_Aint offset, MPI_Count count, MPI_Datatype oldtype, MPI_Datatype * newtype);

/* other datatype creation functions */

int MPIX_Type_contiguous_x(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype * newtype);

/* rename argument "count" to "n" since global search-and-replace on "int count" happens somewhat often. */
int MPIX_Type_create_hvector_x(int n,
                               MPI_Count array_of_blocklengths[],
                               MPI_Aint array_of_displacements[],
                               MPI_Datatype oldtype,
                               MPI_Datatype * newtype);

/* These functions are primarily for internal use but some users may want to use them
 * so they will be in the public API, albeit with a different namespace. */

int BigMPI_Decode_contiguous_x(MPI_Datatype intype, MPI_Count * count, MPI_Datatype * basetype);

/* Requires distributed graph communicators. */
#if MPI_VERSION >= 3
int BigMPI_Create_graph_comm(MPI_Comm comm_old, int root, MPI_Comm * comm_dist_graph);
#endif

/* This is used in tests to query the compile-time setting. */
MPI_Count BigMPI_Get_max_int(void);

/* All of these functions should just be calling MPIX_Type_contiguous_x and
 * then the associated MPI function with count=1 and the newtype if the count
 * is bigger than INT_MAX and dropping into the regular implementation otherwise. */

/* Point-to-point */

int MPIX_Send_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPIX_Recv_x(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

int MPIX_Isend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPIX_Irecv_x(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);

int MPIX_Sendrecv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag,
                    MPI_Comm comm, MPI_Status *status);
int MPIX_Sendrecv_replace_x(void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int sendtag,
                            int source, int recvtag, MPI_Comm comm, MPI_Status *status);

int MPIX_Ssend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPIX_Rsend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPIX_Issend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                  MPI_Request *request);
int MPIX_Irsend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                  MPI_Request *request);

#if MPI_VERSION >= 3
int MPIX_Mrecv_x(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message, MPI_Status *status);
int MPIX_Imrecv_x(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message, MPI_Request *request);
#endif

/* Collectives */

int MPIX_Bcast_x(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm);
int MPIX_Gather_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                  void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPIX_Scatter_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                   void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPIX_Allgather_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                     void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPIX_Alltoall_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm);

#if MPI_VERSION >= 3
int MPIX_Ibcast_x(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request);
int MPIX_Igather_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                   void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
int MPIX_Iscatter_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
int MPIX_Iallgather_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
int MPIX_Ialltoall_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                     void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);

/* Neighborhood collectives */

int MPIX_Neighbor_allgather_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPIX_Neighbor_alltoall_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPIX_Ineighbor_allgather_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Request *request);
int MPIX_Ineighbor_alltoall_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request);
#endif

/* Reductions */

int MPIX_Reduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                  MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int MPIX_Allreduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPIX_Reduce_scatter_block_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
#if MPI_VERSION >= 3
int MPIX_Ireduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request);
int MPIX_Iallreduce_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count,
                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
int MPIX_Ireduce_scatter_block_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
#endif

/* RMA */

int MPIX_Accumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                      MPI_Op op, MPI_Win win);
int MPIX_Get_x(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win);
int MPIX_Put_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win);
#if MPI_VERSION >= 3
int MPIX_Get_accumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                          void *result_addr, MPI_Count result_count, MPI_Datatype result_datatype,
                          int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                          MPI_Op op, MPI_Win win);
int MPIX_Rput_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                MPI_Win win, MPI_Request *request);
int MPIX_Rget_x(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                MPI_Win win, MPI_Request *request);
int MPIX_Raccumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                       int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                       MPI_Op op, MPI_Win win, MPI_Request *request);
int MPIX_Rget_accumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                           void *result_addr, MPI_Count result_count, MPI_Datatype result_datatype,
                           int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                           MPI_Op op, MPI_Win win, MPI_Request *request);
#endif

/* V-collectives */

int MPIX_Gatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                   void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                   int root, MPI_Comm comm);
int MPIX_Scatterv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                    int root, MPI_Comm comm);
int MPIX_Allgatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                      MPI_Comm comm);
int MPIX_Alltoallv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                     void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], MPI_Datatype recvtype,
                     MPI_Comm comm);
int MPIX_Alltoallw_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                     void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                     MPI_Comm comm);

int MPIX_Neighbor_allgatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm);
int MPIX_Neighbor_alltoallv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                              MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                              const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm);
int MPIX_Neighbor_alltoallw_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                              const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                              const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);

#if 0

/* UNSUPPORTED */

/* These are equivalent to reduce+scatterv and will likely be supported in the future. */

int MPIX_Reduce_scatter_x(BIGMPI_CONST void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPIX_Ireduce_scatter_x(BIGMPI_CONST void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                           MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);

/* It is unclear if supporting these is worthwhile... */

int MPIX_Scan_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPIX_Exscan_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPIX_Iscan_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                 MPI_Request *request);
int MPIX_Iexscan_x(BIGMPI_CONST void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                   MPI_Request *request);

/* Nonblocking V-collectives */

/* These are really hard, if not impossible to support, because the argument vectors
 * have to be duplicated for use in Neighborhood_(i)alltoallw, and these vectors cannot
 * be deallocated until the operation has completed, but there is no mechanism to do that. */
/* TODO Implement these using generalized requests (and a thread)... */

int MPIX_Igatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                    void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                    int root, MPI_Comm comm, MPI_Request *request);
int MPIX_Iscatterv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
                     void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
int MPIX_Iallgatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                       MPI_Comm comm, MPI_Request *request);
int MPIX_Ialltoallv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], MPI_Datatype recvtype,
                      MPI_Comm comm, MPI_Request *request);

int MPIX_Ialltoallw_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                      MPI_Comm comm, MPI_Request *request);

int MPIX_Ineighbor_allgatherv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
int MPIX_Ineighbor_alltoallv_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                               MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                               const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
int MPIX_Ineighbor_alltoallw_x(BIGMPI_CONST void *sendbuf, const MPI_Count sendcounts[],
                               const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                               void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                               const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);

#endif // UNSUPPORTED

int MPIX_File_read_at_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_read_at_all_x(MPI_File fh, MPI_Offset offset, void * buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_write_at_x(MPI_File fh, MPI_Offset offset, const void * buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_write_at_all_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_iread_at_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request);
int MPIX_File_iwrite_at_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request);
int MPIX_File_read_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_read_all_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_write_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_write_all_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_iread_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request);
int MPIX_File_iwrite_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request);
int MPIX_File_read_shared_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_write_shared_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_iread_shared_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request);
int MPIX_File_iwrite_shared_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request);
int MPIX_File_read_ordered_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_write_ordered_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
int MPIX_File_read_at_all_begin_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype);
int MPIX_File_write_at_all_begin_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype);
int MPIX_File_read_all_begin_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype);
int MPIX_File_write_all_begin_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype);
int MPIX_File_read_ordered_begin_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype);
int MPIX_File_write_ordered_begin_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype);
int MPIX_File_iread_at_all_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
int MPIX_File_iwrite_at_all_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
int MPIX_File_iread_all_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
int MPIX_File_iwrite_all_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request);


#ifdef __cplusplus
}
#endif

#endif // BIGMPI_H_INCLUDED
