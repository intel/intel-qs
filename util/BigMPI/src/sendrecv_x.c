#include "bigmpi_impl.h"

int MPIX_Send_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Send(buf, (int)count, datatype, dest, tag, comm);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Send(buf, 1, newtype, dest, tag, comm);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Recv_x(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Recv(buf, (int)count, datatype, source, tag, comm, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Recv(buf, 1, newtype, source, tag, comm, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Isend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request * request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Isend(buf, (int)count, datatype, dest, tag, comm, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Isend(buf, 1, newtype, dest, tag, comm, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Irecv_x(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Irecv(buf, (int)count, datatype, source, tag, comm, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Irecv(buf, 1, newtype, source, tag, comm, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Sendrecv_x(BIGMPI_CONST void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag,
                    MPI_Comm comm, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (sendcount <= bigmpi_int_max && recvcount <= bigmpi_int_max )) {
        rc = MPI_Sendrecv(sendbuf, (int)sendcount, sendtype, dest, sendtag,
                          recvbuf, (int)recvcount, recvtype, source, recvtag,
                          comm, status);
    } else if (sendcount <= bigmpi_int_max && recvcount > bigmpi_int_max ) {
        MPI_Datatype newrecvtype;
        BigMPI_Type_contiguous(0,recvcount, recvtype, &newrecvtype);
        MPI_Type_commit(&newrecvtype);
        rc = MPI_Sendrecv(sendbuf, (int)sendcount, sendtype, dest, sendtag,
                          recvbuf, 1, newrecvtype, source, recvtag,
                          comm, status);
        MPI_Type_free(&newrecvtype);
    } else if (sendcount > bigmpi_int_max && recvcount <= bigmpi_int_max ) {
        MPI_Datatype newsendtype;
        BigMPI_Type_contiguous(0,sendcount, sendtype, &newsendtype);
        MPI_Type_commit(&newsendtype);
        rc = MPI_Sendrecv(sendbuf, 1, newsendtype, dest, sendtag,
                          recvbuf, (int)recvcount, recvtype, source, recvtag,
                          comm, status);
        MPI_Type_free(&newsendtype);
    } else {
        MPI_Datatype newsendtype, newrecvtype;
        BigMPI_Type_contiguous(0,sendcount, sendtype, &newsendtype);
        BigMPI_Type_contiguous(0,recvcount, recvtype, &newrecvtype);
        MPI_Type_commit(&newsendtype);
        MPI_Type_commit(&newrecvtype);
        rc = MPI_Sendrecv(sendbuf, 1, newsendtype, dest, sendtag,
                          recvbuf, 1, newrecvtype, source, recvtag,
                          comm, status);
        MPI_Type_free(&newsendtype);
        MPI_Type_free(&newrecvtype);
    }
    return rc;
}

int MPIX_Sendrecv_replace_x(void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int sendtag,
                            int source, int recvtag, MPI_Comm comm, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Sendrecv_replace(buf, (int)count, datatype, dest, sendtag, source, recvtag, comm, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Sendrecv_replace(buf, 1, newtype, dest, sendtag, source, recvtag, comm, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}


int MPIX_Ssend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Ssend(buf, (int)count, datatype, dest, tag, comm);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Ssend(buf, 1, newtype, dest, tag, comm);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Rsend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Rsend(buf, (int)count, datatype, dest, tag, comm);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Rsend(buf, 1, newtype, dest, tag, comm);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Issend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Issend(buf, (int)count, datatype, dest, tag, comm, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Issend(buf, 1, newtype, dest, tag, comm, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Irsend_x(BIGMPI_CONST void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Irsend(buf, (int)count, datatype, dest, tag, comm, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Irsend(buf, 1, newtype, dest, tag, comm, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

#if MPI_VERSION >= 3

int MPIX_Mrecv_x(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Mrecv(buf, (int)count, datatype, message, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Mrecv(buf, 1, newtype, message, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_Imrecv_x(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_Imrecv(buf, (int)count, datatype, message, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_Imrecv(buf, 1, newtype, message, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

#endif
