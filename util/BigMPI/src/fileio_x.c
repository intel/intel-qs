#include "bigmpi_impl.h"

/* Note: MPIO_Request will become either an MPICH grequest or an MPI_Request. */

int MPIX_File_read_at_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_at(fh, offset, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_at(fh, offset, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_at_all_x(MPI_File fh, MPI_Offset offset, void * buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_at_all(fh, offset, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_at_all(fh, offset, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_at_all_begin_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_at_all_begin(fh, offset, buf, (int)count, datatype);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_at_all_begin(fh, offset, buf, 1, newtype);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_all_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_all(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_all(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_shared_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_shared(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_shared(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_ordered_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_ordered(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_ordered(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}


int MPIX_File_read_all_begin_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_all_begin(fh, buf, (int)count, datatype);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_all_begin(fh, buf, 1, newtype);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_read_ordered_begin_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_read_ordered_begin(fh, buf, (int)count, datatype);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_read_ordered_begin(fh, buf, 1, newtype);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iread_at_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iread_at(fh, offset, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iread_at(fh, offset, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iread_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iread(fh, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iread(fh, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iread_shared_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iread_shared(fh, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iread_shared(fh, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iread_at_all_x(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iread_at_all(fh, offset, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iread_at_all(fh, offset, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iread_all_x(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iread_all(fh, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iread_all(fh, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}


int MPIX_File_write_at_x(MPI_File fh, MPI_Offset offset, const void * buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_at(fh, offset, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_at(fh, offset, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_at_all_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_at_all(fh, offset, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_at_all(fh, offset, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_at_all_begin_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_at_all_begin(fh, offset, buf, (int)count, datatype);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_at_all_begin(fh, offset, buf, 1, newtype);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_all_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_all(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_all(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_shared_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_shared(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_shared(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_ordered_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_ordered(fh, buf, (int)count, datatype, status);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_ordered(fh, buf, 1, newtype, status);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_all_begin_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_all_begin(fh, buf, (int)count, datatype);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_all_begin(fh, buf, 1, newtype);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_write_ordered_begin_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_write_ordered_begin(fh, buf, (int)count, datatype);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_write_ordered_begin(fh, buf, 1, newtype);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iwrite_at_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iwrite_at(fh, offset, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iwrite_at(fh, offset, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iwrite_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iwrite(fh, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iwrite(fh, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iwrite_shared_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPIO_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iwrite_shared(fh, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iwrite_shared(fh, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iwrite_at_all_x(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iwrite_at_all(fh, offset, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iwrite_at_all(fh, offset, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}

int MPIX_File_iwrite_all_x(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (count <= bigmpi_int_max )) {
        rc = MPI_File_iwrite_all(fh, buf, (int)count, datatype, request);
    } else {
        MPI_Datatype newtype;
        BigMPI_Type_contiguous(0,count, datatype, &newtype);
        MPI_Type_commit(&newtype);
        rc = MPI_File_iwrite_all(fh, buf, 1, newtype, request);
        MPI_Type_free(&newtype);
    }
    return rc;
}
