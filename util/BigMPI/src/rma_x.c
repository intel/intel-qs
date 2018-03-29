#include "bigmpi_impl.h"

int MPIX_Put_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Put(origin_addr, origin_count, origin_datatype,
                     target_rank, target_disp, target_count, target_datatype, win);
    } else {
        MPI_Datatype neworigin_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Put(origin_addr, 1, neworigin_datatype,
                     target_rank, target_disp, 1, newtarget_datatype, win);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

int MPIX_Get_x(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Get(origin_addr, origin_count, origin_datatype,
                     target_rank, target_disp, target_count, target_datatype, win);
    } else {
        MPI_Datatype neworigin_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Get(origin_addr, 1, neworigin_datatype,
                     target_rank, target_disp, 1, newtarget_datatype, win);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

int MPIX_Accumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                      MPI_Op op, MPI_Win win)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Accumulate(origin_addr, origin_count, origin_datatype,
                            target_rank, target_disp, target_count, target_datatype,
                            op, win);
    } else {
        MPI_Datatype neworigin_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Accumulate(origin_addr, 1, neworigin_datatype,
                            target_rank, target_disp, 1, newtarget_datatype, op, win);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

#if MPI_VERSION >= 3

int MPIX_Get_accumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                          void *result_addr, MPI_Count result_count, MPI_Datatype result_datatype,
                          int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                          MPI_Op op, MPI_Win win)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && result_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Get_accumulate(origin_addr, origin_count, origin_datatype,
                                result_addr, result_count, result_datatype,
                                target_rank, target_disp, target_count, target_datatype,
                                op, win);
    } else {
        MPI_Datatype neworigin_datatype, newresult_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,result_count, result_datatype, &newresult_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newresult_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Get_accumulate(origin_addr, 1, neworigin_datatype,
                                result_addr, 1, newresult_datatype,
                                target_rank, target_disp, 1, newtarget_datatype,
                                op, win);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newresult_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

int MPIX_Rput_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                MPI_Win win, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Rput(origin_addr, origin_count, origin_datatype,
                     target_rank, target_disp, target_count, target_datatype, win, request);
    } else {
        MPI_Datatype neworigin_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Rput(origin_addr, 1, neworigin_datatype,
                     target_rank, target_disp, 1, newtarget_datatype, win, request);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

int MPIX_Rget_x(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                MPI_Win win, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Rget(origin_addr, origin_count, origin_datatype,
                     target_rank, target_disp, target_count, target_datatype, win, request);
    } else {
        MPI_Datatype neworigin_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Rget(origin_addr, 1, neworigin_datatype,
                     target_rank, target_disp, 1, newtarget_datatype, win, request);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

int MPIX_Raccumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                      MPI_Op op, MPI_Win win, MPI_Request *request)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Raccumulate(origin_addr, origin_count, origin_datatype,
                            target_rank, target_disp, target_count, target_datatype,
                            op, win, request);
    } else {
        MPI_Datatype neworigin_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Raccumulate(origin_addr, 1, neworigin_datatype,
                            target_rank, target_disp, 1, newtarget_datatype, op, win, request);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

int MPIX_Rget_accumulate_x(BIGMPI_CONST void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                          void *result_addr, MPI_Count result_count, MPI_Datatype result_datatype,
                          int target_rank, MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                          MPI_Op op, MPI_Win win, MPI_Request * request)
{
    int rc = MPI_SUCCESS;

    if (likely (origin_count <= bigmpi_int_max && result_count <= bigmpi_int_max && target_count <= bigmpi_int_max)) {
        rc = MPI_Rget_accumulate(origin_addr, origin_count, origin_datatype,
                                result_addr, result_count, result_datatype,
                                target_rank, target_disp, target_count, target_datatype,
                                op, win, request);
    } else {
        MPI_Datatype neworigin_datatype, newresult_datatype, newtarget_datatype;
        BigMPI_Type_contiguous(0,origin_count, origin_datatype, &neworigin_datatype);
        BigMPI_Type_contiguous(0,result_count, result_datatype, &newresult_datatype);
        BigMPI_Type_contiguous(0,target_count, target_datatype, &newtarget_datatype);
        MPI_Type_commit(&neworigin_datatype);
        MPI_Type_commit(&newresult_datatype);
        MPI_Type_commit(&newtarget_datatype);
        rc = MPI_Rget_accumulate(origin_addr, 1, neworigin_datatype,
                                result_addr, 1, newresult_datatype,
                                target_rank, target_disp, 1, newtarget_datatype,
                                op, win, request);
        MPI_Type_free(&neworigin_datatype);
        MPI_Type_free(&newresult_datatype);
        MPI_Type_free(&newtarget_datatype);
    }
    return rc;
}

#endif
