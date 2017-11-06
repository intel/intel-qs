//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------
/**
 * @file mpi_test1.cpp
 *
 * This file tests creation/destruction of an MPI session.
 *
 */
#include "../mpi_wrapper.hpp"
#include <omp.h>
#include <stdio.h>
#include <stdexcept>
#include <cassert>


int main(int argc, char **argv) {

    int world_size = 0;
    int world_rank = 0;
    int my_node_size = 0;
    int my_node_rank = 0;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;

    // Initialize the MPI Framework.
    qhipster::MpiWrapper mpi(argc,argv);

    QH_MPI_STATUS_CHECK((MPI_Comm_size(MPI_COMM_WORLD, &world_size)));
    QH_MPI_STATUS_CHECK((MPI_Comm_rank(MPI_COMM_WORLD, &world_rank)));
    QH_MPI_STATUS_CHECK((MPI_Get_processor_name(processor_name, &name_len)));

    printf("World Rank/Size [%2d:%-2d]\n", world_rank, world_size);

    // Split the communicator into subcommunicators which have the capability of creating
    // shared memory.
    MPI_Comm node_comm;
    QH_MPI_STATUS_CHECK(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 
                        world_rank, MPI_INFO_NULL, &node_comm));

    // Determine the number of ranks in the new node we just split and also
    // my rank within it.
    QH_MPI_STATUS_CHECK(MPI_Comm_size(node_comm, &my_node_size));
    QH_MPI_STATUS_CHECK(MPI_Comm_rank(node_comm, &my_node_rank));

    printf("Group Rank/Node#: [%2d/%-2d]\n", my_node_rank, my_node_size / world_size);

    // This should detect non-symmetric partitioning of the processes into
    // uniform sized groups of nodes. For example: World: 16 processes,
    // 4 processes per comm group, 0 left-over processes.
    assert((world_size % my_node_size) == 0);

    // Partition the available global ranks into communication groups that are an even multiple of 2.
    // Left-over ranks stay in the global rank pool.
    





#pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        printf("Node ID/Rank: [%2d/%-2d] Thread [%d] working...\n", my_node_size/world_size, my_node_rank, threadID);

        int a=0,b=0,c=0;
        for(unsigned long i=0;i<2000000000;i++) {
            unsigned long a = (++b/2) + (++c/2);
        }
    }

    return 1;
}
