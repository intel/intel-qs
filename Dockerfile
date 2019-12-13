##------------------------------------------------------------------------------
## Copyright 2019 Intel Corporation
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------

# Use an official Ubuntu linux image as the run time environment.
FROM ubuntu:latest

# Setup the proxy server environment.
ENV http_proxy http://proxy-chain.intel.com:911
ENV https_proxy http://proxy-chain.intel.com:911

# Fetch and install the GNU Make utility.
RUN apt-get update && apt-get install -y build-essential g++ make

# Fetch and install a generic MPI implementation.
RUN apt-get update && apt-get install -y mpich

# Fetch and install OpenSSH (client/server) for interacting between
# nodes of the cluster in a Docker swarm configuration.
RUN apt-get update && apt-get install -y openssh-client
RUN apt-get update && apt-get install -y openssh-server

# Fetch and install CMake 3.15 as required by Intel-QS build process. 
WORKDIR swpkgs/cmake3.15
RUN wget "https://github.com/Kitware/CMake/releases/download/v3.15.2/cmake-3.15.2-Linux-x86_64.tar.gz" 
RUN tar -xzf cmake-3.15.2-Linux-x86_64.tar.gz -C /usr/local/ --strip-components=1 

# Fetch and install the Intel MKL libraries required for building the Intel-QS simulator.
WORKDIR swpkgs/mkl
RUN wget "https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB"
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
RUN sh -c 'echo deb https://apt.repos.intel.com/mpi all main > /etc/apt/sources.list.d/intel-mpi.list'
RUN sh -c 'echo deb https://apt.repos.intel.com/tbb all main > /etc/apt/sources.list.d/intel-tbb.list'
RUN apt-get update
RUN apt-get install -y intel-mkl-64bit-2019.2-057
RUN apt-get install -y intel-mpi-2019.2-057
RUN apt-get install -y intel-tbb-2019.2-057
# Set the (global) environment variable MKLROOT to facilitate the build process.
RUN echo 'export MKLROOT="/opt/intel/mkl"' >> /etc/bash.bashrc
RUN export MKLROOT="/opt/intel/mkl"

# Install libraries for OpenMP.
RUN apt-get install -y libomp-dev

# Install GIT (needed to install the dependency of Googletest).
WORKDIR swpkgs/git
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

# Install VIM to edit files.
RUN apt-get install -y vim

# Setup the local build environment for the simulation framework.
WORKDIR /root/intelqs
COPY . /root/intelqs

#FIXME
# For stability reasons, the container should not be used as 'root'.
# A new user is created, named 'tester':
#RUN useradd --home-dir /home/tester --create-home tester
#WORKDIR /home/tester/intelqs
#COPY . /home/tester/intelqs

#FIXME
# Install Intel Quantum Simulator
#RUN source /opt/intel/mkl/bin/mklvars.sh intel64 ilp64
#RUN mkdir build; cd build; CXX=g++ cmake -DIqsMPI=OFF -DIqsUtest=ON -DIqsPython=OFF ..

