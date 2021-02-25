##------------------------------------------------------------------------------
## Copyright 2021 Intel Corporation
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

# Set timezone persistent
#ENV TZ=Europe/Berlin

# Set apt-get non-interactive for build
ARG DEBIAN_FRONTEND=noninteractive

# Fetch and install the GNU Make utility.
RUN apt-get update && apt-get install -y build-essential g++ make

# Fetch and install a generic MPI implementation.
RUN apt-get update && DEBIAN_FRONTEND=nonitneractive apt-get install -y mpich

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

# Setup the local build environment for the simulation framework.
WORKDIR /root/intelqs
# Copy from docker host cwd everything (the git project files) into the container
COPY . /root/intelqs

# ------------------------------------------------------------------
# If desired, a new user can be created in addition to 'root'.
# Uncomment lines below to create a new user named 'tester':
#RUN useradd --home-dir /home/tester --create-home tester
#WORKDIR /home/tester/intelqs
#COPY . /home/tester/intelqs

# Install Intel Quantum Simulator
RUN /bin/bash -c "source /opt/intel/mkl/bin/mklvars.sh intel64 ilp64"
RUN /bin/bash -c "mkdir build; cd build; CXX=g++ cmake -DIqsMPI=ON -DBuildExamples=ON -DIqsUtest=ON -DIqsPython=OFF .."
WORKDIR /root/intelqs/build
RUN make
WORKDIR /root/intelqs

LABEL mode="MPI" version="1.0" description="intel-qs built with MPI, no py interface, with Examples"

# Install lib for missing pthread module [necessary?]
RUN apt-get -y install libboost-all-dev

# Install and configure conda env
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p $HOME/miniconda
ENV PATH="/root/miniconda/bin:$PATH"
RUN /bin/bash -c ". ~/.bashrc && \
	 	conda install -y notebook && \
		conda install -y pybind11 && \
		conda install -y numpy && \
		conda install -y matplotlib"

# Dockerfile Ends here