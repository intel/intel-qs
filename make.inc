##------------------------------------------------------------------------------
## Copyright 2017 Intel Corporation 
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
## http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------

bigmpi = 0

HOST = $(shell hostname)

CC = g++
CC = CC -mkl -no-ipo
CC = mpicxx
CC = mpiicpc

OPTS = -g -O3
VECS = -qopt-report-phase=vec -qopt-report=2
DEFS = -DUSE_MM_MALLOC -DSTANDALONE -DOPENQU_HAVE_MPI -DMKL_ILP64
DEFS = -DUSE_MM_MALLOC -DSTANDALONE -DINTELQS_HAS_MPI -DMKL_ILP64

BIGMPI_LIB = /opt/bigmpi/lib/libbigmpi.a
BIGMPI_INC = -I/opt/bigmpi/include

CFLAGS = ${OPTS} ${DEFS} \
         -Wall -std=c++11 -qopenmp #-restrict \
# -fasm-blocks

CFLAGS += ${OPTFLAGS}

#ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
INC = -I./ -I../

ifeq ($(bigmpi),1)
 DEFS := ${DEFS} -DBIGMPI
 INC := ${INC} ${BIGMPI_INC}
endif

########### Modified for the ccefw-desktop environment ############
ifeq ($(HOST),ccefw-desktop)
 CC = g++
 CC = /opt/intel/impi/2017.2.174/bin64/mpicxx
 CC = CC -mkl -no-ipo
 CC = /opt/intel/impi/2017.2.174/bin64/mpiicpc
 PATH := ${PATH}:/opt/intel/bin
 INC = -I./ -I../ -I/opt/intel/mkl/include/
 MKLROOT = /opt/intel/mkl/
endif
################################################################################


########### added by Gian for PCL mini-cluster #################################
ifeq ($(HOST),pcl-tiergarten.sc.intel.com)
 PATH := ${PATH}:/swtools/intel/bin
 INC += -I/nfs_home/giangiac/Intel-QS -I/swtools/intel/mkl/include/
 MKLROOT = /swtools/intel/mkl/
endif
################################################################################

