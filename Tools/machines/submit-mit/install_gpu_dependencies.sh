#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail

module load mpi


rm -rf $HOME/.conda/envs/warpx
conda config --set auto_activate_base false
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --yes -n warpx -c conda-forge mamba conda-libmamba-solver
conda activate warpx
conda config --set solver libmamba
mamba install --yes -c conda-forge adios2 \
                                   blaspp \
                                   boost \
                                   ccache \
                                   cmake \
                                   cython \
                                   dask \
                                   dask-jobqueue \
                                   fast-histogram \
                                   fftw \
                                   h5py \
                                   ipykernel \
                                   ipympl \
                                   lapackpp \
                                   libadios2 \
                                   libgomp \
                                   make \
                                   matplotlib \
                                   mpi4py \
                                   numpy \
                                   openpmd-api \
                                   openpmd-viewer \
                                   packaging \
                                   pandas \
                                   pip \
                                   pkg-config \
                                   pyarrow \
                                   python \
                                   python-build \
                                   scipy \
                                   setuptools \
                                   virtualenv \
                                   wheel \
                                   yt \
                                   openmpi=4.1.*=external_*
