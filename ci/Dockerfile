FROM nvidia/cuda:10.1-devel

ENV N_PROCS=8

# Note that libboost-all-dev automatically pull libopenmpi-dev. Thus we cannot
# install our own MPI library.
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
      gcc \
      gfortran \
      build-essential \
      wget \
      curl \
      bison \
      git \
      pkg-config \
      libopenblas-dev \
      libboost-all-dev \
      libscalapack-openmpi-dev \
      libhdf5-openmpi-dev && \
      apt-get clean && rm -rf /var/lib/apt/lists/*

ENV PREFIX=/scratch
ENV ARCHIVE_DIR=${PREFIX}/archive
ENV SOURCE_DIR=${PREFIX}/source
ENV BUILD_DIR=${PREFIX}/build
ENV INSTALL_DIR=/opt

RUN mkdir -p ${PREFIX} && \
    cd ${PREFIX} && \
    mkdir archive && \
    mkdir source && \
    mkdir build

# Install CMake
RUN export CMAKE_VERSION=3.15.0 && \
    export CMAKE_VERSION_SHORT=3.15 && \
    export CMAKE_SHA256=5a681e37b142006ab44aef9b976db610dc5d94124da41743ff53412f4b36a7fe && \
    export CMAKE_URL=https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz && \
    export CMAKE_ARCHIVE=${ARCHIVE_DIR}/cmake.tar.gz && \
    export CMAKE_BUILD_DIR=${BUILD_DIR}/cmake && \
    wget --quiet ${CMAKE_URL} --output-document=${CMAKE_ARCHIVE} && \
    echo "${CMAKE_SHA256} ${CMAKE_ARCHIVE}" | sha256sum -c && \
    mkdir -p ${CMAKE_BUILD_DIR} && \
    tar xf ${CMAKE_ARCHIVE} -C ${CMAKE_BUILD_DIR} --strip-components=1 && \
    mv ${CMAKE_BUILD_DIR} ${INSTALL_DIR} && \
    rm -rf ${CMAKE_ARCHIVE} && \
    rm -rf ${CMAKE_BUILD_DIR}
ENV PATH=${INSTALL_DIR}/cmake/bin:$PATH

# Install MAGMA
RUN export MAGMA_VERSION=2.5.0 && \
    export MAGMA_URL=http://icl.utk.edu/projectsfiles/magma/downloads/magma-${MAGMA_VERSION}.tar.gz && \
    export MAGMA_ARCHIVE=${ARCHIVE_DIR}/magma.tar.gz && \
    export MAGMA_SOURCE_DIR=${SOURCE_DIR}/magma && \
    export MAGMA_INSTALL_DIR=${INSTALL_DIR}/magma && \
    export OPENBLASDIR=/usr/lib/x86_64-linux-gnu && \
    export CUDADIR=/usr/local/cuda && \
    wget --quiet ${MAGMA_URL} --output-document=${MAGMA_ARCHIVE} && \
    mkdir -p ${MAGMA_SOURCE_DIR} && \
    tar -xf ${MAGMA_ARCHIVE} -C ${MAGMA_SOURCE_DIR} --strip-components=1 && \
    cd ${MAGMA_SOURCE_DIR} && \
    cat make.check-openblas make.check-openblas > make.inc && \
    echo "prefix=/opt/magma" >> make.inc && \
    echo "FORT=gfortran" >> make.inc && \
    echo "GPU_TARGET=sm_70" >> make.inc && \
    sed -i "s/#FPIC/FPIC/" Makefile && \
    make -j${N_PROCS} install && \
    rm -rf ${MAGMA_ARCHIVE} && \
    rm -rf ${MAGMA_SOURCE_DIR}
ENV MAGMA_DIR=/opt/magma

# Append the option flag --allow-run-as-root to mpiexec
RUN echo '#!/usr/bin/env bash' > /usr/local/bin/mpiexec && \
    echo '/usr/bin/mpiexec --allow-run-as-root "$@"' >> /usr/local/bin/mpiexec && \
    chmod +x /usr/local/bin/mpiexec && \
    echo '#!/usr/bin/env bash' > /usr/local/bin/mpirun && \
    echo '/usr/bin/mpirun --allow-run-as-root "$@"' >> /usr/local/bin/mpirun && \
    chmod +x /usr/local/bin/mpirun
