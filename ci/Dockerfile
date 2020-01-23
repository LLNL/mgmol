FROM nvidia/cuda:10.1-devel

ENV N_PROCS=8

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
      gcc \
      gfortran \
      build-essential \
      clang-format \
      clang \
      wget \
      curl \
      lcov \
      bison \
      python3 \
      environment-modules \
      tcl-dev \
      libperlio-gzip-perl \
      libjson-pp-perl \
      git && \
      apt-get clean && rm -rf /var/lib/apt/lists/*

ENV PREFIX=/root

RUN export SPACK_HASH=ae52f4144872a281bc97596ecce17d9959c21086 && \
    export SPACK_URL=https://github.com/spack/spack/archive/${SPACK_HASH}.tar.gz && \
    export SPACK_ARCHIVE=${PREFIX}/spack.tar.gz && \
    export SPACK_SOURCE_DIR=${PREFIX}/spack && \
    wget --quiet ${SPACK_URL} --output-document=${SPACK_ARCHIVE} && \
    mkdir -p ${SPACK_SOURCE_DIR} && \
    tar xf ${SPACK_ARCHIVE} -C ${SPACK_SOURCE_DIR} --strip-components=1 && \
    rm -rf ${SPACK_ARCHIVE}
ENV PATH=${PREFIX}/spack/bin:$PATH

# Register the compiler
RUN spack compilers

# Make sure we don't install CUDA
RUN echo 'packages:' >> ${PREFIX}/.spack/packages.yaml && \
    echo '  cuda:' >>  ${PREFIX}/.spack/packages.yaml && \
    echo '    paths:' >>  ${PREFIX}/.spack/packages.yaml && \
    echo '      cuda@10.1%gcc@7.4.0 arch=linux-ubuntu18.04-x86_64: /usr/local/cuda' >>  ${PREFIX}/.spack/packages.yaml && \
    echo '    buildable: False' >>  ${PREFIX}/.spack/packages.yaml

# Install MAGMA with nvcc and gcc 7.4 
RUN spack install magma@2.5.1^openblas@0.3.7+virtual_machine arch=x86_64

# Install gcc 9.2 with offloading support
# Checksum is missing in spack release
RUN spack install gcc@9.2.0+nvptx arch=x86_64 
RUN spack compiler find $(spack location -i gcc)

# Install all the libraries with gcc 9.2
# Install boost
RUN spack install boost@1.70.0 arch=x86_64
# Install openmpi
RUN spack install openmpi@4.0.2+thread_multiple arch=x86_64
# Install hdf5 with mpi for some reason zlib does not respect the arch in this
# case so we set it explicitly
RUN spack install hdf5@1.10.5+hl^openmpi@4.0.2+thread_multiple^zlib@1.2.11 arch=x86_64
# Install scalapack using the version of openblas that we have installed
RUN spack install netlib-scalapack@2.0.2%gcc@9.2.0^openmpi@4.0.2%gcc@9.2.0+thread_multiple^openblas@0.3.7%gcc@7.4.0+virtual_machine arch=x86_64

# Install the development version of lcov because older versions are not
# compatible with gcc 9
RUN cd $PREFIX && git clone https://github.com/linux-test-project/lcov.git && \
    cd lcov && make install && cd ../ && rm -rf lcov

# Set the environment variables to be able to run mpi as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN echo 'module() { eval `/usr/bin/modulecmd sh $*`; }' > ~/.module
RUN echo '. /root/spack/share/spack/setup-env.sh' >> ~/.module
