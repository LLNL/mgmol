FROM ubuntu:22.04

ENV ENVDIR=env

# install sudo
RUN apt-get -yq update && apt-get -yq install sudo

WORKDIR /$ENVDIR

# install packages
RUN sudo apt-get install -yq git
RUN sudo apt-get install --no-install-recommends -yq make gcc gfortran libssl-dev cmake
RUN sudo apt-get install -yq libopenblas-dev libmpich-dev libblas-dev liblapack-dev libscalapack-mpi-dev libhdf5-mpich-dev
RUN sudo apt-get install -yq libboost-all-dev
RUN sudo apt-get install -yq vim
RUN sudo apt-get install -yq git-lfs
RUN sudo apt-get install -yq valgrind
RUN sudo apt-get install -yq wget
RUN sudo apt-get install -yq astyle

# install lldb and gdb for debugging
RUN sudo apt-get install -yq lldb gdb

RUN sudo apt-get clean -q

# NOTE: currently docker does not have a way to set environment variable with a command output.
# The following environment variable should use $(uname -m) for different architecture.
# ENV SCALAPACK_ROOT=/usr/lib/aarch64-linux-gnu/
ENV SCALAPACK_ROOT=/usr/lib/x86_64-linux-gnu/

# create and switch to a user
ENV USERNAME=test
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
RUN useradd --no-log-init -u 1001 --create-home --shell /bin/bash $USERNAME
RUN adduser $USERNAME sudo
USER $USERNAME
WORKDIR /home/$USERNAME

