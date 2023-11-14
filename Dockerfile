# lean base image to start from
FROM alpine:3.17

# add labels
LABEL description="A PEPC build-image that contains all required tools\
GCC/12.2.1 20220924, OpenMPI/4.1.4, JUBE, FORD, fprettify"
LABEL version="1.1"

# install all required packages to build PEPC
# also get and build OpenMPI
RUN apk update && \
    apk add \
    automake \
    autoconf \
    ctags \
    build-base \
    git \
    gcc \
    gfortran \
    openmpi-dev \
    openmpi \
    openssh \
    curl \
    python3 \
    graphviz \
    py3-pip ;\
    pip3 install http://apps.fz-juelich.de/jsc/jube/jube2/download.php?version=latest; \
    pip3 install FORD==6.1.10; \
    pip3 install --upgrade git+https://github.com/dbroemmel/fprettify.git@fixed_relations_and_case

# try and have start dir?
WORKDIR /tmp
