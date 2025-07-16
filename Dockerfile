FROM python:latest

# Preventing apt-get from asking interactive questions during install
ENV DEBIAN_FRONTEND=noninteractive

# Installing dependencies used to compile nthash library and the rest
RUN apt-get update && apt-get install -y \
    pkg-config \
    g++ \
    make \
    git \
    doxygen \
    meson \
    zlib1g-dev \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/* 
# Last command was to clean up to reduce image size

# BUILDING AND INSTALLING NTHASH ONCE
RUN git clone --recursive https://github.com/bcgsc/ntHash.git && \
    cd ntHash && \
    sed -i '2i#include <cstdint>' ./src/internal.hpp && \
    meson setup --buildtype=release --prefix=/usr/local build && \
    meson compile -C build && \
    meson install -C build && \
    cd .. && rm -rf ntHash

# Installing Python build tools
RUN pip3 install build setuptools wheel pybind11

ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
ENV PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:${PKG_CONFIG_PATH}"

WORKDIR /workspace

CMD ["/bin/bash"]



