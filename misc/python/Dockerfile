FROM python:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    pkg-config \
    g++ \
    make \
    git \
    meson \
    zlib1g-dev \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Build and install ntHash (needed by misc/tests benchmark and test suite)
RUN git clone --recursive https://github.com/bcgsc/ntHash.git && \
    cd ntHash && \
    sed -i '2i#include <cstdint>' ./src/internal.hpp && \
    meson setup --buildtype=release --prefix=/usr/local build && \
    meson compile -C build && \
    meson install -C build && \
    cd .. && rm -rf ntHash

RUN pip3 install build setuptools wheel pybind11 numpy pytest

ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib"
ENV PKG_CONFIG_PATH="/usr/local/lib/pkgconfig"

WORKDIR /workspace

CMD ["/bin/bash"]
