FROM python:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    g++ \
    make \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install build setuptools wheel pybind11 numpy pytest

WORKDIR /workspace

CMD ["/bin/bash"]
