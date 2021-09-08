FROM ubuntu:20.04

# tzdata package is required for r-base
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow

RUN apt-get update && apt-get install -y \
    python3.8 python3.8-dev python3-pip \
    bedtools=2.27.1+dfsg-4ubuntu1 \
    r-base=3.6.3-2 \
    libssl-dev=1.1.1f-1ubuntu2.4 \
    curl=7.68.0-1ubuntu2.5 \
    libcurl4-openssl-dev=7.68.0-1ubuntu2.5\
    libxml2-dev=2.9.10+dfsg-5ubuntu0.20.04.1 \
    clustalo=1.2.4-4build1 \
    build-essential libtool-bin automake \
    && rm -rf /var/lib/apt/lists/*


RUN python3.8 -m pip install \
    pybedtools==0.8.2 \
    tqdm==4.60.0 \
    numpy==1.20.1 \
    matplotlib==3.1.3 \
    HTSeq==0.12.3 \
    rpy2==3.3.6 \
    pandas==1.1.5 \
    matplotlib_venn==0.11.6 \
    biopython==1.76 \
    parasail==1.2.2

# requires numpy installed
RUN python3.8 -m pip install \
    pyBigWig==0.3.17


WORKDIR /project
