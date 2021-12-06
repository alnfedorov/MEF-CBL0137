FROM continuumio/miniconda3:4.9.2

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y -c bioconda \
                                r-base=4.1.1 \
                                pandas=1.2.4 \
                                matplotlib=3.4.1 \
                                bioconductor-deseq2=1.32.0 \
                                bioconductor-apeglm=1.14.0 \
                                bioconductor-tximport=1.20.0 \
                                bioconductor-degreport=1.28.0 \
                                pybedtools=0.8.2 \
                                deeptools=3.5.1 \
                                seaborn=0.11.2

# Not available in conda
RUN python3 -m pip install htseq==0.12.3

WORKDIR /project
