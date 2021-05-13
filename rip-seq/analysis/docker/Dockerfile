FROM continuumio/miniconda3:4.9.2

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y -c bioconda \
                                r-base=4.0.3 \
                                pandas=1.2.4 \
                                matplotlib=3.4.1 \
                                xlrd=2.0.1 \
                                openpyxl=3.0.7 \
                                bioconductor-deseq2=1.30.1 \
                                bioconductor-tximport=1.18.0 \
                                bioconductor-gviz=1.34.1 \
                                pybedtools=0.8.2 \
                                deeptools=3.5.1

# Not available in conda
RUN python3 -m pip install htseq==0.12.3

WORKDIR /project
