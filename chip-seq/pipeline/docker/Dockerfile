FROM continuumio/miniconda3:4.9.2

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y -c bioconda \
                                python=3.7 \
                                fastqc=0.11.9 \
                                samtools=1.12 \
                                macs2=2.2.6 \
                                minimap2=2.18 \
                                deeptools=3.5.1 \
                                snakemake-minimal=5.26.0 \
                                multiqc=1.10.1 \
                                bedtools=2.30.0 \
                                ucsc-bedclip=377 \
                                ucsc-bedgraphtobigwig=377 \
                                trim-galore=0.6.6

# samtools fix
RUN cd /opt/conda/lib && ln -s libcrypto.so.1.1 libcrypto.so.1.0.0

# fresh sratoolkit
WORKDIR /tools
RUN wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz && \
    tar -vxzf sratoolkit.tar.gz
ENV PATH "$PATH:/tools/sratoolkit.2.10.8-ubuntu64/bin"

WORKDIR /project
