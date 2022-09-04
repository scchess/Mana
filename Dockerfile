FROM r-base
RUN apt-get update && \
    apt-get install -y git && \
    apt-get install -y wget && \
    apt-get install -y cmake && \
    apt-get install -y libssl-dev && \
    apt-get install -y libfontconfig1-dev && \
    apt-get install -y libbz2-dev && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libncurses5-dev && \
    apt-get install -y build-essential

RUN git clone https://github.com/lh3/minimap2.git && \
    cd minimap2 && \
    make && \
    cp minimap2 /usr/local/bin/

RUN R -e "install.packages('svglite',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggpubr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    bzip2 -d samtools-1.15.1.tar.bz2 && \
    tar -xf samtools-1.15.1.tar && \
    cd samtools-1.15.1 && \
    ./configure && \
    make && \
    cp samtools /usr/local/bin

RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
    bzip2 -d bcftools-1.16.tar.bz2 && \
    tar -xf bcftools-1.16.tar && \
    cd bcftools-1.16 && \
    ./configure && \
    make && \
    cp bcftools /usr/local/bin/

RUN apt-get install -y python2
RUN apt-get install -y python3.9
RUN apt-get install -y python3-pip
RUN apt-get install python-is-python3
RUN pip3 install nanosim-h
RUN pip3 install cython
RUN apt-get install -y libcurl4-openssl-dev
RUN pip3 install "setuptools<58.0.0"
RUN apt-get install -y python3-dev

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
    tar -xzf bedtools-2.30.0.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/bedtools /usr/local/bin

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh && \
    bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -b

RUN /root/miniconda3/bin/conda config --add channels defaults
RUN /root/miniconda3/bin/conda config --add channels conda-forge
RUN /root/miniconda3/bin/conda config --add channels bioconda
RUN /root/miniconda3/bin/conda install -y pysam
RUN /root/miniconda3/bin/conda install -y pysamstats

ENV PATH "$PATH:/root/miniconda3/bin"
RUN pip3 install pytest
RUN pip3 install pandas
RUN pip3 install tabulate

WORKDIR /src
COPY . /src
