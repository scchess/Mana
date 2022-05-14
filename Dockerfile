FROM conda/miniconda3
RUN apt-get update && \
    apt-get install -y git && \
    apt-get install -y wget && \
    apt-get install -y libncurses5-dev && \
    apt-get install -y build-essential

RUN git clone https://github.com/lh3/minimap2.git && \
    cd minimap2 && \
    make && \
    cp minimap2 /usr/local/bin/

RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 &&
    bzip2 -d samtools-1.15.1.tar.bz2 && \
    tar -xf samtools-1.15.1.tar && \
    cd samtools-1.15.1 && \
    make && \
    cp samtools /usr/local/bin

RUN conda install -c bioconda pysamstats
RUN conda install -c pytest