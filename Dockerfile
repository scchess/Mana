FROM conda/miniconda3
RUN apt-get update && \
    apt-get install -y git && \
    apt-get install -y build-essential
RUN git clone https://github.com/lh3/minimap2.git && \
    cd minimap2 && \
    make && \
    cp minimap2 /usr/local/bin/
RUN conda install -c bioconda pysamstats
RUN conda install -c pytest
