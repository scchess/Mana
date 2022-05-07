FROM ubuntu
RUN git clone git@github.com:lh3/minimap2.git && \
    cd minimap2 && \
    make && \
    cp minimap2 /usr/local/bin/
