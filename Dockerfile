FROM ubuntu:20.04

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt update && apt install -y python3-pip wget unzip curl bc

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

ADD . /RFSC

WORKDIR /RFSC

RUN conda install -c conda-forge ncbi-datasets-cli

RUN conda install -c bioconda perl-lwp-protocol-https

RUN conda install -c bioconda sra-tools entrez-direct

RUN conda install perl-io-socket-ssl perl-net-ssleay perl-lwp-protocol-https entrez-direct

RUN ./RFSC.sh --install

CMD tail -f >> /dev/null