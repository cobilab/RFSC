FROM ubuntu:20.04

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt update && apt install -y python3-pip wget unzip curl bc libuv1-dev libdw-dev git

RUN apt-get install -y bc

RUN apt-get install -y curl

RUN apt-get install -y liblwp-protocol-https-perl

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

ADD . /RFSC

WORKDIR /RFSC

RUN ./RFSC.sh --install

CMD tail -f >> /dev/null