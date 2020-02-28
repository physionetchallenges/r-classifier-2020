FROM r-base:latest

## The MAINTAINER instruction sets the Author field of the generated images
MAINTAINER author@sample.com
## DO NOT EDIT THESE 3 lines
RUN mkdir /physionet
COPY ./ /physionet
WORKDIR /physionet

## Install your dependencies here using apt-get etc.
RUN Rscript requirements.txt
