FROM ubuntu

MAINTAINER Marc Chevrette <chevrette@wisc>
RUN apt-get update
RUN apt-get install -y wget git default-jre mafft ncbi-blast+ hmmer bioperl python-pip
RUN pip install --upgrade pip
RUN pip install scipy sklearn numpy
RUN wget http://www.microbesonline.org/fasttree/FastTree
RUN chmod 777 FastTree
RUN git clone https://chevrm@bitbucket.org/chevrm/sandpuma.git
RUN mv FastTree /sandpuma/dependencies
RUN mv /sandpuma/docker/NRPSpredictor2.sh /sandpuma/dependencies/NRPSPredictor2/NRPSpredictor2.sh

ENV PATH=$PATH:/sandpuma/dependencies