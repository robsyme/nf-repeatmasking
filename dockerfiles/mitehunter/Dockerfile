FROM ubuntu:17.10

MAINTAINER Rob Syme <rob.syme@gmail.com>

RUN apt-get update && apt-get install -qqy \
    wget \
    unzip \
	build-essential \
	ncbi-blast+-legacy

WORKDIR /usr/local

# Install mDUST (MITE-Hunter prereq)
RUN wget ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/seqclean/mdust.tar.gz \
	&& tar -xvf mdust.tar.gz \
	&& rm mdust.tar.gz \
	&& cd mdust \
	&& make

# Install Muscle and blast2 (MITE-Hunter prereqs)
RUN apt-get install -qqy muscle 

# RUN wget https://github.com/robsyme/MITE-Hunter/archive/v1.0.1.zip \
# 	&& unzip v1.0.1.zip \
# 	&& mv MITE* mitehunter \
# 	&& cd mitehunter \
# 	&& perl MITE_Hunter_Installer.pl -d /usr/local/mitehunter/ -f /usr/bin/formatdb -b /usr/bin/blastall -m /usr/local/mdust/mdust -M /usr/bin/muscle

ADD code /usr/local/mitehunter

RUN cd /usr/local/mitehunter \
	&& perl MITE_Hunter_Installer.pl -d /usr/local/mitehunter/ -f /usr/bin/formatdb -b /usr/bin/blastall -m /usr/local/mdust/mdust -M /usr/bin/muscle

ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/mitehunter




