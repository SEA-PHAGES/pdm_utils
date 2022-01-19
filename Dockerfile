# pdm_utils
# A commandline utility for building and maintaining Phamerator databases

FROM ubuntu:18.04

MAINTAINER christian.gauthier@pitt.edu
RUN mkdir /home/software
WORKDIR /home/software

# Perform all apt updates and installs, then cleanup to reduce image size
RUN apt update && apt upgrade -y
RUN apt install -y --no-install-recommends software-properties-common && add-apt-repository -y ppa:deadsnakes && apt install -y --no-install-recommends python3.8-venv
RUN apt install mysql-server-5.7 gcc wget infernal make autoconf m4 ncbi-blast+ clustalo -y
RUN apt clean -y && rm -rf /var/lib/apt/lists/*

# Set Python 3.8 as the default (Ubuntu 18.04 default is Python 3.6.9)
RUN python3.8 -m venv /venv
ENV PATH=/venv/bin:$PATH

# Set MySQL password
RUN service mysql start && mysql -u root -e "ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'phage'; FLUSH PRIVILEGES;"

# Install Aragorn
RUN wget http://www.ansikte.se/ARAGORN/Downloads/aragorn1.2.41.c && gcc -O3 -ffast-math -finline-functions -o /usr/local/bin/aragorn aragorn1.2.41.c && rm aragorn1.2.41.c

# Install tRNAscan-SE
RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.9.tar.gz && tar -xzf trnascan-se-2.0.9.tar.gz && rm trnascan-se-2.0.9.tar.gz && cd tRNAscan-SE-2.0 && ./configure && make && make install

# Install MMseqs2
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz && tar -xzf mmseqs-linux-sse41.tar.gz && rm mmseqs-linux-sse41.tar.gz && export PATH=$(pwd)/mmseqs/bin/:$PATH

# Install Markov clustering (MCL)
RUN wget https://micans.org/mcl/src/mcl-14-137.tar.gz && tar -xzf mcl-14-137.tar.gz && rm mcl-14-137.tar.gz && cd mcl-14-137 && ./configure && make && make install

# Get NCBI conserved domain database
RUN mkdir /home/Cdd
WORKDIR /home/Cdd
RUN wget https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz && tar -xzf Cdd_LE.tar.gz && rm Cdd_LE.tar.gz

# Finally, get pdm_utils and perform final setup steps
RUN mkdir /home/pdm_utils
WORKDIR /home/pdm_utils
RUN pip install pdm_utils

# Start mysql service each time we start the container
CMD service mysql start && tail -F /var/log/mysql/error.log