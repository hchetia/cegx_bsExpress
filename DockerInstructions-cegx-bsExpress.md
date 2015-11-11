# bsExpress-0.5cegx Docker Container Set Up #

1. Start Docker Quick Start Terminal 
    * Instructions ([Docker Tutorial](http://docs.docker.com/engine/installation/mac/))
    * Download and install Docker Toolbox [Docker Toolbox](https://www.docker.com/docker-toolbox)

2. Use a ubuntu based container

        docker pull ubuntu:14.04
        docker run -i -t a5a467fddcb8 /bin/bash
    
3. Set up install directory    

        mkdir Prerequisits
        apt-get update
        apt-get install wget
        apt-get install make
        apt-get install build-essential
        apt-get install zlib1g zlib1g-dev
        apt-get install ncurses-dev
        apt-get install unzip
        apt-get install git
        apt-get install oracle-java8-installer
        
## Prerequisits ##

### samtools ###
    wget http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2
    tar -zxvf samtools-1.1.tar.bz2
    cd samtools-1.1
    make
    make install
    
### FastQC ###
Already part of bsExpress-0.5cegx

### Trim_galore ###
Already part of bsExpress-0.5cegx    

### Bedtools ###
    apt-get install bedtools
    
### Cutadapt ###
    apt-get install python-pip
    apt-get install python-dev
    pip install --upgrade cutadapt 
    
### bismark ###
Already part of bsExpress-0.5cegx 

### bowtie2 ###
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip
    unzip bowtie2-2.2.6-linux-x86_64.zip
    cd bowtie2-2.2.6
    cp b* /usr/local/bin/

### BamUtils ###
    wget http://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.1.0.13.tgz
    tar zxvf BamUtilLibStatGen.1.0.13.tgz
    cd bamUtil_1.0.13/
    make all
    make install
    
### R ###
    apt-get -y install r-base
    
## bsExpress_0.5cegx ##

    git clone https://russellshamilton@bitbucket.org/russellshamilton/cegx_bsexpress.git
    cd cegx_bsexpress/
    
    python setup.py install --install-scripts /usr/local/bin
    
# bsExpress-0.5cegx Docker Save/Export #

Starting Attaching:

    docker start dkr_bsexpress-0.5cegx
    docker attach dkr_bsexpress-0.5cegx

Find the ID of the last run docker instance

    docker ps -l
CONTAINER ID | IMAGE | COMMAND | CREATED | STATUS | PORTS | NAMES
:-- | :-- | :-- | :-- | :-- | :-- | :-- |
6b724ce4677a | dkr_bsexpress-0.5cegx | "/bin/bash" | 4 minutes ago | Exited (0) 8 seconds ago | | sharp_raman

Commit the changes made

    docker commit -m "new auto_bsExpress script" 6b724ce4677a cegx_bsexpress_0.5

Exporting the image to a tar file

    docker save cegx_bsexpress_0.5  > cegx_bsexpress_0.5.tar
    #docker export cf039620c3c1 > cegx_bsExpress_dkr_container.tar