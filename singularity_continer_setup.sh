sudo singularity build --sandbox riri docker://ubuntu
sudo singularity shell --writable riri

apt-get update && apt-get install -y build-essential
apt install -y zlib1g-dev libssl-dev
apt-get install -y libcurl4-openssl-dev libncurses-dev libbz2-dev liblzma-dev r-base-core libxml2-dev libcurl4-openssl-dev curl bc openjdk-8-jdk zip wget


wget -O fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc.zip
cd ./FastQC
chmod 755 fastqc
ln -s ${PWD}/fastqc /usr/bin/fastqc
cd


wget -O samtools.tar.bz2 https://iweb.dl.sourceforge.net/project/samtools/samtools/1.14/samtools-1.14.tar.bz2
tar xjf samtools.tar.bz2
cd samtools-1.*
./configure
make
make install
cd 

wget -O bcftools.tar.bz2 https://iweb.dl.sourceforge.net/project/samtools/samtools/1.14/bcftools-1.14.tar.bz2
tar xjf bcftools.tar.bz2
cd bcftools-1.*
./configure
make
make install
cd 

wget -O htslib.tar.bz2 https://iweb.dl.sourceforge.net/project/samtools/samtools/1.14/htslib-1.14.tar.bz2
tar xjf htslib.tar.bz2
cd htslib-1.*
./configure
make
make install
cd 

wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
mv Trimmomatic-0.39/ /usr/bin/Trimmomatic-0.39
rm Trimmomatic-0.39.zip

wget -O subread.tar.gz https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz/download
tar -xzf subread.tar.gz
rm subread.tar.gz
mv subread-2.0.3-Linux-x86_64 /usr/bin/subread

apt-get install -y python2 python3
apt-get -y install python3-pip

pip install multiqc
pip install numpy pysam cython
pip install plastid

wget -O bedtools https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
chmod +x bedtools
mv bedtools /usr/bin/bedtools

wget -O bowtie.zip https://versaweb.dl.sourceforge.net/project/bowtie-bio/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip
unzip bowtie.zip
rm bowtie.zip
mv ./bowtie-1.3.1-linux-x86_64/* /usr/bin/


wget -O bowtie2.zip https://iweb.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip
unzip bowtie2.zip
rm bowtie2.zip
mv ./bowtie2-2.4.5-linux-x86_64/* /usr/bin/

R -e 'install.packages("BiocManager")'
R -e 'BiocManager::install("NOISeq")'
R -e 'BiocManager::install("Repitools")'
wget -O qualimap.tar.gz https://bitbucket.org/kokonech/qualimap/downloads/qualimap-build-31-08-20.tar.gz
tar -xzf qualimap.tar.gz
rm qualimap.tar.gz
Rscript ./qualimap-build-31-08-20/scripts/installDependencies.r
mv ./qualimap-build-31-08-20/* /usr/bin/

wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar
mv picard.jar /usr/bin/picard.jar
##############

#######
sudo singularity build riri.sif riri
singularity push -U riri.sif library://semiquant/default/riri:v0.1









