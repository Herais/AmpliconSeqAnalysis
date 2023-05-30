#!/bin/bash

# Define Path for colab
DPATH_TOOLS = "$1"
echo $DPATH_TOOLS

# Biopython
pip3 install -q jupyter biopython
pip install -q nglview

# logomaker
pip install logomaker

# pysam
pip install pysam

# squarify
pip install squarify

# install weblogo
pip install weblogo


# Install Bowtie2
DPATH_bowtie2="$DPATH_TOOLS/bowtie2-2.4.2-sra-linux-x86_64"
if [ ! -d "$DPATH_bowtie2" ]; then
  cd /content
  url_download="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download"
  wget -nc "$url_download"
  unzip -q -o download -d "$DPATH_TOOLS"
fi

# Install NGmerge
DPATH_NGmerge="$DPATH_TOOLS/NGmerge"
if [ ! -d "$DPATH_TOOLS/NGmerge" ]; then
  git clone https://github.com/harvardinformatics/NGmerge.git -C "$DPATH_TOOLS"
fi
cd "$DPATH_NGmerge"
make

# Install Fastq-join
DPATH_fastqjoin="$DPATH_TOOLS/fastq-join"

if [ ! -d "$DPATH_fastqjoin" ]; then
  git clone https://github.com/brwnj/fastq-join "$DPATH_fastqjoin"
fi
cd "$DPATH_fastqjoin"
make

# Install samtools
DPATH_samtools="$DPATH_TOOLS/samtools-1.17"

if [ ! -d "$DPATH_samtools" ]; then
  cd /content
  url_download='https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2'
  wget -wc "$url_download"
  tar -xjf "/content/samtools-1.17.tar.bz2" -C "$DPATH_TOOLS"
fi

cd "$DPATH_samtools"
./configure
make -s
make -s install

# Install clustal
DPATH_clustalw="$DPATH_TOOLS/clustalw-2.1"

if [ ! -d "$DPATH_clustalw" ]; then
  cd /content
  url_download='http://www.clustal.org/download/current/clustalw-2.1.tar.gz'
  wget "$url_download"
  tar -xzf "/content/clustalw-2.1.tar.gz" -C "$DPATH_TOOLS"
fi

cd "$DPATH_clustalw"
./configure
make
make install