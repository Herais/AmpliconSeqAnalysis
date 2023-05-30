# Define Path for colab
dpath_tools="$1"
echo $dpath_tools

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
dpath_bowtie2="$dpath_tools/bowtie2-2.4.2-sra-linux-x86_64"
if [ ! -d "$dpath_bowtie2" ]; then
  cd /content
  url_download="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download"
  wget -nc "$url_download"
  unzip -q -o download -d "$dpath_tools"
fi

# Install NGmerge
dpath_ngmerge="$dpath_tools/NGmerge"
if [ ! -d "$dpath_tools/NGmerge" ]; then
  git clone https://github.com/harvardinformatics/NGmerge.git -C "$dpath_tools"
fi
cd "$dpath_ngmerge"
make

# Install Fastq-join
dpath_fastqjoin="$dpath_tools/fastq-join"

if [ ! -d "$dpath_fastqjoin" ]; then
  git clone https://github.com/brwnj/fastq-join "$dpath_fastqjoin"
fi
cd "$dpath_fastqjoin"
make

# Install samtools
dpath_samtools="$dpath_tools/samtools-1.17"

if [ ! -d "$dpath_samtools" ]; then
  cd /content
  url_download='https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2'
  wget -wc "$url_download"
  tar -xjf "/content/samtools-1.17.tar.bz2" -C "$dpath_tools"
fi

cd "$dpath_samtools"
./configure
make -s
make -s install

# Install clustal
dpath_clustalw="$dpath_tools/clustalw-2.1"

if [ ! -d "$dpath_clustalw" ]; then
  cd /content
  url_download='http://www.clustal.org/download/current/clustalw-2.1.tar.gz'
  wget "$url_download"
  tar -xzf "/content/clustalw-2.1.tar.gz" -C "$dpath_tools"
fi

cd "$dpath_clustalw"
./configure
make
make install