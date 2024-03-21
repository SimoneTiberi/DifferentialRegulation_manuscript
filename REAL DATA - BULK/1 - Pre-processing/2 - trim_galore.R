cd /home/Shared_sherborne/simone/Diff_Velo/software

# Install Trim Galore
# curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
# tar xvzf trim_galore.tar.gz
# Run Trim Galore
# TrimGalore-0.6.10/trim_galore

trim_galore="/home/Shared_sherborne/simone/Diff_Velo/software/TrimGalore-0.6.10/trim_galore"

# Directory where fastqc files are stored:
# base-directory:
cd /home/Shared_sherborne/simone/Diff_Velo/discovery_bulk
# mkdir 02_trim

$trim_galore -q 20 --phred33 --length 20 \
--illumina -o 02_trim \
--fastqc --gzip \
--paired 01_data/ERR4579267_1.fastq 01_data/ERR4579267_2.fastq

$trim_galore -q 20 --phred33 --length 20 \
--illumina -o 02_trim \
--fastqc --gzip \
--paired 01_data/ERR4579268_1.fastq 01_data/ERR4579268_2.fastq

$trim_galore -q 20 --phred33 --length 20 \
--illumina -o 02_trim \
--fastqc --gzip \
--paired 01_data/ERR4579269_1.fastq 01_data/ERR4579269_2.fastq

$trim_galore -q 20 --phred33 --length 20 \
--illumina -o 02_trim \
--fastqc --gzip \
--paired 01_data/ERR4579270_1.fastq 01_data/ERR4579270_2.fastq

$trim_galore -q 20 --phred33 --length 20 \
--illumina -o 02_trim \
--fastqc --gzip \
--paired 01_data/ERR4579271_1.fastq 01_data/ERR4579271_2.fastq

$trim_galore -q 20 --phred33 --length 20 \
--illumina -o 02_trim \
--fastqc --gzip \
--paired 01_data/ERR4579272_1.fastq 01_data/ERR4579272_2.fastq

####################################################################################################
# in bash:
####################################################################################################
cd trim
multiqc . --interactive
