# get 
cd /home/Shared_sherborne/simone/Diff_Velo/discovery_bulk
grep ">" reference/GRCm39.primary_assembly.genome.fa | cut -d ">" -f 2 | cut -d " " -f 1 > reference/GRCm39.primary_assembly.chromosome_names.txt

## modify path to salmon
salmon="/home/Shared_sherborne/simone/Diff_Velo/software/salmon-1.5.2/bin/salmon"
# typing $salmon is the same as typing /home/stiber/software/salmon-0.12.0_linux_x86_64/bin/salmon

# the reference transcriptome
input_txome="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/reference/gencode.vM33.annotation.gz.expanded.transcripts.fa"
input_genome="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/reference/GRCm39.primary_assembly.genome.fa"
input_chrnames="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/reference/GRCm39.primary_assembly.chromosome_names.txt"

# Directory where I want to build index (Salmon_index is the name of the index to be created)
idx="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/Salmon_index"

# Directory where fastqc files are stored (TRIMMED READS, not original ones!):
fastqDir="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/02_trim"

# Directory where I want to write the Salmon output (alignment and quantification):
salmonDir="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon"

# bash command to build the index (~1 hour):
$salmon index -t <(cat $input_txome $input_genome) -k 31 -i $idx --gencode -p 10 -d $input_chrnames

# each sample will go in a separate directory, created by salmon, but specified by us (e.g. -o $salmonDir/sample_1)
# Run the alignment & quantification step with Salmon (~30-60 mins per sample).
$salmon quant -i $idx -l A \
-1 $fastqDir/ERR4579267_1_val_1.fq.gz -2 $fastqDir/ERR4579267_2_val_2.fq.gz \
-p 10 -o $salmonDir/ERR4579267 --dumpEq --validateMappings \
--seqBias --gcBias

$salmon quant -i $idx -l A \
-1 $fastqDir/ERR4579268_1_val_1.fq.gz -2 $fastqDir/ERR4579268_2_val_2.fq.gz \
-p 10 -o $salmonDir/ERR4579268 --dumpEq --validateMappings \
--seqBias --gcBias

$salmon quant -i $idx -l A \
-1 $fastqDir/ERR4579269_1_val_1.fq.gz -2 $fastqDir/ERR4579269_2_val_2.fq.gz \
-p 10 -o $salmonDir/ERR4579269 --dumpEq --validateMappings \
--seqBias --gcBias

$salmon quant -i $idx -l A \
-1 $fastqDir/ERR4579270_1_val_1.fq.gz -2 $fastqDir/ERR4579270_2_val_2.fq.gz \
-p 10 -o $salmonDir/ERR4579270 --dumpEq --validateMappings \
--seqBias --gcBias

$salmon quant -i $idx -l A \
-1 $fastqDir/ERR4579271_1_val_1.fq.gz -2 $fastqDir/ERR4579271_2_val_2.fq.gz \
-p 10 -o $salmonDir/ERR4579271 --dumpEq --validateMappings \
--seqBias --gcBias

$salmon quant -i $idx -l A \
-1 $fastqDir/ERR4579272_1_val_1.fq.gz -2 $fastqDir/ERR4579272_2_val_2.fq.gz \
-p 10 -o $salmonDir/ERR4579272 --dumpEq --validateMappings \
--seqBias --gcBias
