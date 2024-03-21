cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/

#!/bin/bash
samples=(1 2 3 4)
celltypes=("Adipocytes" "Epithelial_cells" "Hepatocytes")

salmon="../software/salmon-1.5.2/bin/salmon"

################################################################################################
# salmon index
################################################################################################
# use this one: transcriptome.expanded.fa
$salmon index \
-t 05_alevin/transcriptome_splici_fl96.fa \
-i 05_alevin/index \
--gencode -p 32

################################################################################################
# simulation (without DGE)
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim="04_simulation/simulation/normal${sample}_${type}_simulated_reads"
    R1_sim=$input_dir_sim"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim=$input_dir_sim"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim="05_alevin/simulation/normal${sample}_${type}"

    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim -2 $R2_sim -o $output_dir_alevin_sim \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation WITH DGE
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_DGE/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_DGE/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation FC6
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_FC6/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_FC6/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation FC9
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_FC9/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_FC9/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation NULL
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_NULL/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_NULL/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation drop90
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_drop90/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_drop90/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation drop95
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_drop95/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_drop95/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done

################################################################################################
# simulation drop99
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_drop99/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_drop99/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done


################################################################################################
# simulation batch
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    input_dir_sim_DGE="04_simulation/simulation_batch/normal${sample}_${type}_simulated_reads"
    R1_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R1_001.fastq.gz"
    R2_sim_DGE=$input_dir_sim_DGE"/hg_100_S1_L001_R2_001.fastq.gz"
    output_dir_alevin_sim_DGE="05_alevin/simulation_batch/normal${sample}_${type}"
    
    $salmon alevin -l ISR \
    -i 05_alevin/index \
    -1 $R1_sim_DGE -2 $R2_sim_DGE -o $output_dir_alevin_sim_DGE \
    -p 32 --chromium --sketch
  done
done


