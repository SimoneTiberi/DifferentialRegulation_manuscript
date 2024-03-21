cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/

alevin_fry="../software/alevin-fry/target/release/alevin-fry"

#!/bin/bash
samples=(1 2 3 4)
celltypes=("Adipocytes" "Epithelial_cells" "Hepatocytes")

################################################################################################
# simulation (withou DGE)
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim="05_alevin/simulation/normal${sample}_${type}"
    output_dir_alevin_fry_sim="06_alevin_fry_USA/simulation/normal${sample}_${type}"

    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim \
    -o $output_dir_alevin_fry_sim \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim \
    -i $output_dir_alevin_fry_sim \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim \
    -o $output_dir_alevin_fry_sim \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done

################################################################################################
# simulation with DGE
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_DGE/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_DGE/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done

################################################################################################
# simulation FC6
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_FC6/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_FC6/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done


################################################################################################
# simulation FC9
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_FC9/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_FC9/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done


################################################################################################
# simulation NULL
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_NULL/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_NULL/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done

################################################################################################
# simulation drop90
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_drop90/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_drop90/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done

################################################################################################
# simulation drop95
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_drop95/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_drop95/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done

################################################################################################
# simulation drop99
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_drop99/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_drop99/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done

################################################################################################
# simulation batch
################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
    tx2="05_alevin/transcriptome_splici_fl96_t2g_3col.tsv"
    output_dir_alevin_sim_DGE="05_alevin/simulation_batch/normal${sample}_${type}"
    output_dir_alevin_fry_sim_DGE="06_alevin_fry_USA/simulation_batch/normal${sample}_${type}"
    
    $alevin_fry generate-permit-list \
    -i $output_dir_alevin_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -k -d fw
    
    $alevin_fry collate \
    -r $output_dir_alevin_sim_DGE \
    -i $output_dir_alevin_fry_sim_DGE \
    -t 32
    
    $alevin_fry quant \
    -i $output_dir_alevin_fry_sim_DGE \
    -o $output_dir_alevin_fry_sim_DGE \
    -m $tx2 \
    -t 32 --use-mtx -d -r cr-like-em
  done
done
