cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

awk -F "\t" '{if ($1 ~ /-U/) {print $0"-U"} else {print $0}}' 01_annotation_for_minnow/annotation.expanded.tx2gene.tsv > 01_annotation_for_minnow/minnow.tx2gene.tsv

sed -i -e 's/\r//g' 01_annotation_for_minnow/annotation.expanded.fa
sed -i -e 's/\r//g' 01_annotation_for_minnow/minnow.tx2gene.tsv

##################################################################################################
# minnow index
##################################################################################################
minnow="../software/minnow/build/src/minnow"

$minnow index -r 01_annotation_for_minnow/annotation.expanded.fa -k 101 -f 20 \
--tmpdir 04_simulation/tmp -p 20 \
-o 04_simulation/minnow_ind |& stdbuf -oL tr '\r' '\n' > \
04_simulation/minnow_index.log

##################################################################################################
# simulation (without DGE)
##################################################################################################
samples=(1 2 3 4)
celltypes=("Adipocytes" "Epithelial_cells" "Hepatocytes")
minnow_ind="04_simulation/minnow_ind"
g2t="01_annotation_for_minnow/minnow.tx2gene.tsv"
minnow_data="../software/minnow/data"

for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
      input_dir="04_simulation_minnow/simulation/normal$sample/$type"
      output_dir="04_simulation/simulation/normal${sample}_${type}_simulated_reads"
      log_file="04_simulation/simulation/minnow_simulate_normal${sample}_${type}.log"
      cmd="$minnow simulate --splatter-mode --g2t $g2t --inputdir $input_dir --PCR 6 -r $minnow_ind/ref_k101_fixed.fa -e 0.01 -p 20 -o $output_dir --dbg --gfa $minnow_ind/dbg.gfa -w $minnow_data/737K-august-2016.txt --countProb $minnow_data/hg/countProb_pbmc_4k.txt --custom |& stdbuf -oL tr '\r' '\n' > $log_file"
      echo $cmd
      eval $cmd
  done
done


##################################################################################################
# simulation with DGE
##################################################################################################
for sample in ${samples[@]}
do
  for type in ${celltypes[@]}
  do
      input_dir="04_simulation_minnow/simulation_DGE/normal$sample/$type"
      output_dir="04_simulation/simulation_DGE/normal${sample}_${type}_simulated_reads"
      log_file="04_simulation/simulation_DGE/minnow_simulate_normal${sample}_${type}.log"
      cmd="$minnow simulate --splatter-mode --g2t $g2t --inputdir $input_dir --PCR 6 -r $minnow_ind/ref_k101_fixed.fa -e 0.01 -p 20 -o $output_dir --dbg --gfa $minnow_ind/dbg.gfa -w $minnow_data/737K-august-2016.txt --countProb $minnow_data/hg/countProb_pbmc_4k.txt --custom |& stdbuf -oL tr '\r' '\n' > $log_file"
      echo $cmd
      eval $cmd
  done
done
