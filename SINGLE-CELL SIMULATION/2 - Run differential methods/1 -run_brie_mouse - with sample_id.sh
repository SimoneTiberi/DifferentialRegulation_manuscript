# set conda environment:
conda create -n TFProb python=3.7
# yes
# install brie:
# pip install -U brie
# Successfully installed brie-2.2.2

# tensorflow (10* speed-up):
# pip install -U tensorflow
# conda install -c anaconda cudatoolkit

# missing libraries for tensorflow:
#libcudart
#libnvinfer
#libnvinfer_plugin
#libcuda

#!/bin/bash
cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA/BRIE2_USA_sample

celltypes=("Adipocytes" "Epithelial_cells" "Hepatocytes")

#########################################################################################
# NON-DGE:
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation/cell_info_${type}.h5ad \
    -c simulation/cell_info_${type}.tsv \
    -o simulation/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done

#########################################################################################
# DGE:
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation_DGE/cell_info_${type}.h5ad \
    -c simulation_DGE/cell_info_${type}.tsv \
    -o simulation_DGE/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done
  
