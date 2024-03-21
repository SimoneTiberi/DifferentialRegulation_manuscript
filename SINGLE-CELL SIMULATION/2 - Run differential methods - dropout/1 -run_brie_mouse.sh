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
# libcudart
# libnvinfer
# libnvinfer_plugin
# libcuda

#!/bin/bash
cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA/BRIE2_USA

celltypes=("Adipocytes" "Epithelial_cells" "Hepatocytes")

#########################################################################################
# drop95
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation_drop95/cell_info_${type}.h5ad \
    -c simulation_drop95/cell_info_${type}.tsv \
    -o simulation_drop95/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done

real    1590m7.194s
user    8410m42.364s
sys     155m23.081s

#########################################################################################
# drop99
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation_drop99/cell_info_${type}.h5ad \
    -c simulation_drop99/cell_info_${type}.tsv \
    -o simulation_drop99/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done

real    1631m1.530s
user    8625m48.381s
sys     155m5.039s
