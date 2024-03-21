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
# FC 6:
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation_FC6/cell_info_${type}.h5ad \
    -c simulation_FC6/cell_info_${type}.tsv \
    -o simulation_FC6/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done

real    1417m10.365s
user    7464m14.732s
sys     128m16.160s

#########################################################################################
# FC 9:
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation_FC9/cell_info_${type}.h5ad \
    -c simulation_FC9/cell_info_${type}.tsv \
    -o simulation_FC9/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done

real    1439m36.248s
user    7607m36.250s
sys     130m52.843s
