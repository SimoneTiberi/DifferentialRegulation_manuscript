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
cd /home/Shared_sherborne/simone/Diff_Velo/discovery/BRIE2_USA

celltypes=("CPNs" "Cycling" "Immature_CPNs" "Immature_PNs" "RG" "oRG" )

#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i ${type}.h5ad \
    -c ${type}.tsv \
    -o isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done
  
# store time:



