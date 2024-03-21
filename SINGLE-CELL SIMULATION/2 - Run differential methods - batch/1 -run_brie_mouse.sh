# 1) SET CONDA ENVIRONMENT:
source ~/.bashrc
conda update conda

# set conda environment:
conda create -n TFProb python=3.7.9

conda env list
conda activate TFProb

# 2) install python 3.7 if missing (BRIE2 relies on python v 3.7)
# not sure if step 2) should be done before or after step 1).
# install new python version:
# pyenv install 3.7.9
# check python versions:
# pyenv versions
# choose python version:
# pyenv local 3.7.9
# pyenv versions

# 3) BRIE2 if missing
# install brie:
pip install -U brie

# 4) install tensorflow if missing (not necessary, but it gives a 10 speed-up):
pip install -U tensorflow
# conda install -c anaconda cudatoolkit

# this file specified the version to be used to 3.8.2!
# I changed it to 3.7.9
# nano /home/Shared_sherborne/simone/Diff_Velo/.python-version

# 5) run analyses (simulation_batch has open permissions to write in):
#!/bin/bash
cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA/BRIE2_USA

celltypes=("Adipocytes" "Epithelial_cells" "Hepatocytes")

#########################################################################################
# batch
#########################################################################################
time for type in ${celltypes[@]}
  do
    brie-quant -i simulation_batch/cell_info_${type}.h5ad \
    -c simulation_batch/cell_info_${type}.tsv \
    -o simulation_batch/isA_${type}.h5ad \
    --layers=spliced,unspliced,ambiguous \
    --minCell 0 --minCount 0 --minUniqCount 0 --minMIF 0 \
    --interceptMode gene --LRTindex 0
  done
