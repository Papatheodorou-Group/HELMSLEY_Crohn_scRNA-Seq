source ~/.bashrc
conda create -n cellchat-env python=3.9
conda activate cellchat-env
conda install -c conda-forge numpy
conda install -c conda-forge umap-learn
conda deactivate
