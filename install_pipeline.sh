# INSTALL MINIFORGE
wget https://github.com/conda-forge/miniforge/releases/download/24.11.3-0/Miniforge3-Linux-x86_64.sh

# If you want latest version
# wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc

# CLONE PIPELINE REPOSITORY
git clone https://github.com/nuriagaralon/cel-rnaseq.git
cd cel-rnaseq

# CREATE ENVIRONMENT
conda env create -f config/snake_env.yaml