#!/bin/bash
export DEBIAN_FRONTEND=noninteractive
export MPLCONFIGDIR=/root/.matplotlib
export CONFIG_FILE_PATH=/root/.jupyter/jupyter_notebook_config.py
mkdir -p $MPLCONFIGDIR

# Update, install and cleanup of system packages needed at run-time
apt-get update -y
apt-get upgrade -y
apt-get -y install build-essential wget ca-certificates git dos2unix zsh
apt-get -y install dvipng texlive-latex-extra texlive-fonts-recommended cm-super
rm -rf /var/lib/apt/lists/*
apt-get clean

# Install Miniforge3
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
/bin/bash ./Miniforge3-Linux-x86_64.sh -bfp /usr/local
conda config --file /.condarc --add channels defaults
conda config --file /.condarc --add channels conda-forge
conda config --file /.condarc --add channels bioconda
conda config --file /.condarc --add channels salilab
conda update conda -y

# Install Miniforge3 python packages
conda install -c bioconda -c conda-forge -c salilab gromacs_py gromacs==2024.4=nompi_cuda_h5cb645a_0 mdanalysis==2.8.0 scienceplots ambertools acpype openbabel pdb2pqr nglview==3.1.4 mdtraj seaborn dssp ipyparallel -y

# Download the latest CHARMM force field and replace the outdated 2017 one
wget "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz" -O charmm36-jul2022.ff.tgz

# Remove the outdated CHARMM force field
rm -rf /usr/local/lib/python3.12/site-packages/gromacs_py/gmx/template/charmm36-jul2017.ff
# Extract the downloaded tarball
tar -xzf charmm36-jul2022.ff.tgz
# Move the extracted CHARMM force field to the appropriate directory. This is done to fool the buggy gromacs_py module
mv charmm36-jul2022.ff /usr/local/lib/python3.12/site-packages/gromacs_py/gmx/template/charmm36-jul2017.ff
# Clean up
rm charmm36-jul2022.ff.tgz Miniforge3-Linux-x86_64.sh

#Create authentication certificates, self-signed
mkdir -p /root/.jupyter/certificates
openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout /root/.jupyter/certificates/jupyter.key -out /root/.jupyter/certificates/jupyter.pem

# Define the custom path for the Jupyter configuration file
CONFIG_FILE_PATH="/root/.jupyter/jupyter_notebook_config.py"
# Generate Jupyter configuration file at the custom path
jupyter notebook --generate-config -y --config="$CONFIG_FILE_PATH"

# Update the configuration file
echo "c.NotebookApp.ip = '0.0.0.0'" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.port = 8888" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.open_browser = False" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.token = ''" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.allow_origin = '*'" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.disable_check_xsrf = True " >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.allow_root=True" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.password_required=True" >> "$CONFIG_FILE_PATH" 
# Set the password on any python machine by getting the hash and pasting it below
# >>> from notebook.auth import passwd
# >>> passwd()
echo "c.NotebookApp.password = u'<paste hash here>'" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.certfile = u'/root/.jupyter/certificates/jupyter.pem'" >> "$CONFIG_FILE_PATH"
echo "c.NotebookApp.keyfile = u'/root/.jupyter/certificates/jupyter.key'" >> "$CONFIG_FILE_PATH"