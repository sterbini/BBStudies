                                              

# installing conda
mkdir ./Executables
if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform
    if [ "$(uname -m)" == "x86_64" ]; then
        wget -O ./Executables/Miniconda3-latest.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    elif [ "$(uname -m)" == "arm64" ]; then
        wget -O ./Executables/Miniconda3-latest.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
    fi
elif [ "$(uname)" == "Linux" ]; then
    # Do something under Linux platform
    wget -O ./Executables/Miniconda3-latest.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
fi
bash ./Executables/Miniconda3-latest.sh -b  -p ./Executables/miniconda -f


# create your own virtual environment in a new folder
source ./Executables/miniconda/bin/activate
#python -m venv ./Executables/py_BB
#source ./Executables/py_BB/bin/activate


# Install generic python packages
#========================================
pip install jupyterlab
pip install ipywidgets
pip install PyYAML
pip install pyarrow
pip install pandas
pip install matplotlib
pip install scipy
pip install ipympl
pip install ruamel.yaml
pip install rich
pip install lfm
pip install pynaff
pip install NAFFlib

# Adding the jupyter kernel to the list of kernels
python -m ipykernel install --user --name py_BB --display-name "py_BB"
#========================================


# Install CERN packages
#=========================================
pip install cpymad

git clone https://github.com/lhcopt/lhcmask.git ./Executables/miniconda/lhcmask
pip install -e ./Executables/miniconda/lhcmask

git clone https://github.com/xsuite/xobjects ./Executables/miniconda/xobjects
pip install -e ./Executables/miniconda/xobjects

git clone https://github.com/xsuite/xdeps ./Executables/miniconda/xdeps
pip install -e ./Executables/miniconda/xdeps

git clone https://github.com/xsuite/xpart ./Executables/miniconda/xpart
pip install -e ./Executables/miniconda/xpart

git clone https://github.com/xsuite/xtrack ./Executables/miniconda/xtrack
pip install -e ./Executables/miniconda/xtrack

git clone https://github.com/xsuite/xfields ./Executables/miniconda/xfields
pip install -e ./Executables/miniconda/xfields

git clone https://github.com/PyCOMPLETE/FillingPatterns.git ./Executables/miniconda/FillingPatterns
pip install -e ./Executables/miniconda/FillingPatterns

git clone https://github.com/xsuite/tree_maker.git ./Executables/miniconda/tree_maker
python -m pip install -e ./Executables/miniconda/tree_maker

git clone https://github.com/xsuite/xmask.git ./Executables/miniconda/xmask
python -m pip install -e ./Executables/miniconda/xmask

cd ./Executables/miniconda/xmask
git submodule init
git submodule update
cd ../../../

python -m pip install xsuite
xsuite-prebuild

#=========================================