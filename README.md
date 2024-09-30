# EvaluationMaster

EvaluationMaster is a tool focused on chemical data processing and simulation, designed to assist researchers and developers in efficiently screening and evaluating compounds.

## Features

- **Data Processing**: Supports multiple data formats, providing powerful data processing capabilities.
- **Bacth Docking and Evaluation**: Integrates different docking tools to perform complex batch docking evaluation under multiple protein structures.
- **User Interface**: Offers a PyQt-based graphical user interface, easy to operate and visualize data.


## Updated notice （20240924）
We have fixed some errors and constructed an Install file for easier installation.

## Updated notice （20240623）
We discovered that Karmadock's code had been updated since our initial tests, particularly regarding execution scripts. To ensure smooth execution for the user. It is necessary to replace a script in the substitude folder with one in the same path as karmadock.

## Installation Steps
You can just run ./Install to complete installation or you can follow the steps:

- Installation Guide
- 1. Ensure Conda is Installed
- First, ensure that Anaconda or Miniconda is installed on your system. 

- 2. Initialize Conda Environment
- Run the following command in your terminal to initialize Conda shell commands:

bash
eval "$(conda shell.bash hook)"
- 3. Create and Activate Conda Virtual Environment
- Create a Conda virtual environment named "VM" with Python 3.9 and activate it:

bash
conda create -n VM python=3.9 -y
conda activate VM
- 4. Install PyQt and Other Dependencies
- While the virtual environment is activated, execute the following commands to install the required dependencies:

bash
conda install pyqt -y
wget https://files.pythonhosted.org/packages/b6/06/291866f91c573cc637bedbdd008e6a8ca506b589e61de1540ca265cfe7bd/PyQt_Fluent_Widgets-1.6.5-py3-none-any.whl
pip install PyQt_Fluent_Widgets-1.6.5-py3-none-any.whl
pip install qdarkstyle

conda install numpy -c conda-forge -y
conda install vina -c conda-forge -y
conda install rdkit -c conda-forge -y
conda install chembl_webresource_client matplotlib -c conda-forge -y
conda install pyarrow -c conda-forge -y
conda install openpyxl -c conda-forge -y
conda install scipy -c conda-forge -y
conda install scikit-learn==1.3.2 -c conda-forge -y
conda install openbabel -c conda-forge -y
conda install pandas -c conda-forge -y
conda install prody -c conda-forge -y
conda install wandb -c conda-forge -y
conda install scikit-image -c conda-forge -y
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia -y

pip install molgrid
pip install tqdm
pip install seaborn
- 5. Install Autogrid
- Execute the following command to install autogrid:

bash
sudo apt install autogrid
- 6. Install PyMOL
- Install the open-source version of PyMOL:

bash
conda install -c conda-forge pymol-open-source -y
- 7. Create Support_software Directory and Install MGLTools
- Create the Support_software directory and navigate into it, then install MGLTools:

bash
mkdir -p Support_software
cd Support_software

mkdir -p mgltools
cd mgltools
wget https://ccsb.scripps.edu/mgltools/download/491/ -O mgltools_x86_64Linux2_1.5.7p1.tar.gz
tar -zxvf mgltools_x86_64Linux2_1.5.7p1.tar.gz
cd mgltools_x86_64Linux2_1.5.7
bash install.sh
cd ../..
- 8. Install LeDock
- Create the Ledock directory and install LeDock:

bash
mkdir -p Ledock
cd Ledock
wget http://www.lephar.com/download/ledock_linux_x86
wget https://www.lephar.com/download/dok2mol2.cpp
g++ -std=c++11 dok2mol2.cpp -o dok2mol2
mv ledock_linux_x86 ledock
chmod 777 ledock
cd ..
- 9. Install KarmaDock
- Clone the KarmaDock repository and install its environment:

bash
git clone https://github.com/schrojunzhang/KarmaDock
cd KarmaDock
mkdir -p Env
cd Env
wget https://zenodo.org/record/7788732/files/karmadock_env.tar.gz?download=1 -O karmadock_env.tar.gz
tar -zxvf karmadock_env.tar.gz
cd ../..
- 10. Install AutoDock-GPU
- Clone the AutoDock-GPU repository and compile it:

bash
git clone https://github.com/ccsb-scripps/AutoDock-GPU
cd AutoDock-GPU
export GPU_INCLUDE_PATH=/usr/local/cuda/include
export GPU_LIBRARY_PATH=/usr/local/cuda/lib64
make DEVICE=OCLGPU NUMWI=64
cd ..
- 11. Install Deeppocket
- Clone the fpocket and DeepPocket repositories and install dependencies:

bash
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
sudo make install
cd ..
git clone https://github.com/devalab/DeepPocket.git
cd DeepPocket
wget https://zenodo.org/records/13833813/files/first_model_fold1_best_test_auc_85001.pth.tar
wget https://zenodo.org/records/13833813/files/seg0_best_test_IOU_91.pth.tar
cd ..
- 12. Install DeepCoy
- Clone the DeepCoy repository and install dependencies:

bash
git clone https://github.com/fimrie/DeepCoy.git
cd DeepCoy
wget https://zenodo.org/records/13831241/files/Deep_Coy_env.tar.gz
tar -zxvf Deep_Coy_env.tar.gz
cd ..
Download DeepCoy pretrained models:

bash
wget https://opig.stats.ox.ac.uk/data/downloads/DeepCoy_pretrained_models.tar.gz
tar -zxvf DeepCoy_pretrained_models.tar.gz
cd ..
- 13. Set EVALUATIONMASTER Environment Variable
- Add the installation path to .bashrc and set the EVALUATIONMASTER environment variable:

bash
echo "export EVALUATIONMASTER=$(pwd)" >> ~/.bashrc
source ~/.bashrc
- 14. Copy File
- Copy the necessary file to the target directory:

bash
cp $EVALUATIONMASTER/Substitude/KarmaDock/dataset/graph_obj.py $EVALUATIONMASTER/Support_software/KarmaDock/dataset/graph_obj.py
- 15. Complete Installation
- After completing the above steps, the installation process is finished. You can now use the installed software.
