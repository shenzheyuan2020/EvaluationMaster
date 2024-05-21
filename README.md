# EvaluationMaster

EvaluationMaster is a tool focused on chemical data processing and simulation, designed to assist researchers and developers in efficiently screening and evaluating compounds.

## Features

- **Data Processing**: Supports multiple data formats, providing powerful data processing capabilities.
- **Chemical Simulation**: Integrates multiple cheminformatics tools to perform complex chemical simulations.
- **User Interface**: Offers a PyQt-based graphical user interface, easy to operate and visualize data.

## Installation Steps

This project depends on multiple Python libraries and tools. Here are the installation steps:


```bash
conda create -n VM python=3.9
conda activate VM
conda install pyqt
pip install PyQt_Fluent_Widgets-1.4.6-py3-none-any.whl
pip install qdarkstyle
conda install numpy -c conda-forge
conda install vina -c conda-forge
conda install rdkit -c conda-forge
conda install chembl_webresource_client matplotlib -c conda-forge
conda install pyarrow -c conda-forge
conda install openpyxl -c conda-forge
conda install scipy -c conda-forge
conda install scikit-learn==1.3.2 -c conda-forge
conda install openbabel -c conda-forge
conda install pandas -c conda-forge
```
![image](https://github.com/shenzheyuan2020/EvaluationMaster/assets/73147896/283104e4-b812-4c51-bb4c-e520116ca0ac)

![image](https://github.com/shenzheyuan2020/EvaluationMaster/assets/73147896/e1cf75a5-0499-480f-bc80-addfab20a919)


## SupportSoftware could be download here
- KarmaDock github link  ——  https://github.com/schrojunzhang/KarmaDock
- AutoDock-GPU github link ——  https://github.com/ccsb-scripps/AutoDock-GPU
- LeDock download link  ——  http://www.lephar.com/download.htm
- mgltools download link  —— https://ccsb.scripps.edu/mgltools/downloads/ 
- DeepCoy github link  ——  https://github.com/fimrie/DeepCoy
- DeepPocket github link  ——    https://github.com/devalab/DeepPocket


## Details
    mkdir Support_software
    cd Support_software
    
    #For mgltools(20240521):
    mkdir mgltools
    cd mgltools
    wget https://ccsb.scripps.edu/mgltools/download/491/
    tar -zxvf mgltools_x86_64Linux2_1.5.7p1.tar.gz
    cd ..
    
    #For LeDock:
    mkdir Ledock
    cd Ledock
    wget http://www.lephar.com/download/ledock_linux_x86
    cd ..
    
    #For KarmaDock
    git clone https://github.com/schrojunzhang/KarmaDock

    #For AutodockGPU
    git clone https://github.com/ccsb-scripps/AutoDock-GPU
    cd AutoDock-GPU
    make DEVICE=GPU NUMWI=64
