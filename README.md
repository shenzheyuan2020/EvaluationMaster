# EvaluationMaster

EvaluationMaster 是一个专注于化学数据处理和模拟的工具，旨在帮助科研人员和开发者高效地进行化合物的筛选和评估。

## 功能特点

- **数据处理**: 支持多种数据格式，提供强大的数据处理能力。
- **化学模拟**: 集成多个化学信息学工具，进行复杂的化学模拟。
- **用户界面**: 提供基于PyQt的图形用户界面，易于操作和可视化数据。

## 安装步骤

本项目依赖于多个Python库和工具，以下是安装步骤：

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

![图片1](https://github.com/shenzheyuan2020/EvaluationMaster/assets/73147896/cb0ba3e5-7620-4733-b8dc-0dbf21337d8e)

KarmaDock github link  ——  https://github.com/schrojunzhang/KarmaDock
AutoDock-GPU github link ——  https://github.com/ccsb-scripps/AutoDock-GPU
LeDock download link  ——  http://www.lephar.com/download.htm
mgltools download link  —— https://ccsb.scripps.edu/mgltools/downloads/
DeepCoy github link  ——  https://github.com/fimrie/DeepCoy
DeepPocket github link  ——    https://github.com/devalab/DeepPocket
