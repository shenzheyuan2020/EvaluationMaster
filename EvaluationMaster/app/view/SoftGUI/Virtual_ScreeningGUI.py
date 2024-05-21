import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QLabel, QPushButton, QCheckBox
from PyQt5.QtCore import Qt, QPoint
import shutil
import subprocess
from PyQt5.QtWidgets import (QHBoxLayout,  QLineEdit, 
                             QFileDialog, QMessageBox, QComboBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
from PyQt5.QtGui import QMovie, QPixmap
import os
import pandas as pd

from PyQt5.QtWidgets import QDialog, QVBoxLayout 
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtGui import QMovie, QPixmap
from ...Script.ScreeningScript.file_conversion_script import convert_file_to_csv, process_csv_file

from PyQt5.QtWidgets import QSpacerItem, QSizePolicy


from .utilsGUI import *
class basicdialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(1)
        self.setStyleSheet("QDialog { margin-top: 0px; }")
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()
        leftLayout.setContentsMargins(15, 5, 10, 10)  # 设置左、上（减小这个值）、右、下边距

        leftLayout.setSpacing(0)
        leftLayout.setContentsMargins(15, 0, 10, 10)  # 设置左、上、右、下边距
        # 创建勾选框
        self.chkStage1 = QCheckBox("Stage1", self)
        self.chkStage2 = QCheckBox("Stage2", self)
        self.chkStage3 = QCheckBox("Stage3", self)
        self.chkStage4 = QCheckBox("Stage4", self)

        # For Stage1 - 已存在的代码
        self.cmbStage1 = QComboBox(self)
        self.cmbStage1.addItems(["AutodockGPU", "LeDock", "AutodockVina", "KarmaDock"])


        # For Stage2
        self.cmbStage2 = QComboBox(self)
        self.cmbStage2.addItems(["AutodockGPU", "LeDock", "AutodockVina", "KarmaDock"])
        self.numInputStage2 = QLineEdit(self)



        # For Stage3
        self.cmbStage3 = QComboBox(self)
        self.cmbStage3.addItems(["AutodockGPU", "LeDock", "AutodockVina", "KarmaDock"])
        self.numInputStage3 = QLineEdit(self)


        
        # For Stage4
        self.cmbStage4 = QComboBox(self)
        self.cmbStage4.addItems(["AutodockGPU", "LeDock", "AutodockVina", "KarmaDock"])
        self.numInputStage4 = QLineEdit(self)



        # 在类定义中添加新的元素
        self.ligandCsvFileEdit = QLineEdit(self)
        self.ligandCsvFileEdit.setReadOnly(True)
        self.ligandCsvFileButton = QPushButton("Choose Ligand CSV File", self)
        self.ligandCsvFileButton.clicked.connect(self.chooseLigandCsvFile)
        # 添加新的配体CSV文件选择器
        leftLayout.addWidget(QLabel("ligand_csv_file"))
        leftLayout.addWidget(self.ligandCsvFileEdit)
        leftLayout.addWidget(self.ligandCsvFileButton)

         # csv_file input (file chooser for CSV files)
        self.csvFileEdit = QLineEdit(self)
        self.csvFileEdit.setReadOnly(True)
        self.csvFileButton = QPushButton("Choose CSV File", self)
        self.csvFileButton.clicked.connect(self.chooseCsvFile)
        leftLayout.addWidget(QLabel("protein_csv_file"))
        leftLayout.addWidget(self.csvFileEdit)
        leftLayout.addWidget(self.csvFileButton)
        
        # save_path input
        self.saveDirEdit = QLineEdit(self)
        self.saveDirEdit.setReadOnly(True)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("save_path"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)


        # 左侧容器
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # 右侧布局（动画区）
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # 动画路径
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # 右侧容器
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # 将左右两侧添加到主布局
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)

        # 添加勾选框和下拉选项框到布局中
        leftLayout.addWidget(self.chkStage1)
        leftLayout.addWidget(self.cmbStage1)

        leftLayout.addWidget(self.chkStage2)
        leftLayout.addWidget(self.cmbStage2)
        leftLayout.addWidget(self.numInputStage2)  # 添加Stage2数字输入
        leftLayout.addWidget(self.chkStage3)
        leftLayout.addWidget(self.cmbStage3)

        leftLayout.addWidget(self.numInputStage3)  # 添加Stage3数字输入
        leftLayout.addWidget(self.chkStage4)
        leftLayout.addWidget(self.cmbStage4)
        leftLayout.addWidget(self.numInputStage4)  # 添加Stage2数字输入
        # 设置窗口的初始大小
        self.setFixedSize(800, 500)  # 可以根据需要调整尺寸

    def chooseCsvFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.csvFileEdit.setText(file_name)

    # 实现选择配体CSV文件的槽函数
    def chooseLigandCsvFile(self):
        fileName, _ = QFileDialog.getOpenFileName(self, "Open Ligand CSV File", "", "CSV Files (*.csv)")
        if fileName:
            self.ligandCsvFileEdit.setText(fileName)

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)


class autodockgpudialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(1)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # mgltool_path input
        self.mglToolPathEdit = QLineEdit(self)
        self.mglToolPathEdit.setReadOnly(True)
        self.mglToolPathButton = QPushButton("Choose MGLTools Path", self)
        self.mglToolPathButton.clicked.connect(self.chooseMGLToolsPath)
        leftLayout.addWidget(QLabel("mgltool_path"))
        leftLayout.addWidget(self.mglToolPathEdit)
        leftLayout.addWidget(self.mglToolPathButton)

        # Autodock_GPU_file input
        self.gpuFileEdit = QLineEdit(self)
        self.gpuFileEdit.setReadOnly(True)
        self.gpuFileButton = QPushButton("Choose AutoDock GPU File", self)
        self.gpuFileButton.clicked.connect(self.chooseGPUFile)
        leftLayout.addWidget(QLabel("Autodock_GPU_file"))
        leftLayout.addWidget(self.gpuFileEdit)
        leftLayout.addWidget(self.gpuFileButton)

        # GPU_num input
        self.gpuNumEdit = QLineEdit(self)
        self.gpuNumEdit.setPlaceholderText("Enter GPU number")
        leftLayout.addWidget(QLabel("GPU_num"))
        leftLayout.addWidget(self.gpuNumEdit)

        # n_run input
        self.nRunEdit = QLineEdit(self)
        self.nRunEdit.setPlaceholderText("Enter number of runs")
        leftLayout.addWidget(QLabel("n_run"))
        leftLayout.addWidget(self.nRunEdit)


        # 左侧容器
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # 右侧布局（动画区）
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # 动画路径
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # 右侧容器
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # 将左右两侧添加到主布局
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)

        # 设置窗口的初始大小
        self.setFixedSize(800, 400)  # 可以根据需要调整尺寸

    def chooseLigFile(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Ligand File")
        if directory:
            self.ligFileEdit.setText(directory)

    def chooseMGLToolsPath(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose MGLTools Path")
        if directory:
            self.mglToolPathEdit.setText(directory)

    def chooseGPUFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose AutoDock GPU File", "", "All Files (*)")
        if file_name:
            self.gpuFileEdit.setText(file_name)


class ledockdialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(1)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # input total_threads
        self.total_threadsEdit = QLineEdit(self)
        self.total_threadsEdit.setPlaceholderText("Enter total threads")
        leftLayout.addWidget(QLabel("total_threads"))
        leftLayout.addWidget(self.total_threadsEdit)


        # ledock_path input
        self.ledockFileEdit = QLineEdit(self)
        self.ledockFileEdit.setReadOnly(True)
        self.ledockFileButton = QPushButton("Choose LeDock Executable path", self)
        self.ledockFileButton.clicked.connect(self.chooseLeDockFile)
        leftLayout.addWidget(QLabel("ledock_path"))
        leftLayout.addWidget(self.ledockFileEdit)
        leftLayout.addWidget(self.ledockFileButton)

        
        # 左侧容器
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # 右侧布局（动画区）
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # 动画路径
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # 右侧容器
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # 将左右两侧添加到主布局
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)

        # 设置窗口的初始大小
        self.setFixedSize(800, 300)  # 可以根据需要调整尺寸

    def chooseLigFile(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Ligand File path")
        if directory:
            self.ligFileEdit.setText(directory)
    def chooseLeDockFile(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose LeDockpath")
        if directory:
            self.ledockFileEdit.setText(directory)


class vinadialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(1)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # tool_path input
        self.MGLToolsDirEdit = QLineEdit(self)
        self.MGLToolsDirEdit.setReadOnly(True)
        self.MGLToolsDirButton = QPushButton("Choose MGLTools Directory", self)
        self.MGLToolsDirButton.clicked.connect(self.chooseMGLToolsDir)
        leftLayout.addWidget(QLabel("tool_path"))
        leftLayout.addWidget(self.MGLToolsDirEdit)
        leftLayout.addWidget(self.MGLToolsDirButton)



        # 左侧容器
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # 右侧布局（动画区）
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # 动画路径
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # 右侧容器
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # 将左右两侧添加到主布局
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)

        # 设置窗口的初始大小
        self.setFixedSize(800, 200)  # 可以根据需要调整尺寸

    def chooseMGLToolsDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose MGLTools Directory")
        if directory:
            self.MGLToolsDirEdit.setText(directory)

    def chooseLigFile(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Ligand File path")
        if directory:
            self.ligFileEdit.setText(directory)



# class glidedialogs(QDialog):
#     def __init__(self, parent=None):
#         super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
#         self.setWindowOpacity(1)
#         self.initUI()

#     def initUI(self):
#         # Main layout
#         mainLayout = QHBoxLayout(self)

#         # Left layout (Input area)
#         leftLayout = QVBoxLayout()

#         # UI elements for schro_dir
#         self.schroDirEdit = QLineEdit(self)
#         self.schroDirButton = QPushButton("Choose Schrodir Directory", self)
#         self.schroDirButton.clicked.connect(self.chooseSchroDir)
#         leftLayout.addWidget(QLabel("Schrodinger Directory"))
#         leftLayout.addWidget(self.schroDirEdit)
#         leftLayout.addWidget(self.schroDirButton)

#         # UI elements for mm_sahre_dir
#         self.mmShareDirEdit = QLineEdit(self)
#         self.mmShareDirButton = QPushButton("Choose MMshare Directory", self)
#         self.mmShareDirButton.clicked.connect(self.chooseMMshareDir)
#         leftLayout.addWidget(QLabel("MMshare Directory"))
#         leftLayout.addWidget(self.mmShareDirEdit)
#         leftLayout.addWidget(self.mmShareDirButton)

#         # UI elements for KarmaDock_mode
#         self.glideModeComboBox = QComboBox(self)
#         self.glideModeComboBox.addItems(["SP", "XP"])
#         leftLayout.addWidget(QLabel("KarmaDock Mode"))
#         leftLayout.addWidget(self.glideModeComboBox)


#         # UI elements for threads_num
#         self.threadsNumEdit = QLineEdit(self)
#         self.threadsNumEdit.setPlaceholderText("Enter threads number")
#         leftLayout.addWidget(QLabel("Threads Number"))
#         leftLayout.addWidget(self.threadsNumEdit)                                               
      
#         # 左侧容器
#         leftWidget = QWidget()
#         leftWidget.setLayout(leftLayout)

#         # 右侧布局（动画区）
#         rightLayout = QVBoxLayout()
#         self.animationLabel = QLabel(self)
#         current_path = os.getcwd()
#         movie_dir = current_path + "/images/loading.gif"
#         self.animationMovie = QMovie(movie_dir)  # 动画路径
#         self.animationLabel.setMovie(self.animationMovie)
#         rightLayout.addWidget(self.animationLabel)

#         # 右侧容器
#         rightWidget = QWidget()
#         rightWidget.setLayout(rightLayout)

#         # 将左右两侧添加到主布局
#         mainLayout.addWidget(leftWidget)
#         mainLayout.addWidget(rightWidget)

#         # 设置窗口的初始大小
#         self.setFixedSize(800, 400)  # 可以根据需要调整尺寸


#     def chooseLigFilePath(self):
#         file_name, _ = QFileDialog.getOpenFileName(self, "Choose Ligand File", "", "All Files (*)")
#         if file_name:
#             self.ligFilePathEdit.setText(file_name)

#     def chooseSchroDir(self):
#         directory = QFileDialog.getExistingDirectory(self, "Choose Schrodinger Directory")
#         if directory:
#             self.schroDirEdit.setText(directory)

#     def chooseMMshareDir(self):
#         directory = QFileDialog.getExistingDirectory(self, "Choose MMshare Directory")
#         if directory:
#             self.mmShareDirEdit.setText(directory)










class glidedialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(1)
        self.initUI()

    def initUI(self):
        # Main layout
        mainLayout = QHBoxLayout(self)

        # Left layout (Input area)
        leftLayout = QVBoxLayout()

        # UI elements for schro_dir
        self.schroDirEdit = QLineEdit(self)
        self.schroDirButton = QPushButton("Choose KarmaDock Directory", self)
        self.schroDirButton.clicked.connect(self.chooseSchroDir)
        leftLayout.addWidget(QLabel("KarmaDock Directory"))
        leftLayout.addWidget(self.schroDirEdit)
        leftLayout.addWidget(self.schroDirButton)

        # # UI elements for mm_sahre_dir
        # self.mmShareDirEdit = QLineEdit(self)
        # self.mmShareDirButton = QPushButton("Choose MMshare Directory", self)
        # self.mmShareDirButton.clicked.connect(self.chooseMMshareDir)
        # leftLayout.addWidget(QLabel("MMshare Directory"))
        # leftLayout.addWidget(self.mmShareDirEdit)
        # leftLayout.addWidget(self.mmShareDirButton)

        # # UI elements for KarmaDock_mode
        # self.glideModeComboBox = QComboBox(self)
        # self.glideModeComboBox.addItems(["SP", "XP"])
        # leftLayout.addWidget(QLabel("KarmaDock Mode"))
        # leftLayout.addWidget(self.glideModeComboBox)


        # UI elements for threads_num
        self.threadsNumEdit = QLineEdit(self)
        self.threadsNumEdit.setPlaceholderText("Enter Score Threshold")
        leftLayout.addWidget(QLabel("Threshold"))
        leftLayout.addWidget(self.threadsNumEdit)                                               
      
        # 左侧容器
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # 右侧布局（动画区）
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # 动画路径
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # 右侧容器
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # 将左右两侧添加到主布局
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)

        # 设置窗口的初始大小
        self.setFixedSize(800, 400)  # 可以根据需要调整尺寸


    # def chooseLigFilePath(self):
    #     file_name, _ = QFileDialog.getOpenFileName(self, "Choose Ligand File", "", "All Files (*)")
    #     if file_name:
    #         self.ligFilePathEdit.setText(file_name)

    def chooseSchroDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Schrodinger Directory")
        if directory:
            self.schroDirEdit.setText(directory)

    # def chooseMMshareDir(self):
    #     directory = QFileDialog.getExistingDirectory(self, "Choose MMshare Directory")
    #     if directory:
    #         self.mmShareDirEdit.setText(directory)






class TabbedGUI2(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        self.threadQueue = []  # 线程队列初始化为类成员变量
        self.currentThread = None  # 当前正在运行的线程
        self.oldPos = self.pos()

    def enqueueThread(self, thread):
        """将线程添加到队列"""
        self.threadQueue.append(thread)

    def startNextThread(self):
        """开始执行队列中的下一个线程"""
        if self.threadQueue:
            self.currentThread = self.threadQueue.pop(0)  # 取出下一个线程
            self.currentThread.finished.connect(self.threadFinished)  # 连接完成信号到槽
            self.currentThread.start()
        else:
            self.currentThread = None  # 没有更多线程时清理当前线程引用
    def threadFinished(self):
        """当前线程完成时调用"""
        self.startNextThread()  # 尝试启动队列中的下一个线程


    def initUI(self):
        mainLayout = QHBoxLayout(self)

        self.setWindowTitle('Docking Tools Configuration')
        self.setGeometry(100, 100, 800, 600)
        self.setWindowFlags(Qt.FramelessWindowHint)  # 去掉窗口边框
        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()
        # 创建QTabWidget对象
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # 创建每个工具的标签页
        self.tabs.addTab(self.createTab('Basic settings'), 'Basic settings')
        self.tabs.addTab(self.createTab('AutodockGPU'), 'AutodockGPU')
        self.tabs.addTab(self.createTab('LeDock'), 'LeDock')
        self.tabs.addTab(self.createTab('AutodockVina'), 'AutodockVina')
        self.tabs.addTab(self.createTab('KarmaDock'), 'KarmaDock')

        # 添加一个关闭按钮，并将其放置在窗口底部中央
        self.btnClose = QPushButton('Close', self)
        self.btnClose.clicked.connect(self.close)
        buttonWidth = 80
        buttonHeight = 30
        self.btnClose.setGeometry((self.width() - buttonWidth) / 2  + 60, self.height() - buttonHeight - 10, buttonWidth, buttonHeight)

        # Save and Close buttons
        self.saveButton = QPushButton("Screening", self)
        self.saveButton.clicked.connect(self.save_data)
        buttonWidth = 80
        buttonHeight = 30
        self.saveButton.setGeometry((self.width() - buttonWidth) / 2 - 60, self.height() - buttonHeight - 10, buttonWidth, buttonHeight)
        
        # 左侧容器
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # 右侧布局（动画区）
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # 动画路径
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # 右侧容器
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # 将左右两侧添加到主布局
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)




    def save_data(self):
        csv_file = self.basicdialogs.csvFileEdit.text()
        save_path = self.basicdialogs.saveDirEdit.text()
        # ligand_csv_file = self.basicdialogs.ligandCsvFileEdit.text()
        stages = [
            (self.basicdialogs.chkStage1, self.basicdialogs.cmbStage1, '100'),
            (self.basicdialogs.chkStage2, self.basicdialogs.cmbStage2, self.basicdialogs.numInputStage2.text() if self.basicdialogs.numInputStage2 else '100'),
            (self.basicdialogs.chkStage3, self.basicdialogs.cmbStage3, self.basicdialogs.numInputStage3.text() if self.basicdialogs.numInputStage3 else '100'),
            (self.basicdialogs.chkStage4, self.basicdialogs.cmbStage4, self.basicdialogs.numInputStage4.text() if self.basicdialogs.numInputStage4 else '100'),
        ]
        temp_save_path = save_path
        for index, (chkStage, cmbStage, numInputValue) in enumerate(stages, start=1):
            if chkStage.isChecked():
                selectedTool = cmbStage.currentText()
                preserve_num = float(numInputValue)  # 确保preserve_num是浮点数以支持百分比计算

                # 创建ligands文件夹
                # lig_path = f"{save_path}/stage{index}/ligands/"

                # if not os.path.exists(lig_path):
                #     os.makedirs(lig_path)

                # 对于第2、3、4阶段，处理ligand_csv_file以只保留前preserve_num%的数据
                if index == 1:
                    # 创建ligands文件夹
                    lig_path = f"{save_path}/stage{index}/ligands/"

                    if not os.path.exists(lig_path):
                        os.makedirs(lig_path)

                    save_path = temp_save_path
                    print(f"Using {save_path} as the new save_path for stage 1")
                    orig_save_path = save_path ### Keep original save_path for later use
                    save_path = f"{save_path}/stage{index}/"  #### 第1阶段，Use save_path + "stage1" as the new save_path
                    ligand_csv_file = self.basicdialogs.ligandCsvFileEdit.text() ### It is the Key
                    # ligand_csv_basename = os.path.basename(ligand_csv_file)
                    # ligand_csv_file = os.path.dirname(ligand_csv_file)
                    backup_name = os.path.join(orig_save_path, "ref.csv" )
                    shutil.copy(ligand_csv_file, backup_name)
                    ligand_csv_file = process_csv_file(backup_name)

                    # ligand_csv_file = hand_ligand_csv_filename  # 更新ligand_csv_file为处理后的文件路径       
                # elif index > 1:  # 当在第2、3、4阶段时
                #     save_path = temp_save_path
                #     print(f"Using {save_path} as the new save_path for stage {index}")
                #     temp_csv_file = f"{save_path}/stage{index-1}/result.csv"
                #     save_path = f"{save_path}/stage{index}/"

                #     df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
                #     num_rows_to_keep = int(len(df) * (preserve_num / 100.0))
                #     df_subset = df.head(num_rows_to_keep)
                #     df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
                #     ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
                else:
                    save_path = temp_save_path                                                                                                                                                                                                             
                    # 创建ligands文件夹
                    lig_path = f"{save_path}/stage{index}/ligands/"

                    if not os.path.exists(lig_path):
                        os.makedirs(lig_path)
                   
                if selectedTool == "AutodockGPU":
                    mgl_tools_path = self.AutodockGPUDialog.mglToolPathEdit.text()
                    gpu_file = self.AutodockGPUDialog.gpuFileEdit.text()
                    gpu_num = self.AutodockGPUDialog.gpuNumEdit.text()
                    n_run = self.AutodockGPUDialog.nRunEdit.text()
                    threadInstance = ScriptThreadadgpu(csv_file, ligand_csv_file, lig_path, mgl_tools_path, gpu_file, save_path, gpu_num, n_run, index, preserve_num)

                elif selectedTool == "LeDock":
                    total_threads = self.ledockDialog.total_threadsEdit.text()
                    ledock_path = self.ledockDialog.ledockFileEdit.text()
                    threadInstance = ScriptThreadledock(total_threads, csv_file, ligand_csv_file, lig_path, save_path, ledock_path, index, preserve_num)

                elif selectedTool == "AutodockVina":
                    tool_path = self.vinadialogs.MGLToolsDirEdit.text()
                    threadInstance = ScriptThreadvina(csv_file, ligand_csv_file, lig_path, tool_path, save_path, index, preserve_num)

                # elif selectedTool == "KarmaDock":
                #     schro_dir = self.glideDialog.schroDirEdit.text()
                #     mm_share_dir = self.glideDialog.mmShareDirEdit.text()
                #     glide_mode = self.glideDialog.glideModeComboBox.currentText()
                #     threads_num = self.glideDialog.threadsNumEdit.text()
                #     threadInstance = ScriptThreadglide(csv_file, ligand_csv_file ,lig_path, schro_dir, mm_share_dir, save_path, glide_mode, threads_num, index, preserve_num)
                elif selectedTool == "KarmaDock":
                    schro_dir = self.glideDialog.schroDirEdit.text()
                    # mm_share_dir = self.glideDialog.mmShareDirEdit.text()
                    # glide_mode = self.glideDialog.glideModeComboBox.currentText()
                    threads_num = self.glideDialog.threadsNumEdit.text()
                    file_path = os.path.dirname(csv_file)
                    threadInstance = ScriptThreadglide(schro_dir, csv_file, file_path, ligand_csv_file, save_path, threads_num,  index, preserve_num)
                self.enqueueThread(threadInstance)

        # Begin executing the first thread if none are currently running
        if self.currentThread is None:
            self.startNextThread()

    def run_autodockgpu_script(self, csv_file, lig_path, mgltool_path, gpu_file, save_path, gpu_num, n_run):
        self.thread = ScriptThreadadgpu(csv_file, lig_path, mgltool_path, gpu_file, save_path, gpu_num, n_run)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)
        self.thread.start()

    def run_ledock_script(self,total_threads, csv_file, lig_path, save_path, ledock_path):
        self.thread = ScriptThreadledock(total_threads,csv_file, lig_path, save_path, ledock_path)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def run_vina_script(self,csv_file, lig_path, tool_path, save_path):
        self.thread = ScriptThreadvina(csv_file, lig_path, tool_path, save_path)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()

    # def run_glide_script(self, csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, KarmaDock_mode, threads_num):
    #     # Create and start the ScriptThread with the collected parameters
    #     self.thread = ScriptThreadglide(csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, KarmaDock_mode, threads_num)
    #     self.thread.started.connect(self.on_script_start)
    #     self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
    #     self.thread.start()

    def run_glide_script(self, schro_dir, csv_file, file_path, ligand_csv_file, save_path,threads_num,):
        # Create and start the ScriptThread with the collected parameters
        self.thread = ScriptThreadglide(schro_dir, csv_file, file_path, ligand_csv_file,save_path,threads_num)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, csv_file):
        self.animationMovie.stop()
        # current_path = os.getcwd()
        hand = os.path.splitext(csv_file)[0]
        current_path = os.getcwd()
        finishmap_path = f"{current_path}/images/logo.png"
        pixmap = QPixmap(finishmap_path)  # 结束时显示的图片
        scaled_pixmap = pixmap.scaled(self.animationLabel.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)  # 缩放图片
        self.animationLabel.setPixmap(scaled_pixmap)  # 显示缩放后的图片



    def createTab(self, tool_name):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        if tool_name == "Basic settings":
            self.basicdialogs=basicdialogs()
            layout.addWidget(self.basicdialogs)
        elif tool_name == "AutodockGPU":
            self.AutodockGPUDialog = autodockgpudialogs()
            layout.addWidget(self.AutodockGPUDialog)
        elif tool_name == "LeDock":
            self.ledockDialog = ledockdialogs()
            layout.addWidget(self.ledockDialog)
        elif tool_name == "AutodockVina": 
            self.vinadialogs = vinadialogs()
            layout.addWidget(self.vinadialogs)
        elif tool_name == "KarmaDock":
            # 对于 KarmaDock 标签，显示 glidedialogs 界面
            self.glideDialog = glidedialogs()
            layout.addWidget(self.glideDialog)

        return tab
    
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.oldPos = event.globalPos()

    def mouseMoveEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            delta = QPoint(event.globalPos() - self.oldPos)
            self.move(self.x() + delta.x(), self.y() + delta.y())
            self.oldPos = event.globalPos()

    def apply_stylesheet(app, stylesheet_path):
        with open(stylesheet_path, "r") as file:
            app.setStyleSheet(file.read())
if __name__ == "__main__":
    app = QApplication(sys.argv)
    stylesheet_path = "style.qss"  # 样式表文件的路径
    TabbedGUI2.apply_stylesheet(app, stylesheet_path)
    main_gui = TabbedGUI2()
    main_gui.show()
    sys.exit(app.exec_())
