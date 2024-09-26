###If you want to add a new docking script UI, you need to create a new class that inherits from QDialog, and then add the UI elements to the class.  Also, the ScriptThread should be created normally


import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QLabel, QPushButton, QCheckBox
from PyQt5.QtCore import Qt, QPoint

import subprocess
from PyQt5.QtWidgets import (QHBoxLayout,  QLineEdit, 
                             QFileDialog, QMessageBox, QComboBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
from PyQt5.QtGui import QMovie, QPixmap
import os


from PyQt5.QtWidgets import QDialog, QVBoxLayout 
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtGui import QMovie, QPixmap

class ScriptKarmadockThread(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str, str)

    def __init__(self, KarmaDock_dir, multi_protein_csv, lig_csv,  save_path):
        super().__init__()
        self.KarmaDock_dir = KarmaDock_dir
        self.multi_protein_csv = multi_protein_csv
        self.multi_protein_dir = os.path.dirname(multi_protein_csv)  # Update multi_protein_dir to the directory of multi_protein_csv
        self.lig_csv = lig_csv
        self.out_dir = save_path
        self.score_threshold = 0


    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = os.path.join(current_path, "app/Script/DockingScript/Karmadock.py")
        process = subprocess.Popen(
            [f"{self.KarmaDock_dir}/Env/bin/python", "-u", script_path, self.KarmaDock_dir, self.multi_protein_csv, self.multi_protein_dir, self.lig_csv, self.out_dir, str(self.score_threshold)],
            stdout=subprocess.PIPE,
            text=True)
        output, _ = process.communicate()
        print("DEBUG OUTPUT:")
        print(output)
        # 初始化变量，以防找不到RESULT行
        mol_name = "Not found"
        out_dir = "Not found"
        # 查找以"RESULT,"开头的行
        for line in output.strip().split("\n"):
            if line.startswith("RESULT,"):
                _, mol_name, out_dir = line.split(",")
                break  # 找到RESULT行后退出循环
        
        self.finished.emit(mol_name, out_dir)



class ScriptThreadadgpu(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, csv_file, lig_path, mgltool_path, gpu_file, save_path, gpu_num, n_run):
        super().__init__()
        self.csv_file = csv_file
        self.lig_path = lig_path
        self.mgltool_path = mgltool_path
        self.gpu_file = gpu_file
        self.save_path = save_path
        self.gpu_num = gpu_num
        self.n_run = n_run

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/DockingScript/AutodockGPU.py"
        subprocess.call([f"python", script_path, self.csv_file, self.lig_path, self.mgltool_path, self.gpu_file, self.save_path, str(self.gpu_num), str(self.n_run)])
        self.finished.emit(self.csv_file)  # Emit with parameter


class ScriptThreadledock(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, total_threads, csv_file, lig_path, save_path,ledock_path):
        super().__init__()
        self.total_threads = total_threads
        self.csv_file = csv_file
        self.lig_path = lig_path
        self.save_path = save_path
        self.ledock_path = ledock_path

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/DockingScript/LeDock.py"
        subprocess.call([f"python", script_path, self.total_threads, self.csv_file, self.lig_path, self.save_path, self.ledock_path])
        self.finished.emit(self.csv_file)  # Emit with parameters



class ScriptThreadvina(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, csv_file, lig_path, tool_path, save_path):
        super().__init__()

        self.csv_file = csv_file
        self.lig_path = lig_path
        self.tool_path = tool_path
        self.save_path = save_path

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/DockingScript/AutoDock_Vina.py"
        subprocess.call([f"python", script_path, self.csv_file, self.lig_path, self.tool_path, self.save_path])
        self.finished.emit(self.csv_file)  # Emit with parameters

# class ScriptThreadglide(QThread):
#     started = pyqtSignal()
#     finished = pyqtSignal(str)

#     def __init__(self, csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode,threads_num):
#         super().__init__()
#         self.csv_file_path = csv_file_path
#         self.lig_file = lig_file
#         self.schro_dir = schro_dir
#         self.mm_sahre_dir = mm_sahre_dir
#         self.save_dir = save_dir
#         self.Glide_mode = Glide_mode
#         self.threads_num = threads_num
#     def run(self):
#         self.started.emit()
#         current_path = os.getcwd()
#         script_path = current_path + "/app/Script/DockingScript/GlideG.py"
#         subprocess.call([f"python", script_path, self.csv_file_path, self.lig_file, self.schro_dir, self.mm_sahre_dir, self.save_dir, self.Glide_mode, self.threads_num])
#         self.finished.emit(self.csv_file_path)  # Emit with parameters


def setupAnimation(self, gif_path):
    label = QLabel(self)
    animation = QMovie(gif_path)
    label.setMovie(animation)
    animation.start()
    return label

class basicdialogs(QDialog):
    def __init__(self, parent=None):
        super(basicdialogs, self).__init__(parent)  # 保留parent参数以允许设置父窗口
        self.setWindowFlags(self.windowFlags() | Qt.Window | Qt.FramelessWindowHint)
        self.setWindowOpacity(1)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # 创建勾选框
        self.chkAutodockGPU = QCheckBox("Run AutodockGPU", self)
        self.chkLeDock = QCheckBox("Run LeDock", self)
        self.chkAutodockVina = QCheckBox("Run AutodockVina", self)
        self.chkKarmaodock = QCheckBox("Run Karmadock", self)


         # csv_file input (file chooser for CSV files)
        self.csvFileEdit = QLineEdit(self)
        self.csvFileEdit.setGeometry(QtCore.QRect(40, 340, 251, 25))
        self.csvFileEdit.setPlaceholderText("Choose the csv file which contain your protein name and coordinates information")
        self.csvFileButton = QPushButton("Choose CSV File", self)
        self.csvFileButton.clicked.connect(self.chooseCsvFile)
        leftLayout.addWidget(QLabel("protein csv file"))
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
        # 将勾选框添加到布局中
        leftLayout.addWidget(self.chkAutodockGPU)
        leftLayout.addWidget(self.chkLeDock)
        leftLayout.addWidget(self.chkAutodockVina)
        leftLayout.addWidget(self.chkKarmaodock)

        # 设置窗口的初始大小
        self.setFixedSize(800, 300)  # 可以根据需要调整尺寸

    def chooseCsvFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.csvFileEdit.setText(file_name)

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)


class autodockgpudialogs(QDialog):
    def __init__(self, parent=None):
        super(autodockgpudialogs, self).__init__(parent)  # 保留parent参数以允许设置父窗口
        self.setWindowFlags(self.windowFlags() | Qt.Window | Qt.FramelessWindowHint)
        self.setWindowOpacity(1)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()
      
        # lig_path input
        self.ligFileEdit = QLineEdit(self)
        self.ligFileEdit.setPlaceholderText("Choose the ligand file path (where you stored the prepared pdbqt file) for docking")
        self.ligFileButton = QPushButton("Choose Ligand file path", self)
        self.ligFileButton.clicked.connect(self.chooseLigFile)
        leftLayout.addWidget(QLabel("lig_path"))
        leftLayout.addWidget(self.ligFileEdit)
        leftLayout.addWidget(self.ligFileButton)

        # mgltool_path input
        self.mglToolPathEdit = QLineEdit(self)
        evaluation_master = os.getenv('EVALUATIONMASTER', '')  # Get environment variable
        self.mglToolPathEdit.setText(os.path.join(evaluation_master, "Support_software/mgltools/mgltools_x86_64Linux2_1.5.7")) 
        self.mglToolPathButton = QPushButton("Choose MGLTools Path", self)
        self.mglToolPathButton.clicked.connect(self.chooseMGLToolsPath)
        leftLayout.addWidget(QLabel("mgltool_path"))
        leftLayout.addWidget(self.mglToolPathEdit)
        leftLayout.addWidget(self.mglToolPathButton)

        # Autodock_GPU_file input
        self.gpuFileEdit = QLineEdit(self)
        evaluation_master = os.getenv('EVALUATIONMASTER', '')  # Get environment variable
        self.gpuFileEdit.setText(os.path.join(evaluation_master, "Support_software/AutoDock-GPU/bin/autodock_gpu_64wi")) 
        self.gpuFileButton = QPushButton("Choose AutoDock GPU File", self)
        self.gpuFileButton.clicked.connect(self.chooseGPUFile)
        leftLayout.addWidget(QLabel("Autodock_GPU_file"))
        leftLayout.addWidget(self.gpuFileEdit)
        leftLayout.addWidget(self.gpuFileButton)

        # GPU_num input
        self.gpuNumEdit = QLineEdit(self)
        self.gpuNumEdit.setPlaceholderText("Enter the GPU number for your use (recomended 1)")
        leftLayout.addWidget(QLabel("GPU_num"))
        leftLayout.addWidget(self.gpuNumEdit)

        # n_run input
        self.nRunEdit = QLineEdit(self)
        self.nRunEdit.setPlaceholderText("Enter number of runs (recomended 100)")
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
        super(ledockdialogs, self).__init__(parent)  # 保留parent参数以允许设置父窗口
        self.setWindowFlags(self.windowFlags() | Qt.Window | Qt.FramelessWindowHint)
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
        
        #lig_path input
        self.ligFileEdit = QLineEdit(self)
        self.ligFileEdit.setPlaceholderText("Choose the ligand file path (where you stored the prepared mol2 file) for docking")
        self.ligFileButton = QPushButton("Choose Ligand file path", self)
        self.ligFileButton.clicked.connect(self.chooseLigFile)
        leftLayout.addWidget(QLabel("lig_path"))
        leftLayout.addWidget(self.ligFileEdit)
        leftLayout.addWidget(self.ligFileButton)


        # ledock_path input
        self.ledockFileEdit = QLineEdit(self)
        evaluation_master = os.getenv('EVALUATIONMASTER', '')  # Get environment variable
        self.ledockFileEdit.setText(os.path.join(evaluation_master, "Support_software/Ledock")) 
        self.ledockFileButton = QPushButton("Choose Ledock Executable path", self)
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
        super(vinadialogs, self).__init__(parent)  # 保留parent参数以允许设置父窗口
        self.setWindowFlags(self.windowFlags() | Qt.Window | Qt.FramelessWindowHint)
        self.setWindowOpacity(1)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        #lig_path input
        self.ligFileEdit = QLineEdit(self)
        self.ligFileEdit.setPlaceholderText("Choose the ligand file path (where you stored the prepared pdbqt file) for docking")
        self.ligFileButton = QPushButton("Choose Ligand file path", self)
        self.ligFileButton.clicked.connect(self.chooseLigFile)
        leftLayout.addWidget(QLabel("lig_path"))
        leftLayout.addWidget(self.ligFileEdit)
        leftLayout.addWidget(self.ligFileButton)

        # tool_path input
        self.MGLToolsDirEdit = QLineEdit(self)
        evaluation_master = os.getenv('EVALUATIONMASTER', '')  # Get environment variable
        self.MGLToolsDirEdit.setText(os.path.join(evaluation_master, "Support_software/mgltools/mgltools_x86_64Linux2_1.5.7")) 
        self.MGLToolsDirButton = QPushButton("Choose MGLTools Path", self)
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


class karmadockdialogs(QDialog):
    def __init__(self, parent=None):
        super(karmadockdialogs, self).__init__(parent)  # 保留parent参数以允许设置父窗口
        self.setWindowFlags(self.windowFlags() | Qt.Window | Qt.FramelessWindowHint)
        self.setWindowOpacity(1)
        self.initUI()

    def initUI(self):
        # Main layout
        mainLayout = QHBoxLayout(self)

        # Left layout (Input area)
        leftLayout = QVBoxLayout()


        # UI elements for karmadock_dir
        self.karmadock_dir = QLineEdit(self)

        evaluation_master = os.getenv('EVALUATIONMASTER', '')  # Get environment variable
        self.karmadock_dir.setText(os.path.join(evaluation_master, "Support_software/KarmaDock/")) 
        self.Button_karmadockdir = QPushButton("KarmaDock dir", self)
        self.Button_karmadockdir.clicked.connect(self.choosekarmaDir)
        leftLayout.addWidget(QLabel("KarmaDock dir"))
        leftLayout.addWidget(self.karmadock_dir)
        leftLayout.addWidget(self.Button_karmadockdir)

        self.Line_pro_coor_file = QLineEdit(self)
        self.Line_pro_coor_file.setPlaceholderText("Choose coordinate csv file")        
        self.Button_pro_coor_file = QPushButton("Choose coordinate csv file", self)
        self.Button_pro_coor_file.clicked.connect(lambda: self.chooseFile(self.Line_pro_coor_file, "CSV Files (*.csv)"))
        leftLayout.addWidget(QLabel("Choose protein coordinate csv file"))
        leftLayout.addWidget(self.Line_pro_coor_file)
        leftLayout.addWidget(self.Button_pro_coor_file)


        self.Line_ligand_path = QLineEdit(self)
        self.Line_ligand_path.setPlaceholderText("Choose ligand file (.csv) for docking")
        self.Button_ligand = QPushButton("Choose ligand csv path", self)
        self.Button_ligand.clicked.connect(lambda: self.chooseFile(self.Line_ligand_path, "CSV Files (*.csv)"))
        leftLayout.addWidget(QLabel("Choose ligand csv file"))
        leftLayout.addWidget(self.Line_ligand_path)
        leftLayout.addWidget(self.Button_ligand)

        self.Line_scorethread = QLineEdit(self)
        self.Line_scorethread.setPlaceholderText("Enter score threshold")
        leftLayout.addWidget(QLabel("Set score threshold"))
        leftLayout.addWidget(self.Line_scorethread)

        # Left container
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)

        # Right layout (Animation area)
        rightLayout = QVBoxLayout()
        self.animationLabel = QLabel(self)
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        self.animationMovie = QMovie(movie_dir)  # Animation path
        self.animationLabel.setMovie(self.animationMovie)
        rightLayout.addWidget(self.animationLabel)

        # Right container
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)

        # Add both sides to the main layout
        mainLayout.addWidget(leftWidget)
        mainLayout.addWidget(rightWidget)

        # Set the window's initial size
        self.setFixedSize(800, 400)

    def choosekarmaDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.karmadock_dir.setText(dir)  

    # def chooseproteinDir(self):
    #     dir = QFileDialog.getExistingDirectory(self, "Select Directory")
    #     if dir:
    #         self.Line_protein_path.setText(dir)  

    def chooseFile(self, lineEdit, fileType):
        file, _ = QFileDialog.getOpenFileName(self, "Select File", "", fileType)
        if file:
            lineEdit.setText(file)



from PyQt5 import QtWidgets, QtCore

class TabbedGUI(QMainWindow):
    def __init__(self, parent=None):
        super(TabbedGUI, self).__init__(parent)  # 保留parent参数以允许设置父窗口
        self.setWindowFlags(self.windowFlags() | Qt.Window | Qt.FramelessWindowHint)
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
        self.tabs.addTab(self.createTab('Karmadock'), 'Karmadock')

        # 添加一个关闭按钮，并将其放置在窗口底部中央
        self.btnClose = QPushButton('Close', self)
        self.btnClose.clicked.connect(self.close)
        buttonWidth = 80
        buttonHeight = 30
        self.btnClose.setGeometry((self.width() - buttonWidth) / 2  + 60, self.height() - buttonHeight - 10, buttonWidth, buttonHeight)

        # Save and Close buttons
        self.saveButton = QPushButton("Docking", self)
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

        # 检查 AutodockGPU 是否被选中
        if self.basicdialogs.chkAutodockGPU.isChecked():
            lig_path = self.AutodockGPUDialog.ligFileEdit.text()
            mgl_tools_path = self.AutodockGPUDialog.mglToolPathEdit.text()
            gpu_file = self.AutodockGPUDialog.gpuFileEdit.text()
            gpu_num = self.AutodockGPUDialog.gpuNumEdit.text()
            nrun = self.AutodockGPUDialog.nRunEdit.text()
            # 创建 AutodockGPU 线程并加入队列
            threadADGPU = ScriptThreadadgpu(csv_file, lig_path, mgl_tools_path, gpu_file, save_path, gpu_num, nrun)
            self.enqueueThread(threadADGPU)

        # 检查 LeDock 是否被选中
        if self.basicdialogs.chkLeDock.isChecked():
            total_threads = self.ledockDialog.total_threadsEdit.text()
            lig_path = self.ledockDialog.ligFileEdit.text()
            ledock_path = self.ledockDialog.ledockFileEdit.text()
            # 创建 LeDock 线程并加入队列
            threadLeDock = ScriptThreadledock(total_threads, csv_file, lig_path, save_path, ledock_path)
            self.enqueueThread(threadLeDock)

        # 检查 AutodockVina 是否被选中
        if self.basicdialogs.chkAutodockVina.isChecked():
            lig_path = self.vinadialogs.ligFileEdit.text()
            mgl_tools_path = self.vinadialogs.MGLToolsDirEdit.text()
            # 创建 AutodockVina 线程并加入队列
            threadVina = ScriptThreadvina(csv_file, lig_path, mgl_tools_path, save_path)
            self.enqueueThread(threadVina)

        # 检查 Glide 是否被选中
        # if self.basicdialogs.chkKarmadock.isChecked():
        #     lig_file = self.glideDialog.ligFilePathEdit.text()
        #     schro_dir = self.glideDialog.schroDirEdit.text()
        #     mm_share_dir = self.glideDialog.mmShareDirEdit.text()
        #     glide_mode = self.glideDialog.glideModeComboBox.currentText()
        #     threads_num = self.glideDialog.threadsNumEdit.text()
        #     # 创建 Glide 线程并加入队列
        #     threadGlide = ScriptThreadglide(csv_file, lig_file, schro_dir, mm_share_dir, save_path, glide_mode, threads_num)
        #     self.enqueueThread(threadGlide)
    
        # 检查Karadock是否被选中
        if self.basicdialogs.chkKarmaodock.isChecked():
            KarmaDock_dir = self.karmadockdialogs.karmadock_dir.text()
            multi_protein_csv = self.karmadockdialogs.Line_pro_coor_file.text()
            lig_csv = self.karmadockdialogs.Line_ligand_path.text()
            #score_threshold = self.karmadockdialogs.Line_scorethread.text()
            # 创建 Glide 线程并加入队列
            threadKarmadock = ScriptKarmadockThread( KarmaDock_dir, multi_protein_csv, lig_csv, save_path)
            self.enqueueThread(threadKarmadock)

        # 开始执行队列中的第一个线程（如果有）
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

    # def run_glide_script(self, csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode, threads_num):
    #     # Create and start the ScriptThread with the collected parameters
    #     self.thread = ScriptThreadglide(csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode, threads_num)
    #     self.thread.started.connect(self.on_script_start)
    #     self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
    #     self.thread.start()

    def run_karmadock_script(self, KarmaDock_dir, multi_protein_csv, lig_csv,  save_path):
        self.thread = ScriptKarmadockThread(KarmaDock_dir, multi_protein_csv, lig_csv, save_path)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)
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
        # elif tool_name == "Glide":
        #     # 对于 Glide 标签，显示 glidedialogs 界面
        #     self.glidedialogs = glidedialogs()
        #     layout.addWidget(self.glidedialogs )
        elif tool_name == "Karmadock":
            # 对于 Glide 标签，显示 glidedialogs 界面
            self.karmadockdialogs = karmadockdialogs()
            layout.addWidget(self.karmadockdialogs)

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
    TabbedGUI.apply_stylesheet(app, stylesheet_path)
    main_gui = TabbedGUI()
    main_gui.show()
    sys.exit(app.exec_())
