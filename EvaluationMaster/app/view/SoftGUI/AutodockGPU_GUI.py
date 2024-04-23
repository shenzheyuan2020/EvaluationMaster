import sys
import subprocess
from PyQt5.QtWidgets import (QDialog, QApplication, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                             QFileDialog, QMessageBox, QWidget)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
from PyQt5.QtGui import QMovie, QPixmap
import os



class ScriptThread(QThread):
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

class autodockgpudialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

         # csv_file input (file chooser for CSV files)
        self.csvFileEdit = QLineEdit(self)
        self.csvFileEdit.setReadOnly(True)
        self.csvFileButton = QPushButton("Choose CSV File", self)
        self.csvFileButton.clicked.connect(self.chooseCsvFile)
        leftLayout.addWidget(QLabel("csv_file"))
        leftLayout.addWidget(self.csvFileEdit)
        leftLayout.addWidget(self.csvFileButton)
        
        # lig_path input
        self.ligFileEdit = QLineEdit(self)
        self.ligFileEdit.setReadOnly(True)
        self.ligFileButton = QPushButton("Choose Ligand File", self)
        self.ligFileButton.clicked.connect(self.chooseLigFile)
        leftLayout.addWidget(QLabel("lig_path"))
        leftLayout.addWidget(self.ligFileEdit)
        leftLayout.addWidget(self.ligFileButton)

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

        # save_path input
        self.saveDirEdit = QLineEdit(self)
        self.saveDirEdit.setReadOnly(True)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("save_path"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)

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

        # Save and Close buttons
        self.saveButton = QPushButton("AutodockGPU Docking", self)
        self.saveButton.clicked.connect(self.save_data)
        leftLayout.addWidget(self.saveButton)

        self.closeButton = QPushButton("Close", self)
        self.closeButton.clicked.connect(self.close)
        leftLayout.addWidget(self.closeButton)

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
        self.setFixedSize(600, 400)  # 可以根据需要调整尺寸

        # 设置样式
        self.setStyleSheet("""
            QDialog {
                background-color: rgba(255, 255, 255, 200);
            }
            QLabel, QLineEdit, QPushButton {
                /* 样式设置 */
            }
        """)

    # def chooseSavePlace(self):
    #     file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
    #     if file_name:
    #         self.savePlaceEdit.setText(file_name)

    # def chooseInputFile(self):
    #     file_name, _ = QFileDialog.getOpenFileName(self, "Choose Input File", "", "CSV Files (*.csv)")
    #     if file_name:
    #         self.inputFileEdit.setText(file_name)

    def chooseCsvFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.csvFileEdit.setText(file_name)

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

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)

    def save_data(self):
        # Gather inputs
        csv_file = self.csvFileEdit.text()
        lig_path = self.ligFileEdit.text()
        mgltool_path = self.mglToolPathEdit.text()
        gpu_file = self.gpuFileEdit.text()
        save_path = self.saveDirEdit.text()
        gpu_num = self.gpuNumEdit.text()
        n_run = self.nRunEdit.text()

        if csv_file and lig_path and mgltool_path and gpu_file and save_path and gpu_num and n_run:
            self.run_script(csv_file, lig_path, mgltool_path, gpu_file, save_path, int(gpu_num), int(n_run))
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields with valid data")

    def run_script(self, csv_file, lig_path, mgltool_path, gpu_file, save_path, gpu_num, n_run):
        self.thread = ScriptThread(csv_file, lig_path, mgltool_path, gpu_file, save_path, gpu_num, n_run)
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





if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = autodockgpudialogs()
    dialog.show()
    sys.exit(app.exec_())
