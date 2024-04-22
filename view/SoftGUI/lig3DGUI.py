import sys
import subprocess
from PyQt5.QtWidgets import (QDialog, QApplication, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                             QFileDialog, QMessageBox, QWidget)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
from PyQt5.QtGui import QMovie, QPixmap
import os



class ScriptThread(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str, str)

    def __init__(self, lig_start_name, csv_file_path , save_dir, mgltool_path):
        super().__init__()
        self.lig_start_name = lig_start_name
        self.csv_file_path = csv_file_path
        self.save_dir = save_dir
        self.mgltool_path = mgltool_path

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/LigScript/Gen3D.py"
        subprocess.call(["python", script_path, self.lig_start_name, self.csv_file_path, self.save_dir, self.mgltool_path])
        self.finished.emit(self.lig_start_name, self.save_dir)  # Emit with parameters

class Input3DligDialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # lig_start_name input
        self.lig_start_nameEdit = QLineEdit(self)
        self.lig_start_nameEdit.setPlaceholderText("Enter lig_start_name, eg: inh and you will get inh_1 to inh_n")
        leftLayout.addWidget(QLabel("lig_start_name"))
        leftLayout.addWidget(self.lig_start_nameEdit)

        # csv_file_path input
        self.inputFileEdit = QLineEdit(self)
        self.inputFileEdit.setReadOnly(True)
        self.inputFileButton = QPushButton("Choose Input File", self)
        self.inputFileButton.clicked.connect(self.chooseInputFile)
        leftLayout.addWidget(QLabel("csv_file_path"))
        leftLayout.addWidget(self.inputFileEdit)
        leftLayout.addWidget(self.inputFileButton)

        # save_dir input
        self.saveDirEdit = QLineEdit(self)
        self.saveDirEdit.setReadOnly(True)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("save_dir"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)

        # mgltool_path input
        self.MGLTOOLDirEdit = QLineEdit(self)
        self.MGLTOOLDirEdit.setReadOnly(True)
        self.MGLTOOLDirButton = QPushButton("Choose MGLTOOLS Directory", self)
        self.MGLTOOLDirButton.clicked.connect(self.chooseSchroDir)
        leftLayout.addWidget(QLabel("mgltool_path"))
        leftLayout.addWidget(self.MGLTOOLDirEdit)
        leftLayout.addWidget(self.MGLTOOLDirButton)

        # Save and Close buttons
        self.saveButton = QPushButton("Lig_3Dgen", self)
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

    def chooseInputFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose Input File", "", "CSV Files (*.csv)")
        if file_name:
            self.inputFileEdit.setText(file_name)

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)

    def chooseSchroDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Schro Directory")
        if directory:
            self.MGLTOOLDirEdit.setText(directory)

    def save_data(self):
        lig_start_name = self.lig_start_nameEdit.text()
        csv_file_path = self.inputFileEdit.text()
        save_place = self.saveDirEdit.text()
        mgltool_path = self.MGLTOOLDirEdit.text()

        if lig_start_name and csv_file_path and save_place:
            self.run_script(lig_start_name, csv_file_path, save_place, mgltool_path)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, lig_start_name, csv_file_path,  save_place, mgltool_path):
        self.thread = ScriptThread(lig_start_name, csv_file_path,  save_place, mgltool_path)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, lig_start_name, save_dir):
        self.animationMovie.stop()
        # current_path = os.getcwd()
        hand = os.path.splitext(save_dir)[0]
        current_path = os.getcwd()
        finishmap_path = f"{current_path}/images/logo.png"
        pixmap = QPixmap(finishmap_path)  # 结束时显示的图片
        scaled_pixmap = pixmap.scaled(self.animationLabel.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)  # 缩放图片
        self.animationLabel.setPixmap(scaled_pixmap)  # 显示缩放后的图片

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = Input3DligDialogs()
    dialog.show()
    sys.exit(app.exec_())
