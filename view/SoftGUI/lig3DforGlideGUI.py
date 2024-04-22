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

    def __init__(self, njobs, input_file , save_dir, Schro_dir):
        super().__init__()
        self.njobs = njobs
        self.input_file = input_file
        self.save_dir = save_dir
        self.Schro_dir = Schro_dir


    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/LigScript/ligprepG.py"
        subprocess.call(["python", script_path, self.njobs, self.input_file, self.save_dir, self.Schro_dir])
        self.finished.emit(self.njobs, self.save_dir)  # Emit with parameters

class Input3DligGDialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # njobs input
        self.njobsEdit = QLineEdit(self)
        self.njobsEdit.setPlaceholderText("Enter number of jobs")
        leftLayout.addWidget(QLabel("njobs"))
        leftLayout.addWidget(self.njobsEdit)

        # input_file input
        self.inputFileEdit = QLineEdit(self)
        self.inputFileEdit.setReadOnly(True)
        self.inputFileButton = QPushButton("Choose Input File", self)
        self.inputFileButton.clicked.connect(self.chooseInputFile)
        leftLayout.addWidget(QLabel("input_file"))
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

        # Schro_dir input
        self.SchroDirEdit = QLineEdit(self)
        self.SchroDirEdit.setReadOnly(True)
        self.SchroDirButton = QPushButton("Choose Schro Directory", self)
        self.SchroDirButton.clicked.connect(self.chooseSchroDir)
        leftLayout.addWidget(QLabel("Schro_dir"))
        leftLayout.addWidget(self.SchroDirEdit)
        leftLayout.addWidget(self.SchroDirButton)

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
            self.SchroDirEdit.setText(directory)

    def save_data(self):
        njobs = self.njobsEdit.text()
        input_file = self.inputFileEdit.text()
        save_place = self.saveDirEdit.text()
        Schro_dir = self.SchroDirEdit.text()

        if njobs and input_file and save_place:
            self.run_script(njobs, input_file, save_place, Schro_dir)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, njobs, input_file,  save_place, Schro_dir):
        self.thread = ScriptThread(njobs, input_file,  save_place, Schro_dir)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, njobs, save_dir):
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
    dialog = Input3DligGDialogs()
    dialog.show()
    sys.exit(app.exec_())
