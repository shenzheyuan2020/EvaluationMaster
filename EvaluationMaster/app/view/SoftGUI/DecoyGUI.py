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

    def __init__(self, generate_number, valid_number , csv_file, tool_path, save_dir):
        super().__init__()
        self.generate_number = generate_number
        self.valid_number = valid_number
        self.csv_file = csv_file
        self.tool_path = tool_path
        self.save_dir = save_dir

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/LigScript/decoy.py"
        subprocess.call([f"{self.tool_path}/Deepcoy/bin/python", script_path, self.generate_number, self.valid_number, self.csv_file, self.tool_path, self.save_dir])
        self.finished.emit(self.generate_number, self.csv_file)  # Emit with parameters

class decoydialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # generate_number input
        self.generate_numberEdit = QLineEdit(self)
        self.generate_numberEdit.setPlaceholderText("Enter generate number")
        leftLayout.addWidget(QLabel("generate_number"))
        leftLayout.addWidget(self.generate_numberEdit)

        # valid_number input
        self.valid_numberEdit = QLineEdit(self)
        self.valid_numberEdit.setPlaceholderText("Enter valid number")
        leftLayout.addWidget(QLabel("valid_number"))
        leftLayout.addWidget(self.valid_numberEdit)

        # csv_file input (changed to file chooser for CSV files)
        self.csvFileEdit = QLineEdit(self)
        self.csvFileEdit.setReadOnly(True)
        self.csvFileButton = QPushButton("Choose CSV File", self)
        self.csvFileButton.clicked.connect(self.chooseCsvFile)
        leftLayout.addWidget(QLabel("csv_file"))
        leftLayout.addWidget(self.csvFileEdit)
        leftLayout.addWidget(self.csvFileButton)

        # tool_path input
        self.DeepcoyDirEdit = QLineEdit(self)
        self.DeepcoyDirEdit.setReadOnly(True)
        self.DeepcoyDirButton = QPushButton("Choose Deepcoy Directory", self)
        self.DeepcoyDirButton.clicked.connect(self.chooseDeepcoyDir)
        leftLayout.addWidget(QLabel("tool_path"))
        leftLayout.addWidget(self.DeepcoyDirEdit)
        leftLayout.addWidget(self.DeepcoyDirButton)

        # save_dir input
        self.saveDirEdit = QLineEdit(self)
        self.saveDirEdit.setReadOnly(True)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("save_dir"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)

        # Save and Close buttons
        self.saveButton = QPushButton("Decoy_Generation", self)
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

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)

    def chooseDeepcoyDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Deepcoy Directory")
        if directory:
            self.DeepcoyDirEdit.setText(directory)
    def chooseCsvFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.csvFileEdit.setText(file_name)

    def save_data(self):
        generate_number = self.generate_numberEdit.text()
        valid_number = self.valid_numberEdit.text()
        save_place = self.csvFileEdit.text()
        tool_path = self.DeepcoyDirEdit.text()
        save_dir = self.saveDirEdit.text()
        if generate_number and valid_number and save_place and save_dir:
            self.run_script(generate_number, valid_number, save_place, tool_path, save_dir)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, generate_number, valid_number, save_place, tool_path, save_dir):
        self.thread = ScriptThread(generate_number, valid_number,  save_place, tool_path, save_dir)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, generate_number, csv_file):
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
    dialog = decoydialogs()
    dialog.show()
    sys.exit(app.exec_())
