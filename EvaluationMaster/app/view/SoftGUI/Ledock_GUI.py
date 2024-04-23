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

    def __init__(self, total_threads, csv_file, lig_path, save_path,ledock_path):
        super().__init__()
        # self.generate_number = generate_number
        # self.valid_number = valid_number
        self.total_threads = total_threads
        self.csv_file = csv_file
        self.lig_path = lig_path
        self.save_path = save_path
        self.ledock_path = ledock_path

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/DockingScript/Ledock.py"
        subprocess.call([f"python", script_path, self.total_threads, self.csv_file, self.lig_path, self.save_path, self.ledock_path])
        self.finished.emit(self.csv_file)  # Emit with parameters

class ledockdialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()

        # # generate_number input
        # self.generate_numberEdit = QLineEdit(self)
        # self.generate_numberEdit.setPlaceholderText("Enter generate number")
        # leftLayout.addWidget(QLabel("generate_number"))
        # leftLayout.addWidget(self.generate_numberEdit)

        # valid_number input
        # self.valid_numberEdit = QLineEdit(self)
        # self.valid_numberEdit.setPlaceholderText("Enter valid number")
        # leftLayout.addWidget(QLabel("valid_number"))
        # leftLayout.addWidget(self.valid_numberEdit)

        # input total_threads
        self.total_threadsEdit = QLineEdit(self)
        self.total_threadsEdit.setPlaceholderText("Enter total threads")
        leftLayout.addWidget(QLabel("total_threads"))
        leftLayout.addWidget(self.total_threadsEdit)


        # csv_file input (changed to file chooser for CSV files)
        self.csvFileEdit = QLineEdit(self)
        self.csvFileEdit.setReadOnly(True)
        self.csvFileButton = QPushButton("Choose CSV File", self)
        self.csvFileButton.clicked.connect(self.chooseCsvFile)
        leftLayout.addWidget(QLabel("csv_file"))
        leftLayout.addWidget(self.csvFileEdit)
        leftLayout.addWidget(self.csvFileButton)
        
        #lig_path input
        self.ligFileEdit = QLineEdit(self)
        self.ligFileEdit.setReadOnly(True)
        self.ligFileButton = QPushButton("Choose Ligand File path", self)
        self.ligFileButton.clicked.connect(self.chooseLigFile)
        leftLayout.addWidget(QLabel("lig_path"))
        leftLayout.addWidget(self.ligFileEdit)
        leftLayout.addWidget(self.ligFileButton)

        # save_path input
        self.saveDirEdit = QLineEdit(self)
        self.saveDirEdit.setReadOnly(True)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("save_path"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)

        # ledock_path input
        self.ledockFileEdit = QLineEdit(self)
        self.ledockFileEdit.setReadOnly(True)
        self.ledockFileButton = QPushButton("Choose Ledock Executable path", self)
        self.ledockFileButton.clicked.connect(self.chooseLedockFile)
        leftLayout.addWidget(QLabel("ledock_path"))
        leftLayout.addWidget(self.ledockFileEdit)
        leftLayout.addWidget(self.ledockFileButton)

        # Save and Close buttons
        self.saveButton = QPushButton("Ledock Docking", self)
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
    def chooseCsvFile(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.csvFileEdit.setText(file_name)
    def chooseLigFile(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Ligand File path")
        if directory:
            self.ligFileEdit.setText(directory)
    def chooseLedockFile(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Ledockpath")
        if directory:
            self.ledockFileEdit.setText(directory)

    def save_data(self):
        # generate_number = self.generate_numberEdit.text()
        # valid_number = self.valid_numberEdit.text()
        total_threads = self.total_threadsEdit.text()
        csv_file = self.csvFileEdit.text()
        lig_path = self.ligFileEdit.text()
        save_path = self.saveDirEdit.text()
        ledock_path = self.ledockFileEdit.text()
        if total_threads and csv_file and lig_path and save_path and ledock_path:
            self.run_script( total_threads,csv_file, lig_path, save_path,ledock_path)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self,total_threads, csv_file, lig_path, save_path, ledock_path):
        self.thread = ScriptThread(total_threads,csv_file, lig_path, save_path, ledock_path)
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

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = ledockdialogs()
    dialog.show()
    sys.exit(app.exec_())
