import sys
import subprocess
from PyQt5.QtWidgets import (QDialog, QApplication, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                             QFileDialog, QMessageBox, QWidget, QComboBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
from PyQt5.QtGui import QMovie, QPixmap
import os



class ScriptThread(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, csv_file_path, save_dir, Glide_mode):
        super().__init__()
        self.csv_file_path = csv_file_path
        self.save_dir = save_dir
        self.Glide_mode = Glide_mode
    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path =f"/home/hoo/Install/Evaluation_X/Support_software/BK-score-model-main/{self.Glide_mode}.py"
        subprocess.call([f"python", script_path, self.csv_file_path, self.save_dir])
        self.finished.emit(self.csv_file_path)  # Emit with parameters

class BKdialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # Main layout
        mainLayout = QHBoxLayout(self)

        # Left layout (Input area)
        leftLayout = QVBoxLayout()

           # Add new UI elements for csv_file_path
        self.csvFilePathEdit = QLineEdit(self)
        self.csvFilePathButton = QPushButton("Choose CSV File Path", self)
        self.csvFilePathButton.clicked.connect(self.chooseCsvFilePath)
        leftLayout.addWidget(QLabel("CSV File Path"))
        leftLayout.addWidget(self.csvFilePathEdit)
        leftLayout.addWidget(self.csvFilePathButton)

        # # UI elements for save_dir
        # self.ligFilePathEdit = QLineEdit(self)
        # self.ligFilePathButton = QPushButton("Choose Ligand File", self)
        # self.ligFilePathButton.clicked.connect(self.chooseLigFilePath)
        # leftLayout.addWidget(QLabel("Ligand File"))
        # leftLayout.addWidget(self.ligFilePathEdit)
        # leftLayout.addWidget(self.ligFilePathButton)


        # UI elements for save_dir
        self.saveDirEdit = QLineEdit(self)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("Save Directory"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)


        # UI elements for Glide_mode
        self.glideModeComboBox = QComboBox(self)
        self.glideModeComboBox.addItems(["B_SCORE", "K_SCORE"])
        leftLayout.addWidget(QLabel("B or K"))
        leftLayout.addWidget(self.glideModeComboBox)


        # Save and Close buttons
        self.saveButton = QPushButton("Start prediction", self)
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

    def chooseCsvFilePath(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose excel File Path", "", "xlsx Files (*.xlsx)")
        if file_name:
            self.csvFilePathEdit.setText(file_name)

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)

            
    def save_data(self):
        # Collect values from UI elements
        csv_file_path = self.csvFilePathEdit.text()
        save_dir = self.saveDirEdit.text()
        Glide_mode = self.glideModeComboBox.currentText()

        # Check if all required fields are filled
        if csv_file_path and save_dir and Glide_mode:
            self.run_script(csv_file_path, save_dir,  Glide_mode)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, csv_file_path, save_dir,  Glide_mode):
        # Create and start the ScriptThread with the collected parameters
        self.thread = ScriptThread(csv_file_path, save_dir,  Glide_mode)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()

    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, csv_file_path):
        self.animationMovie.stop()
        # current_path = os.getcwd()
        hand = os.path.splitext(csv_file_path)[0]
        current_path = os.getcwd()
        finishmap_path = f"{current_path}/images/logo.png"
        pixmap = QPixmap(finishmap_path)  # 结束时显示的图片
        scaled_pixmap = pixmap.scaled(self.animationLabel.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)  # 缩放图片
        self.animationLabel.setPixmap(scaled_pixmap)  # 显示缩放后的图片

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = BKdialogs()
    dialog.show()
    sys.exit(app.exec_())