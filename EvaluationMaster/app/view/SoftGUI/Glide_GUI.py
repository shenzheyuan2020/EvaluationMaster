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

    def __init__(self, csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode,threads_num):
        super().__init__()
        self.csv_file_path = csv_file_path
        self.lig_file = lig_file
        self.schro_dir = schro_dir
        self.mm_sahre_dir = mm_sahre_dir
        self.save_dir = save_dir
        self.Glide_mode = Glide_mode
        self.threads_num = threads_num
    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/DockingScript/GlideG.py"
        subprocess.call([f"python", script_path, self.csv_file_path, self.lig_file, self.schro_dir, self.mm_sahre_dir, self.save_dir, self.Glide_mode, self.threads_num])
        self.finished.emit(self.csv_file_path)  # Emit with parameters

class glidedialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # Main layout
        mainLayout = QHBoxLayout(self)

        # Left layout (Input area)
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

           # Add new UI elements for csv_file_path
        self.csvFilePathEdit = QLineEdit(self)
        self.csvFilePathButton = QPushButton("Choose CSV File Path", self)
        self.csvFilePathButton.clicked.connect(self.chooseCsvFilePath)
        leftLayout.addWidget(QLabel("CSV File Path"))
        leftLayout.addWidget(self.csvFilePathEdit)
        leftLayout.addWidget(self.csvFilePathButton)

        # UI elements for lig_file
        self.ligFilePathEdit = QLineEdit(self)
        self.ligFilePathButton = QPushButton("Choose Ligand File", self)
        self.ligFilePathButton.clicked.connect(self.chooseLigFilePath)
        leftLayout.addWidget(QLabel("Ligand File"))
        leftLayout.addWidget(self.ligFilePathEdit)
        leftLayout.addWidget(self.ligFilePathButton)

        # UI elements for schro_dir
        self.schroDirEdit = QLineEdit(self)
        self.schroDirButton = QPushButton("Choose Schrodir Directory", self)
        self.schroDirButton.clicked.connect(self.chooseSchroDir)
        leftLayout.addWidget(QLabel("Schrodinger Directory"))
        leftLayout.addWidget(self.schroDirEdit)
        leftLayout.addWidget(self.schroDirButton)

        # UI elements for mm_sahre_dir
        self.mmShareDirEdit = QLineEdit(self)
        self.mmShareDirButton = QPushButton("Choose MMshare Directory", self)
        self.mmShareDirButton.clicked.connect(self.chooseMMshareDir)
        leftLayout.addWidget(QLabel("MMshare Directory"))
        leftLayout.addWidget(self.mmShareDirEdit)
        leftLayout.addWidget(self.mmShareDirButton)

        # UI elements for save_dir
        self.saveDirEdit = QLineEdit(self)
        self.saveDirButton = QPushButton("Choose Save Directory", self)
        self.saveDirButton.clicked.connect(self.chooseSaveDir)
        leftLayout.addWidget(QLabel("Save Directory"))
        leftLayout.addWidget(self.saveDirEdit)
        leftLayout.addWidget(self.saveDirButton)

        # UI elements for Glide_mode
        self.glideModeComboBox = QComboBox(self)
        self.glideModeComboBox.addItems(["SP", "XP"])
        leftLayout.addWidget(QLabel("Glide Mode"))
        leftLayout.addWidget(self.glideModeComboBox)


        # UI elements for threads_num
        self.threadsNumEdit = QLineEdit(self)
        self.threadsNumEdit.setPlaceholderText("Enter threads number")
        leftLayout.addWidget(QLabel("Threads Number"))
        leftLayout.addWidget(self.threadsNumEdit)                                               

        # Save and Close buttons
        self.saveButton = QPushButton("Start Docking", self)
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
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File Path", "", "CSV Files (*.csv)")
        if file_name:
            self.csvFilePathEdit.setText(file_name)

    def chooseLigFilePath(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose Ligand File", "", "All Files (*)")
        if file_name:
            self.ligFilePathEdit.setText(file_name)

    def chooseSchroDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Schrodinger Directory")
        if directory:
            self.schroDirEdit.setText(directory)

    def chooseMMshareDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose MMshare Directory")
        if directory:
            self.mmShareDirEdit.setText(directory)

    def chooseSaveDir(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Save Directory")
        if directory:
            self.saveDirEdit.setText(directory)

            
    def save_data(self):
        # Collect values from UI elements
        csv_file_path = self.csvFilePathEdit.text()
        lig_file = self.ligFilePathEdit.text()
        schro_dir = self.schroDirEdit.text()
        mm_sahre_dir = self.mmShareDirEdit.text()
        save_dir = self.saveDirEdit.text()
        Glide_mode = self.glideModeComboBox.currentText()
        threads_num = self.threadsNumEdit.text()
        # Check if all required fields are filled
        if csv_file_path and lig_file and schro_dir and mm_sahre_dir and save_dir:
            self.run_script(csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode, threads_num)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode, threads_num):
        # Create and start the ScriptThread with the collected parameters
        self.thread = ScriptThread(csv_file_path, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode, threads_num)
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
    dialog = glidedialogs()
    dialog.show()
    sys.exit(app.exec_())