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

    def __init__(self, clus_num, pIC50threshold , handle_file):
        super().__init__()
        self.clus_num = clus_num
        self.pIC50threshold = pIC50threshold
        self.handle_file = handle_file

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/LigScript/clus.py"


        subprocess.call(["python", script_path, self.clus_num, self.pIC50threshold, self.handle_file])
        print ("python", script_path, self.clus_num, self.pIC50threshold, self.handle_file)
        self.finished.emit(self.clus_num, self.handle_file)  # Emit with parameters

class InputclusDialogs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent, Qt.FramelessWindowHint | Qt.WindowSystemMenuHint)
        self.setWindowOpacity(0.9)
        self.initUI()


    def initUI(self):
        # 主布局
        mainLayout = QHBoxLayout(self)

        # 左侧布局（输入区）
        leftLayout = QVBoxLayout()
        # Name input
        self.nameEdit = QLineEdit(self)
        self.nameEdit.setPlaceholderText("Enter clus_num, recommended 50")
        leftLayout.addWidget(QLabel("clus_num"))
        leftLayout.addWidget(self.nameEdit)


        # Uniprot_ID input
        self.PDBIdEdit = QLineEdit(self)
        self.PDBIdEdit.setPlaceholderText("Enter pIC50threshold Attention 5 is 10uM and 6 is 1uM")
        leftLayout.addWidget(QLabel("pIC50threshold"))
        leftLayout.addWidget(self.PDBIdEdit)

        # Save_Place input
        self.savePlaceEdit = QLineEdit(self)
        self.savePlaceEdit.setPlaceholderText("Choose The Ligand csv file for clustering") 
        self.savePlaceButton = QPushButton("Choose handle_file", self)
        self.savePlaceButton.clicked.connect(self.chooseSavePlace)
        leftLayout.addWidget(QLabel("handle_file"))
        leftLayout.addWidget(self.savePlaceEdit)
        leftLayout.addWidget(self.savePlaceButton)

        # Save and Close buttons
        self.saveButton = QPushButton("Lig_Clus", self)
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

    def chooseSavePlace(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", "", "CSV Files (*.csv)")
        if file_name:
            self.savePlaceEdit.setText(file_name)


    def save_data(self):
        clus_num = self.nameEdit.text()
        pIC50threshold = self.PDBIdEdit.text()
        save_place = self.savePlaceEdit.text()

        if clus_num and pIC50threshold and save_place:
            self.run_script(clus_num, pIC50threshold, save_place)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, clus_num, pIC50threshold,  save_place):
        self.thread = ScriptThread(clus_num, pIC50threshold,  save_place)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, clus_num, handle_file):
        self.animationMovie.stop()
        # current_path = os.getcwd()
        hand = os.path.splitext(handle_file)[0]
        finishmap_path = f"{hand}_tSNE.png"
        pixmap = QPixmap(finishmap_path)  # 结束时显示的图片
        scaled_pixmap = pixmap.scaled(self.animationLabel.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)  # 缩放图片
        self.animationLabel.setPixmap(scaled_pixmap)  # 显示缩放后的图片

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = InputclusDialogs()
    dialog.show()
    sys.exit(app.exec_())
