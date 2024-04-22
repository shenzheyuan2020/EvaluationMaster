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

    def __init__(self, pdb_id, name, save_place):
        super().__init__()
        self.pdb_id = pdb_id
        self.name = name
        self.save_place = save_place

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/ProScript/downpdb.py"

        subprocess.call(["python", script_path, self.pdb_id, self.name, self.save_place])
        self.finished.emit(self.name, self.save_place)  # Emit with parameters

class InputprodownDialogs(QDialog):
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
        self.nameEdit.setPlaceholderText("Enter Name")
        leftLayout.addWidget(QLabel("Name"))
        leftLayout.addWidget(self.nameEdit)

        # Uniprot_ID input
        self.PDBIdEdit = QLineEdit(self)
        self.PDBIdEdit.setPlaceholderText("Enter PDB_ID")
        leftLayout.addWidget(QLabel("PDB_ID"))
        leftLayout.addWidget(self.PDBIdEdit)

        # Save_Place input
        self.savePlaceEdit = QLineEdit(self)
        self.savePlaceEdit.setReadOnly(True)
        self.savePlaceButton = QPushButton("Choose Save Place", self)
        self.savePlaceButton.clicked.connect(self.chooseSavePlace)
        leftLayout.addWidget(QLabel("Save_Place"))
        leftLayout.addWidget(self.savePlaceEdit)
        leftLayout.addWidget(self.savePlaceButton)

        # Save and Close buttons
        self.saveButton = QPushButton("Protein_download", self)
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

        # # 设置样式
        # self.setStyleSheet("""
        #     QDialog {
        #         background-color: rgba(255, 255, 255, 200);
        #     }
        #     QLabel, QLineEdit, QPushButton {
        #         /* 样式设置 */
        #     }
        # """)

    def chooseSavePlace(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.savePlaceEdit.setText(directory)

    def save_data(self):
        name = self.nameEdit.text()
        uniprot_id = self.PDBIdEdit.text()
        save_place = self.savePlaceEdit.text()

        if name and uniprot_id and save_place:
            self.run_script(uniprot_id, name, save_place)
        else:
            QMessageBox.warning(self, "Warning", "Please fill all fields")

    def run_script(self, uniprot_id, name, save_place):
        self.thread = ScriptThread(uniprot_id, name, save_place)
        self.thread.started.connect(self.on_script_start)
        self.thread.finished.connect(self.on_script_finish)  # Connect to the slot
        self.thread.start()


    def on_script_start(self):
        self.animationMovie.start()

    def on_script_finish(self, name, save_place):
        self.animationMovie.stop()
        current_path = os.getcwd()
        finishmap_path = f"{current_path}/images/logo.png"
        pixmap = QPixmap(finishmap_path)  # 结束时显示的图片
        scaled_pixmap = pixmap.scaled(self.animationLabel.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)  # 缩放图片
        self.animationLabel.setPixmap(scaled_pixmap)  # 显示缩放后的图片
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
    InputprodownDialogs.apply_stylesheet(app, stylesheet_path)
    dialog = InputprodownDialogs()
    dialog.show()
    sys.exit(app.exec_())
