from PyQt5 import QtCore, QtGui, QtWidgets
# from qt_material import apply_stylesheet
from qdarkstyle.light.palette import LightPalette
import qdarkstyle
from PyQt5.QtWidgets import QFileDialog,QMessageBox
from PyQt5.QtCore import Qt, QPoint, Qt, QThread, pyqtSignal
import sys
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QTableWidget, QTableWidgetItem, QLineEdit, QGraphicsView, QFileDialog
import os
import subprocess
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QGraphicsPixmapItem
import csv
from PyQt5.QtGui import QMovie


class ScriptThread(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, deep_pocket_path, input_folder, excel_file_path, gpu_num):
        super().__init__()
        self.deep_pocket_path = deep_pocket_path
        self.input_folder = input_folder
        self.excel_file_path = excel_file_path
        self.gpu_num = gpu_num

    def run(self):
        self.started.emit()
        current_path = os.getcwd() 
        script_path = os.path.join(current_path, "app/Script/ProScript/predict_pocket.py")
        subprocess.call([sys.executable, script_path, self.deep_pocket_path, self.input_folder, self.excel_file_path, str(self.gpu_num)])
        self.finished.emit(self.excel_file_path)


class DeeppocketDialog(QDialog):
    def __init__(self, parent=None):
        super(DeeppocketDialog, self).__init__(parent)
        self.setObjectName("Dialog")
        self.resize(844, 603)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)
        self.initUI()

    def initUI(self):
        self.showexcel = QTableWidget(self)
        self.showexcel.setGeometry(QtCore.QRect(40, 20, 761, 181))
        self.showexcel.setColumnCount(15)
        self.showexcel.setRowCount(12)

        self.show_graphicsView = QGraphicsView(self)
        self.show_graphicsView.setGeometry(QtCore.QRect(410, 220, 391, 311))

        self.Deeppocket_dir = QLineEdit(self)
        self.Deeppocket_dir.setGeometry(QtCore.QRect(40, 220, 351, 25))

        self.Line_protein_path = QLineEdit(self)
        self.Line_protein_path.setGeometry(QtCore.QRect(40, 300, 351, 25))

        self.Line_outdir = QLineEdit(self)
        self.Line_outdir.setGeometry(QtCore.QRect(40, 380, 351, 25))

        self.GPUnumber = QLineEdit(self)
        self.GPUnumber.setGeometry(QtCore.QRect(40, 460, 351, 25))

        self.button_close = QPushButton("Close", self)
        self.button_close.setGeometry(QtCore.QRect(350, 550, 80, 25))
        self.button_close.clicked.connect(self.close)

        self.button_run = QPushButton("Run DeepPocket", self)
        self.button_run.setGeometry(QtCore.QRect(40, 500, 351, 25))

        self.button_outdir = QPushButton("Choose", self)
        self.button_outdir.setGeometry(QtCore.QRect(40, 420, 351, 25))

        self.BButton_Deeppocketdir = QPushButton("DeepPocket dir", self)
        self.BButton_Deeppocketdir.setGeometry(QtCore.QRect(40, 260, 351, 25))

        self.Button_protein = QPushButton("Choose", self)
        self.Button_protein.setGeometry(QtCore.QRect(40, 340, 351, 25))

        self.setupConnections()
        self.setupLineEditWidgets()
        self.displayLogo()
        self.button_run.clicked.connect(self.onRunClicked)

    def setupLineEditWidgets(self):
        current_path = os.getcwd()
        evaluation_master = os.getenv('EVALUATIONMASTER', '')  # Get environment variable
        self.Deeppocket_dir.setText(os.path.join(evaluation_master, "Support_software/DeepPocket/"))
        self.Line_protein_path.setPlaceholderText("PATH_TO_WHERE_YOU_SAVE_PDB_FILE")
        self.Line_outdir.setPlaceholderText("PATH_TO_WHERE_YOU_SAVE_THE_OUTPUT")
        self.GPUnumber.setText("1")

    def setupConnections(self):
        self.button_outdir.clicked.connect(self.chooseoutDir)
        self.BButton_Deeppocketdir.clicked.connect(self.choosedeeppocketDir)
        self.Button_protein.clicked.connect(self.chooseproteinDir)

    def choosedeeppocketDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.Deeppocket_dir.setText(dir)  

    def chooseoutDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.Line_outdir.setText(dir)  

    def chooseproteinDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.Line_protein_path.setText(dir)  

    def displayLogo(self):
        pixmap = QPixmap("images/logo.png")
        if not pixmap.isNull():
            scaledPixmap = pixmap.scaled(self.show_graphicsView.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            item = QGraphicsPixmapItem(scaledPixmap)
            scene = QtWidgets.QGraphicsScene(self)
            scene.addItem(item)
            self.show_graphicsView.setScene(scene)

    def resizeEvent(self, event):
        super(DeeppocketDialog, self).resizeEvent(event)
        self.displayLogo()

    def onRunClicked(self):
        Deeppocket_dir = self.Deeppocket_dir.text().strip()
        protein_path = self.Line_protein_path.text().strip()
        outdir = self.Line_outdir.text().strip()
        gpu_num = self.GPUnumber.text().strip()

        if not Deeppocket_dir or not protein_path or not outdir or not gpu_num:
            QMessageBox.information(self, "Script Finished", "The DeepPocket script has completed.")
            return

        self.displayLoadingGif()
        self.thread = ScriptThread(Deeppocket_dir, protein_path, outdir, gpu_num)
        self.thread.started.connect(self.onScriptStarted)
        self.thread.finished.connect(self.onScriptFinished)
        self.thread.start()

    def onScriptStarted(self):
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished(self, out_dir):
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()
        QMessageBox.information(self, "Script Finished", f"The DeepPocket script has completed. Output directory: {out_dir}")

        csv_file_path = os.path.join(out_dir, "coordinates.csv")
        self.displayCSV(csv_file_path)

    def displayCSV(self, csv_file_path):
        with open(csv_file_path, 'r') as file:
            reader = csv.reader(file)
            data = list(reader)
            if data:
                self.showexcel.setRowCount(len(data))
                self.showexcel.setColumnCount(len(data[0]))
                for row_index, row in enumerate(data):
                    for column_index, item in enumerate(row):
                        self.showexcel.setItem(row_index, column_index, QTableWidgetItem(item))

    def displayLoadingGif(self):
        if hasattr(self, 'loadingLabel') and self.loadingLabel:
            self.loadingLabel.deleteLater()
        
        self.loadingLabel = QtWidgets.QLabel(self.show_graphicsView)
        self.loadingLabel.setAlignment(Qt.AlignCenter)
        self.loadingLabel.setGeometry(0, 0, self.show_graphicsView.width(), self.show_graphicsView.height())
        current_path = os.getcwd()
        movie_dir = os.path.join(current_path, "images/loading.gif")
        movie = QMovie(movie_dir)
        self.loadingLabel.setMovie(movie)
        movie.setScaledSize(self.loadingLabel.size())
        movie.start()

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

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    Dialog = DeeppocketDialog()
    stylesheet_path = "style.qss"
    Dialog.setWindowFlags(Qt.FramelessWindowHint)
    DeeppocketDialog.apply_stylesheet(app, stylesheet_path)
    Dialog.show()
    sys.exit(app.exec_())
