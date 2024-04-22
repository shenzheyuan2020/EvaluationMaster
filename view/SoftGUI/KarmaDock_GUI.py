import sys
from PyQt5 import QtCore, QtGui, QtWidgets
# from qt_material import apply_stylesheet
from qdarkstyle.light.palette import LightPalette
import qdarkstyle
from PyQt5.QtWidgets import QFileDialog,QMessageBox
from PyQt5.QtCore import Qt, QPoint, Qt, QThread, pyqtSignal

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
    finished = pyqtSignal(str, str)

    def __init__(self, KarmaDock_dir, multi_protein_csv, multi_protein_dir, lig_csv, out_dir, score_threshold):
        super().__init__()
        self.KarmaDock_dir = KarmaDock_dir
        self.multi_protein_csv = multi_protein_csv
        self.multi_protein_dir = multi_protein_dir
        self.lig_csv = lig_csv
        self.out_dir = out_dir
        self.score_threshold = score_threshold

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = os.path.join(current_path, "app/Script/DockingScript/Karmadock.py")
        process = subprocess.Popen(
            [f"{self.KarmaDock_dir}/Env/karmadock_env/bin/python", "-u", script_path, self.KarmaDock_dir, self.multi_protein_csv, self.multi_protein_dir, self.lig_csv, self.out_dir, str(self.score_threshold)],
            stdout=subprocess.PIPE,
            text=True)
        output, _ = process.communicate()
        print("DEBUG OUTPUT:")
        print(output)
        # 初始化变量，以防找不到RESULT行
        mol_name = "Not found"
        out_dir = "Not found"
        # 查找以"RESULT,"开头的行
        for line in output.strip().split("\n"):
            if line.startswith("RESULT,"):
                _, mol_name, out_dir = line.split(",")
                break  # 找到RESULT行后退出循环
        
        self.finished.emit(mol_name, out_dir)



#######################################
class KarmaDock_Dialog(QDialog):
    def __init__(self, parent=None):
        super(KarmaDock_Dialog, self).__init__(parent)
        self.setObjectName("Dialog")
        self.resize(840, 590)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)
        #self.setWindowFlags(QtCore.Qt.FramelessWindowHint)  # Removes the title bar
        self.initUI()

    def initUI(self):
        # Table Widget
        self.showexcel = QTableWidget(self)
        self.showexcel.setGeometry(QtCore.QRect(40, 20, 761, 181))
        self.showexcel.setColumnCount(15)
        self.showexcel.setRowCount(12)
  
        # Graphics View
        self.show_graphicsView =  QGraphicsView(self)
        self.show_graphicsView.setGeometry(QtCore.QRect(410, 220, 391, 311))

        # Line Edits
        self.Line_pro_coor_file = QLineEdit(self)
        self.Line_pro_coor_file.setGeometry(QtCore.QRect(40, 300, 251, 25))

        self.Line_ligand_path = QLineEdit(self)
        self.Line_ligand_path.setGeometry(QtCore.QRect(40, 380, 251, 25))

        self.Line_outdir = QLineEdit(self)
        self.Line_outdir.setGeometry(QtCore.QRect(40, 420, 251, 25))

        self.Line_scorethread = QLineEdit(self)
        self.Line_scorethread.setGeometry(QtCore.QRect(40, 460, 351, 25))

        self.Line_protein_path = QLineEdit(self)
        self.Line_protein_path.setGeometry(QtCore.QRect(40, 340, 251, 25))

        self.karmadock_dir = QLineEdit(self)
        self.karmadock_dir.setGeometry(QtCore.QRect(40, 220, 351, 25))

        # Buttons
        self.Button_pro_coor_file = QPushButton("Choose", self)
        self.Button_pro_coor_file.setGeometry(QtCore.QRect(290, 300, 101, 25))
   
        self.button_close = QPushButton("Close", self)
        self.button_close.setGeometry(QtCore.QRect(350, 550, 80, 25))
        self.button_close.clicked.connect(self.close)

        self.Button_ligand = QPushButton("Choose", self)
        self.Button_ligand.setGeometry(QtCore.QRect(290, 380, 101, 25))
    
        self.button_outdir = QPushButton("Choose", self)
        self.button_outdir.setGeometry(QtCore.QRect(290, 420, 101, 25))

        self.button_run = QPushButton("Run KarmaDock", self)
        self.button_run.setGeometry(QtCore.QRect(40, 500, 351, 25))

        self.Button_karmadockdir = QPushButton("KarmaDock dir", self)
        self.Button_karmadockdir.setGeometry(QtCore.QRect(40, 260, 351, 25))

        self.Button_protein = QPushButton("Choose", self)
        self.Button_protein.setGeometry(QtCore.QRect(290, 340, 101, 25))

        # Setup Connections
        self.setupConnections()
        self.setupLineEditWidgets()
        self.displayLogo()
        self.button_run.clicked.connect(self.onRunClicked)

###############################
    def setupLineEditWidgets(self):
        self.Line_pro_coor_file.setPlaceholderText("Coordinate file (.csv) path")
        self.Line_ligand_path.setPlaceholderText("Choose ligand path (.csv)")
        self.Line_outdir.setPlaceholderText("Choose out dir")
        self.Line_scorethread.setPlaceholderText("Set Score thread")
        self.karmadock_dir.setPlaceholderText( "KarmaDock dir")
        self.Line_protein_path.setPlaceholderText("Choose protein(.pdb) path")

    def setupConnections(self):
        self.button_outdir.clicked.connect(self.chooseoutDir)
        self.Button_karmadockdir.clicked.connect(self.choosekarmaDir)
        self.Button_protein.clicked.connect(self.chooseproteinDir)
        self.Button_pro_coor_file.clicked.connect(lambda: self.chooseFile(self.Line_pro_coor_file, "CSV Files (*.csv)"))
        self.Button_ligand.clicked.connect(lambda: self.chooseFile(self.Line_ligand_path, "CSV Files (*.csv)"))

    def choosekarmaDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.karmadock_dir.setText(dir)  

    def chooseproteinDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.Line_protein_path.setText(dir)  

    def chooseoutDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.Line_outdir.setText(dir)  

    def chooseFile(self, lineEdit, fileType):
        file, _ = QFileDialog.getOpenFileName(self, "Select File", "", fileType)
        if file:
            lineEdit.setText(file)

    def displayLogo(self):
        # Load the image
        pixmap = QPixmap("images/logo.png")
        
        # Check if the image was successfully loaded
        if not pixmap.isNull():
            # Scale the QPixmap to fit the size of the show_graphicsView, maintaining the aspect ratio
            scaledPixmap = pixmap.scaled(self.show_graphicsView.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            
            # Create a QGraphicsPixmapItem with the scaled QPixmap
            item = QGraphicsPixmapItem(scaledPixmap)
            
            # Create a new QGraphicsScene for displaying the item
            scene = QtWidgets.QGraphicsScene(self)
            scene.addItem(item)
            
            # Set the QGraphicsScene to the show_graphicsView
            self.show_graphicsView.setScene(scene)
        else:
            print("Failed to load the image.")

    def resizeEvent(self, event):
        super(KarmaDock_Dialog, self).resizeEvent(event)
        self.displayLogo()

    def onRunClicked(self):###################
        coor_file = self.Line_pro_coor_file.text().strip()
        lig_file = self.Line_ligand_path.text().strip()
        outdir = self.Line_outdir.text().strip()
        scorethread = self.Line_scorethread.text().strip()
        protein_path = self.Line_protein_path.text().strip()
        karmadock_dir = self.karmadock_dir.text().strip()

        # 验证输入
        if not karmadock_dir or not coor_file or not protein_path or not lig_file or not outdir or not scorethread:
            QMessageBox.information(self, "Script Finished", "The KarmaDock script has completed.")
            return

        # 显示加载GIF
        self.displayLoadingGif()

        # 初始化并启动下载线程
        self.thread = ScriptThread(karmadock_dir, coor_file, protein_path, lig_file, outdir, scorethread)
        self.thread.started.connect(self.onScriptStarted)
        self.thread.finished.connect(self.onScriptFinished)  # 完成时连接到onScriptFinished槽
        self.thread.start()

    def onScriptStarted(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished(self, mol_name, out_dir):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()

        # 显示完成信息
        QMessageBox.information(self, "Script Finished", f"The KarmaDock script has completed. Output directory: {out_dir}")

        # 使用mol_name变量构建CSV文件路径
        csv_file_path = os.path.join(out_dir, f"{mol_name}_score.csv")
        self.displayCSV(csv_file_path)


#####################################
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

    def displayImage(self, image_file_path):
        # 加载图片
        pixmap = QPixmap(image_file_path)

        if not pixmap.isNull():
            # 获取QGraphicsView的尺寸
            rect = self.graphicsView.rect()
            
            # 缩放图片以适应QGraphicsView的尺寸，同时保持纵横比
            scaledPixmap = pixmap.scaled(rect.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            
            # 创建QGraphicsPixmapItem并使用缩放后的图片
            item = QGraphicsPixmapItem(scaledPixmap)
            
            # 创建QGraphicsScene并添加item
            scene = QtWidgets.QGraphicsScene(self)
            scene.addItem(item)
            
            # 将创建的scene设置给QGraphicsView
            self.graphicsView.setScene(scene)
        else:
            QMessageBox.warning(self, "Display Error", "Failed to load the image.")

    def displayLoadingGif(self):
        # 如果已有，则先清除现有的场景
        if hasattr(self, 'loadingLabel') and self.loadingLabel:
            self.loadingLabel.deleteLater()
        
        # 在show_graphicsView中显示加载GIF
        self.loadingLabel = QtWidgets.QLabel(self.show_graphicsView)
        self.loadingLabel.setAlignment(Qt.AlignCenter)
        # 注意这里的宽度和高度设置可能需要根据实际情况进行调整
        self.loadingLabel.setGeometry(0, 0, self.show_graphicsView.width(), self.show_graphicsView.height())
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
        movie = QMovie(movie_dir)
        self.loadingLabel.setMovie(movie)
        # 根据实际需要调整尺寸
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
    Dialog = KarmaDock_Dialog()
    stylesheet_path = "style.qss"  # 样式表文件的路径
    Dialog.setWindowFlags(Qt.FramelessWindowHint)
    Dialog.apply_stylesheet(app, stylesheet_path)
    Dialog.show()
    # sys.exit(app.exec_())

