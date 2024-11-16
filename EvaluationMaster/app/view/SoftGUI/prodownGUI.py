from PyQt5 import QtCore, QtGui, QtWidgets
# from qt_material import apply_stylesheet
from qdarkstyle.light.palette import LightPalette
import qdarkstyle
from PyQt5.QtWidgets import QFileDialog,QMessageBox
from PyQt5.QtCore import Qt, QPoint, Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import QComboBox, QFileDialog, QMessageBox
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

    def __init__(self, uniprot_id, name, resolution_threshold, method_filter, save_path):
        super().__init__()
        self.uniprot_id = uniprot_id
        self.name = name
        self.resolution_threshold = resolution_threshold
        self.method_filter = method_filter
        self.save_path = save_path

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/ProScript/downpdb.py"

        subprocess.call(["python", script_path, self.uniprot_id, self.name, self.resolution_threshold, self.method_filter, self.save_path])
        self.finished.emit(self.name, self.save_path)  # Emit with parameters

class InputprodownDialogs(QDialog):
    def __init__(self, parent=None):
        super(InputprodownDialogs, self).__init__(parent)
        self.setObjectName("Dialog")
        self.resize(765, 502)
        
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)  ###windows property
        # self.setWindowFlags(QtCore.Qt.FramelessWindowHint)  # Removes the title bar
        self.initUI()

    def initUI(self):
        # Table Widget
        self.tableWidget = QTableWidget(self)
        self.tableWidget.setGeometry(QtCore.QRect(20, 20, 721, 151))
        self.tableWidget.setColumnCount(15)
        self.tableWidget.setRowCount(12)

        # Graphics View
        self.graphicsView = QGraphicsView(self)
        self.graphicsView.setGeometry(QtCore.QRect(350, 180, 391, 271))

        # Line Edits
        self.UniprotID = QLineEdit(self)
        self.UniprotID.setGeometry(QtCore.QRect(20, 180, 321, 25))

        self.NameforDown = QLineEdit(self)
        self.NameforDown.setGeometry(QtCore.QRect(20, 230, 321, 25))

        self.SaveDirdown = QLineEdit(self)
        self.SaveDirdown.setGeometry(QtCore.QRect(20, 380, 221, 25))

        self.Resolution_thread = QLineEdit(self)
        self.Resolution_thread.setGeometry(QtCore.QRect(20, 280, 321, 25))

        # Buttons
        self.ChooseDir1 = QPushButton("Choose", self)
        self.ChooseDir1.setGeometry(QtCore.QRect(240, 380, 101, 25))

        self.Download = QPushButton("Download", self)
        self.Download.setGeometry(QtCore.QRect(20, 420, 321, 31))

        self.closeButton = QPushButton("Close", self)
        self.closeButton.setGeometry(QtCore.QRect(260, 460, 191, 25))
        self.closeButton.clicked.connect(self.close)

        # QComboBox 
        self.method = QComboBox(self)
        self.method.setGeometry(QtCore.QRect(20, 330, 321, 25))
        self.method.setObjectName("method")

        # Adding items to method QComboBox
        self.method.addItem("x-ray")  # First item
        self.method.addItem("nmr")    # Second item
        self.method.addItem("alphaFold")    # Third item

        # Setup Connections
        self.setupConnections()
        self.setupLineEditWidgets()
        self.displayLogo()
        self.Download.clicked.connect(self.onDownloadClicked)

#########################

    def setupLineEditWidgets(self):
            self.UniprotID.setPlaceholderText("Enter UniProt ID")
            self.NameforDown.setPlaceholderText("Enter Name for Download")
            self.SaveDirdown.setPlaceholderText("Directory Path")
            self.Resolution_thread.setPlaceholderText("Enter Resolution Threshold Value")


    def setupConnections(self):
        self.ChooseDir1.clicked.connect(self.chooseSaveDir)

    def chooseSaveDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.SaveDirdown.setText(dir)

    def displayLogo(self):
        # Load the image
        pixmap = QPixmap("images/logo.png")
        
        # Check if the image was successfully loaded
        if not pixmap.isNull():
            # Scale the QPixmap to fit the size of the graphicsView, maintaining the aspect ratio
            scaledPixmap = pixmap.scaled(self.graphicsView.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            
            # Create a QGraphicsPixmapItem with the scaled QPixmap
            item = QGraphicsPixmapItem(scaledPixmap)
            
            # Create a new QGraphicsScene for displaying the item
            scene = QtWidgets.QGraphicsScene(self)
            scene.addItem(item)
            
            # Set the QGraphicsScene to the graphicsView
            self.graphicsView.setScene(scene)
        else:
            print("Failed to load the image.")

    def resizeEvent(self, event):
        super(InputprodownDialogs, self).resizeEvent(event)
        # Call displayLogo again to resize the image with the new dimensions
        self.displayLogo()

    def onDownloadClicked(self):
        # 从界面读取数据
        uniprot_id = self.UniprotID.text().strip()
        name = self.NameforDown.text().strip()
        save_dir = self.SaveDirdown.text().strip()
        resolution_threshold = self.Resolution_thread.text().strip()

        # 获取 IC50kdki 的值
        method_filter = self.method.currentText()  # 获取选中的文本

        # 验证输入
        if not uniprot_id or not name or not save_dir:
            QMessageBox.warning(self, "Input Error", "Please fill all the fields.")
            return

        # 显示加载GIF
        self.displayLoadingGif()

        self.thread = ScriptThread(uniprot_id, name, resolution_threshold, method_filter, save_dir)
        self.thread.started.connect(self.onScriptStarted)
        self.thread.finished.connect(self.onScriptFinished)  # 完成时连接到 onScriptFinished 槽
        self.thread.start()

    def onScriptStarted(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished(self, name, save_dir):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()
        QMessageBox.information(self, "Download Complete", f"The data for {name} has been saved to {save_dir}.")

        # 使用 selected_type 动态生成文件名
        csv_file_path = os.path.join(save_dir, f"{name}_pdb_info.csv")
        self.displayCSV(csv_file_path)

    def displayCSV(self, csv_file_path):
        with open(csv_file_path, 'r') as file:
            reader = csv.reader(file)
            data = list(reader)
            if data:
                self.tableWidget.setRowCount(len(data))
                self.tableWidget.setColumnCount(len(data[0]))
                for row_index, row in enumerate(data):
                    for column_index, item in enumerate(row):
                        self.tableWidget.setItem(row_index, column_index, QTableWidgetItem(item))

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
        
        # 在QGraphicsView中显示加载GIF
        self.loadingLabel = QtWidgets.QLabel(self.graphicsView)
        self.loadingLabel.setAlignment(Qt.AlignCenter)
        self.loadingLabel.setGeometry(0, 0, self.graphicsView.width(), self.graphicsView.height())
        current_path = os.getcwd()
        movie_dir = current_path + "/images/loading.gif"
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
            
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = InputprodownDialogs()
    stylesheet_path = "style.qss"  # 样式表文件的路径
    # Set the window to be frameless
    Dialog.setWindowFlags(Qt.FramelessWindowHint)
    InputprodownDialogs.apply_stylesheet(app, stylesheet_path)
    #apply_stylesheet(app, theme='light_teal.xml')
    # setup stylesheet
    # app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
 
    # or in new API
    # app.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyqt5', palette=LightPalette))

    Dialog.show()
    sys.exit(app.exec_())

