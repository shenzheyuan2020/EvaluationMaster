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
    finished = pyqtSignal(str, str, str)

    def __init__(self, uniprot_id, name, save_dir, selected_type):
        super().__init__()
        self.uniprot_id = uniprot_id
        self.name = name
        self.save_dir = save_dir
        self.selected_type = selected_type

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/LigScript/down-step1.py"

        subprocess.call(["python", script_path, self.uniprot_id, self.name, self.save_dir, self.selected_type])
        self.finished.emit(self.name, self.save_dir, self.selected_type)  # Emit with parameters

### script Thread
class ScriptThread_filter(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str, str, str)

    def __init__(self, csv_file_path, filter_name, save_dir, selected_type_f, threshold):
        super().__init__()
        self.csv_file_path = csv_file_path
        self.filter_name = filter_name
        self.save_dir = save_dir
        self.selected_type_f = selected_type_f
        self.threshold = threshold

    def run(self):  
        self.started.emit()
        current_path = os.getcwd()
        filter_script_path = current_path + "/app/Script/LigScript/filter-step2.py"

        subprocess.call(["python", filter_script_path, self.csv_file_path, self.filter_name,  self.save_dir, self.selected_type_f, self.threshold])
        self.finished.emit(self.filter_name, self.save_dir, self.selected_type_f)

### clean script Thread#############
class ScriptThread_clean(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str, str)

    def __init__(self, csv_file_path, filter_name, save_dir):
        super().__init__()
        self.csv_file_path = csv_file_path
        self.filter_name = filter_name
        self.save_dir = save_dir

    def run(self):  
        self.started.emit()
        current_path = os.getcwd()
        clean_script_path = current_path + "/app/Script/LigScript/clean-step3.py"

        subprocess.call(["python", clean_script_path, self.csv_file_path, self.filter_name, self.save_dir])
        self.finished.emit(self.filter_name, self.save_dir)
        

#######################################        

class LigprepDialog(QDialog):
    def __init__(self, parent=None):
        super(LigprepDialog, self).__init__(parent)
        self.setObjectName("Dialog")
        self.resize(1034, 590)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)  ###windows property
        # self.setWindowFlags(QtCore.Qt.FramelessWindowHint)  # Removes the title bar
        self.initUI()


    def initUI(self):
        # Table Widget
        self.tableWidget = QTableWidget(self)
        self.tableWidget.setGeometry(QtCore.QRect(20, 50, 481, 151))
        self.tableWidget.setColumnCount(15)
        self.tableWidget.setRowCount(12)

        # Graphics View
        self.graphicsView = QGraphicsView(self)
        self.graphicsView.setGeometry(QtCore.QRect(510, 50, 501, 471))

        # Line Edits
        self.UniprotID = QLineEdit(self)
        self.UniprotID.setGeometry(QtCore.QRect(20, 220, 181, 25))
        
        self.NameforDown = QLineEdit(self)
        self.NameforDown.setGeometry(QtCore.QRect(210, 220, 151, 25))
        
        self.SaveDirdown = QLineEdit(self)
        self.SaveDirdown.setGeometry(QtCore.QRect(20, 260, 241, 25))
        
        self.LigforPre = QLineEdit(self)
        self.LigforPre.setGeometry(QtCore.QRect(20, 340, 241, 25))
        
        self.Threshold = QLineEdit(self)
        self.Threshold.setGeometry(QtCore.QRect(20, 380, 241, 25))
        
        self.Name = QLineEdit(self)
        self.Name.setGeometry(QtCore.QRect(270, 380, 91, 25))
        
        self.LigforClean = QLineEdit(self)
        self.LigforClean.setGeometry(QtCore.QRect(20, 460, 241, 25))

        # Buttons
        self.ChooseDir1 = QPushButton("Choose Directory 1", self)
        self.ChooseDir1.setGeometry(QtCore.QRect(270, 260, 231, 25))
        
        self.Download = QPushButton("Download", self)
        self.Download.setGeometry(QtCore.QRect(20, 300, 481, 25))
        
        self.ChooseDir2 = QPushButton("Choose Directory 2", self)
        self.ChooseDir2.setGeometry(QtCore.QRect(270, 340, 231, 25))
        
        self.filter = QPushButton("Filter", self)
        self.filter.setGeometry(QtCore.QRect(20, 420, 481, 25))
        
        self.ChooseDir3 = QPushButton("Choose Directory 3", self)
        self.ChooseDir3.setGeometry(QtCore.QRect(270, 460, 231, 25))
        
        self.CleanButton = QPushButton("Clean", self)
        self.CleanButton.setGeometry(QtCore.QRect(20, 500, 481, 25))
        
        self.closeButton = QPushButton("Close", self)
        self.closeButton.setGeometry(QtCore.QRect(250, 550, 481, 25))
        self.closeButton.clicked.connect(self.close)

        # QComboBox for IC50kdki
        self.IC50kdki = QComboBox(self)
        self.IC50kdki.setGeometry(QtCore.QRect(371, 220, 131, 25))
        self.IC50kdki.setObjectName("IC50kdki")

        # Adding items to IC50kdki QComboBox
        self.IC50kdki.addItem("IC50")  # First item
        self.IC50kdki.addItem("Kd")    # Second item
        self.IC50kdki.addItem("Ki")    # Third item

        # QComboBox for IC50kdki
        self.IC50kdki_2 = QComboBox(self)
        self.IC50kdki_2.setGeometry(QtCore.QRect(370, 380, 131, 25))
        self.IC50kdki_2.setObjectName("IC50kdki_2")

        # Adding items to IC50kdki QComboBox
        self.IC50kdki_2.addItem("pIC50")  # First item
        self.IC50kdki_2.addItem("pKd")    # Second item
        self.IC50kdki_2.addItem("pKi")    # Third item

        # Setup Connections
        self.setupConnections()
        self.setupLineEditWidgets()
        self.displayLogo()
        self.Download.clicked.connect(self.onDownloadClicked)

###  Filter function
        self.filter.clicked.connect(self.filter_Clicked)
###  clean function
        self.CleanButton.clicked.connect(self.clean_Clicked)

###############
    def setupLineEditWidgets(self):
            self.UniprotID.setPlaceholderText("Enter UniProt ID")
            self.NameforDown.setPlaceholderText("Enter Name")
            self.SaveDirdown.setPlaceholderText("Output Directory")
            self.LigforPre.setPlaceholderText("Ligand CSV File Path for Filtering")
            self.Threshold.setPlaceholderText("Enter Threshold Value")
            self.Name.setPlaceholderText("Enter Name")
            self.LigforClean.setPlaceholderText("Ligand CSV File Path for Cleaning")

    def setupConnections(self):
        self.ChooseDir1.clicked.connect(self.chooseSaveDir)
        # Pass the file type filter as an argument using lambda
        self.ChooseDir2.clicked.connect(lambda: self.chooseFile(self.LigforPre, "CSV Files (*.csv)"))
        self.ChooseDir3.clicked.connect(lambda: self.chooseFile(self.LigforClean, "CSV Files (*.csv)"))
        # Repeat for other buttons as needed

    def chooseSaveDir(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.SaveDirdown.setText(dir)

    def chooseFile(self, lineEdit, fileType):
        file, _ = QFileDialog.getOpenFileName(self, "Select File", "", fileType)
        if file:
            lineEdit.setText(file)

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
        super(LigprepDialog, self).resizeEvent(event)
        # Call displayLogo again to resize the image with the new dimensions
        self.displayLogo()

    def onDownloadClicked(self):
        # 从界面读取数据
        uniprot_id = self.UniprotID.text().strip()
        name = self.NameforDown.text().strip()
        save_dir = self.SaveDirdown.text().strip()

        # 获取 IC50kdki 的值
        selected_type = self.IC50kdki.currentText()  # 获取选中的文本

        # 验证输入
        if not uniprot_id or not name or not save_dir:
            QMessageBox.warning(self, "Input Error", "Please fill all the fields.")
            return

        # 显示加载GIF
        self.displayLoadingGif()

        # 初始化并启动下载线程，传入 selected_type 参数
        self.thread = ScriptThread(uniprot_id, name, save_dir, selected_type)
        self.thread.started.connect(self.onScriptStarted)
        self.thread.finished.connect(self.onScriptFinished)  # 完成时连接到 onScriptFinished 槽
        self.thread.start()


    def onScriptStarted(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished(self, name, save_dir, selected_type):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()
        QMessageBox.information(self, "Download Complete", f"The data for {name} has been saved to {save_dir}.")

        # 使用 selected_type 动态生成文件名
        csv_file_path = os.path.join(save_dir, f"{name}_{selected_type}.csv")
        self.displayCSV(csv_file_path)

        # 读取图片文件
        image_file_path = os.path.join(save_dir, f"{name}_p{selected_type}_downlig.png")
        self.displayImage(image_file_path)

###  Filter function
    def filter_Clicked(self):
            # 从界面读取数据
            csv_file_path = self.LigforPre.text().strip()
            filter_name = self.Name.text().strip()
            threshold = self.Threshold.text().strip()
            filter_save_dir = self.SaveDirdown.text().strip()

            # Debug prints
            print(f"'{csv_file_path}', '{filter_name}', '{threshold}'")
            
            selected_type_f = self.IC50kdki_2.currentText()

            # 验证输入
            if not csv_file_path or not filter_name or not threshold:
                QMessageBox.warning(self, "Input Error", "Please fill all the fields.")
                return

            # 显示加载GIF
            self.displayLoadingGif()

            # 初始化并启动
            self.filter_thread = ScriptThread_filter(csv_file_path, filter_name, filter_save_dir, selected_type_f, threshold)
            self.filter_thread.started.connect(self.onScriptStarted_filter)
            self.filter_thread.finished.connect(self.onScriptFinished_filter)
            self.filter_thread.start()

###  A new onScriptFinished function (remember to rename a new one )
    def onScriptStarted_filter(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished_filter(self, filter_name, filter_save_dir, selected_type_f):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()
        QMessageBox.information(self, "Download Complete", f"The data for {filter_name} has been saved to {filter_save_dir}.")

        # 显示CSV文件内容
        csv_file_path = os.path.join(filter_save_dir, f"{filter_name}_{selected_type_f}_filter.csv")
        self.displayCSV(csv_file_path)
        
        # 显示图片
        image_file_path = os.path.join(filter_save_dir, f"{filter_name}_{selected_type_f}_filter.png")
        self.displayImage(image_file_path)

################clean function#####################################
    def clean_Clicked(self):
        # 从界面读取数据
        csv_file_path = self.LigforClean.text().strip()
        filter_name = self.Name.text().strip()
        clean_save_dir = self.SaveDirdown.text().strip()

        # Debug prints
        print(f"'{csv_file_path}', '{filter_name}'")

        # 验证输入
        if not csv_file_path or not filter_name:
            QMessageBox.warning(self, "Input Error", "Please fill all the fields.")
            return

        # 显示加载GIF
        self.displayLoadingGif()

        # 初始化并启动
        self.clean_thread = ScriptThread_clean(csv_file_path, filter_name, clean_save_dir)
        self.clean_thread.started.connect(self.onScriptStarted_clean)
        self.clean_thread.finished.connect(self.onScriptFinished_clean)
        self.clean_thread.start()


    
    def onScriptStarted_clean(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished_clean(self, filter_name, clean_save_dir):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()
        QMessageBox.information(self, "Download Complete", f"The data for {filter_name} has been saved to {clean_save_dir}.")

        # 显示CSV文件内容
        csv_file_path = os.path.join(clean_save_dir, f"{filter_name}_clean.csv")
        self.displayCSV(csv_file_path)
        
        # 显示图片
        image_file_path = os.path.join(clean_save_dir, f"{filter_name}_clean.png")
        self.displayImage(image_file_path)


#############################        
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
    Dialog = LigprepDialog()
    stylesheet_path = "style.qss"  # 样式表文件的路径
    # Set the window to be frameless
    Dialog.setWindowFlags(Qt.FramelessWindowHint)
    LigprepDialog.apply_stylesheet(app, stylesheet_path)
    #apply_stylesheet(app, theme='light_teal.xml')
    # setup stylesheet
    # app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
 
    # or in new API
    # app.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyqt5', palette=LightPalette))

    Dialog.show()
    sys.exit(app.exec_())

