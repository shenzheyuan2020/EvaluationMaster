from PyQt5 import QtCore, QtGui, QtWidgets
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


class ScriptThread_sum(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, decoy_dir, inh_dir, out_dir):
        super().__init__()
        self.decoy_dir = decoy_dir
        self.inh_dir = inh_dir
        self.out_dir = out_dir

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = os.path.join(current_path, "app/Script/Analysis/anal_sum.py")

        subprocess.call(["python", script_path, self.decoy_dir, self.inh_dir, self.out_dir])
        self.finished.emit(self.out_dir)


###script2 plot
class ScriptThread_plot(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, active_path, decoy_path, output_path):
        super().__init__()
        self.active_path = active_path
        self.decoy_path = decoy_path
        self. output_path = output_path

    def run(self):
        self.started.emit()
        current_path = os.getcwd()
        script_path = os.path.join(current_path, "app/Script/Analysis/anal_plot.py")

        subprocess.call(["python", script_path, self.active_path, self.decoy_path, self.output_path])
        self.finished.emit(self.output_path)

########################################
class AnalDialog(QDialog):
    def __init__(self, parent=None):
        super(AnalDialog, self).__init__(parent)
        self.setObjectName("Dialog")
        self.resize(801, 572)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)  ###windows property
        # self.setWindowFlags(QtCore.Qt.FramelessWindowHint)  # Removes the title bar
        self.initUI()

    def initUI(self):
        # Table Widget
        self.showexcel = QTableWidget(self)
        self.showexcel.setGeometry(QtCore.QRect(20, 20, 761, 151))
        self.showexcel.setColumnCount(15)
        self.showexcel.setRowCount(12)

        # Graphics View
        self.show_graphicsView = QGraphicsView(self)
        self.show_graphicsView.setGeometry(QtCore.QRect(340, 190, 441, 311))

        # Line Edits
        self.line_activedir_step1 = QLineEdit(self)
        self.line_activedir_step1.setGeometry(QtCore.QRect(20, 240, 191, 25))

        self.line_decoydir_step1 = QLineEdit(self)
        self.line_decoydir_step1.setGeometry(QtCore.QRect(20, 200, 191, 25))

        self.line_outdir_step1 = QLineEdit(self)
        self.line_outdir_step1.setGeometry(QtCore.QRect(20, 280, 191, 25))

        self.line_decoycsv_2 = QLineEdit(self)
        self.line_decoycsv_2.setGeometry(QtCore.QRect(20, 360, 191, 25))

        self.line_activecsv_2 = QLineEdit(self)
        self.line_activecsv_2.setGeometry(QtCore.QRect(20, 400, 191, 25))

        self.line_outdir_2 = QLineEdit(self)
        self.line_outdir_2.setGeometry(QtCore.QRect(20, 440, 191, 25))
        #Buttons
        self.button_close = QPushButton("Close", self)
        self.button_close.setGeometry(QtCore.QRect(330, 510, 80, 25))
        self.button_close.clicked.connect(self.close)
        
        self.Button_decoy_dir_step1 = QPushButton("decoy_dir", self)
        self.Button_decoy_dir_step1.setGeometry(QtCore.QRect(210, 200, 111, 25))

        self.button_acitvedir_step1 = QPushButton("active_dir", self)
        self.button_acitvedir_step1.setGeometry(QtCore.QRect(210, 240, 111, 25))

        self.button_outdir_step1 = QPushButton("out_dir", self)
        self.button_outdir_step1.setGeometry(QtCore.QRect(210, 280, 111, 25))

        self.button_summarizeresult_step1 = QPushButton("run_summarizeresult", self)
        self.button_summarizeresult_step1.setGeometry(QtCore.QRect(20, 320, 301, 25))

        self.button_runanalysis = QPushButton("run_Analysis", self)
        self.button_runanalysis.setGeometry(QtCore.QRect(20, 480, 301, 25))
     
        self.button_decoycsv_2 =  QPushButton("decoycsv", self)
        self.button_decoycsv_2.setGeometry(QtCore.QRect(210, 360, 111, 25))

        self.button_activecsv_2 = QPushButton("activecsv", self)
        self.button_activecsv_2.setGeometry(QtCore.QRect(210, 400, 111, 25))

        self.button_outdir_2 = QPushButton("plotoutdir", self)
        self.button_outdir_2.setGeometry(QtCore.QRect(210, 440, 111, 25))

        # Setup Connections
        self.setupConnections()
        self.setupLineEditWidgets()
        self.displayLogo()
        self.button_summarizeresult_step1.clicked.connect(self.SummarizeClicked)

### analysis_plot function
        self.button_runanalysis.clicked.connect(self.analysis_Clicked)

###############
    def setupLineEditWidgets(self):
        self.line_activedir_step1.setPlaceholderText("Path to active_dir")
        self.line_decoydir_step1.setPlaceholderText("Path to decoy_dir")
        self.line_outdir_step1.setPlaceholderText("Choose out dir")
        self.line_decoycsv_2.setPlaceholderText("Path to CSV for decoy")
        self.line_activecsv_2.setPlaceholderText("Path to CSV for active")
        self.line_outdir_2.setPlaceholderText("Choose anal_plot outdir")

    def setupConnections(self):
        self.Button_decoy_dir_step1.clicked.connect(self.choosedecoy_dir_step1)
        self.button_acitvedir_step1.clicked.connect(self.chooseacitvedir_step1)
        self.button_outdir_step1.clicked.connect(self.chooseoutDir1)
        self.button_decoycsv_2.clicked.connect(lambda: self.chooseFile(self.line_decoycsv_2, "CSV Files (*.csv)"))
        self.button_activecsv_2.clicked.connect(lambda: self.chooseFile(self.line_activecsv_2, "CSV Files (*.csv)"))
        self.button_outdir_2.clicked.connect(self.chooseoutDir2)

    def chooseoutDir1(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.line_outdir_step1.setText(dir)  
     
    def choosedecoy_dir_step1(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.line_decoydir_step1.setText(dir)
     
    def chooseacitvedir_step1(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.line_activedir_step1.setText(dir)

    def chooseoutDir2(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir:
            self.line_outdir_2.setText(dir)  

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
        super(AnalDialog, self).resizeEvent(event)
        self.displayLogo()

    def SummarizeClicked(self):
        decoy_file = self.line_decoydir_step1.text().strip()
        active_file = self.line_activedir_step1.text().strip()
        outdir = self.line_outdir_step1.text().strip()


        # 验证输入
        if not decoy_file or not active_file or not outdir:
            QMessageBox.warning(self, "Input Error", "Please fill all the fields.")
            return

        # 显示加载GIF
        self.displayLoadingGif()

        # 初始化并启动下载线程
        self.sum_thread = ScriptThread_sum(decoy_file, active_file, outdir)
        self.sum_thread.started.connect(self.onScriptStarted)
        self.sum_thread.finished.connect(self.onScriptFinished)  # 完成时连接到onScriptFinished槽
        self.sum_thread.start()

    def onScriptStarted(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished(self, out_dir):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()

        # 显示完成信息
        QMessageBox.information(self, "Script Finished", f"The summarize result script has completed. Output directory: {out_dir}")

        # 使用mol_name变量构建CSV文件路径
        csv_file_path = os.path.join(out_dir, f"inh_results.csv")
        self.displayCSV(csv_file_path)

#######analplot function#####################################
    def analysis_Clicked(self):
        decoy_csv = self.line_decoycsv_2.text().strip()
        active_csv = self.line_activecsv_2.text().strip()
        plot_out_dir = self.line_outdir_2.text().strip()
            # Debug prints
        print(f"'{decoy_csv}', '{active_csv}', '{plot_out_dir}'")

            # 验证输入
        if not decoy_csv or not active_csv or not plot_out_dir:
            QMessageBox.warning(self, "Input Error", "Please fill all the fields.")
            return

            # 显示加载GIF
        self.displayLoadingGif()

        # 初始化并启动
        self.filter_thread = ScriptThread_plot(decoy_csv, active_csv, plot_out_dir)
        self.filter_thread.started.connect(self.onScriptStarted_plot)
        self.filter_thread.finished.connect(self.onScriptFinished_plot)
        self.filter_thread.start()

###  A new onScriptFinished function 
    def onScriptStarted_plot(self):
        # 开始并显示GIF
        self.loadingLabel.show()
        self.loadingLabel.movie().start()

    def onScriptFinished_plot(self, plot_out_dir):
        # 停止并隐藏GIF
        self.loadingLabel.movie().stop()
        self.loadingLabel.hide()
        QMessageBox.information(self, "Analysis Complete", f"The Plots have been saved to {plot_out_dir}.")

        # 显示图片
        image_file_path1 = os.path.join(plot_out_dir, f"roc_plot.png")
        self.displayImage(image_file_path1)

        # 使用mol_name变量构建CSV文件路径
        csv_file_path = os.path.join(plot_out_dir, f"ttest_results.csv")
        self.displayCSV(csv_file_path)

        # image_file_path2 = os.path.join(plot_out_dir, f"boxplot.png")
        # self.displayImage(image_file_path2)

#############################        
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
            rect = self.show_graphicsView.rect()
            
            # 缩放图片以适应QGraphicsView的尺寸，同时保持纵横比
            scaledPixmap = pixmap.scaled(rect.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            
            # 创建QGraphicsPixmapItem并使用缩放后的图片
            item = QGraphicsPixmapItem(scaledPixmap)
            
            # 创建QGraphicsScene并添加item
            scene = QtWidgets.QGraphicsScene(self)
            scene.addItem(item)
            
            # 将创建的scene设置给QGraphicsView
            self.show_graphicsView.setScene(scene)
        else:
            QMessageBox.warning(self, "Display Error", "Failed to load the image.")

    def displayLoadingGif(self):
        # 如果已有，则先清除现有的场景
        if hasattr(self, 'loadingLabel') and self.loadingLabel:
            self.loadingLabel.deleteLater()
        
        # 在QGraphicsView中显示加载GIF
        self.loadingLabel = QtWidgets.QLabel(self.show_graphicsView)
        self.loadingLabel.setAlignment(Qt.AlignCenter)
        self.loadingLabel.setGeometry(0, 0, self.show_graphicsView.width(), self.show_graphicsView.height())
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
            
        

if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = AnalDialog()
    stylesheet_path = "style.qss"
    Dialog.setWindowFlags(Qt.FramelessWindowHint)
    AnalDialog.apply_stylesheet(app, stylesheet_path)
    Dialog.show()
    sys.exit(app.exec_())
