import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QLabel, QPushButton, QCheckBox
import subprocess
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
import os
import pandas as pd
from PyQt5.QtCore import QThread, pyqtSignal




# class ScriptligprepGThread(QThread):
#     started = pyqtSignal()
#     finished = pyqtSignal(str, str)

#     def __init__(self, njobs, input_file , save_dir, Schro_dir):
#         super().__init__()
#         self.njobs = njobs
#         self.input_file = input_file
#         self.save_dir = save_dir
#         self.Schro_dir = Schro_dir


#     def run(self):
#         self.started.emit()
#         current_path = os.getcwd()
#         script_path = current_path + "/app/Script/LigScript/ligprepG.py"
#         subprocess.call(["python", script_path, self.njobs, self.input_file, self.save_dir, self.Schro_dir])
#         self.finished.emit(self.njobs, self.save_dir)  # Emit with parameters

class ScriptThreadadgpu(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, csv_file, ligand_csv_file, lig_path, mgltool_path, gpu_file, save_path, gpu_num, n_run, index, preserve_num):
        super().__init__()
        self.csv_file = csv_file
        self.ligand_csv_file = ligand_csv_file
        self.lig_path = lig_path
        self.mgltool_path = mgltool_path
        self.gpu_file = gpu_file
        self.save_path = save_path
        self.gpu_num = gpu_num
        self.n_run = n_run
        self.index = index
        self.preserve_num = preserve_num


    def run(self):
        self.started.emit()
        if self.index > 1:  # 当在第2、3、4阶段时
            temp_csv_file = f"{self.save_path}/stage{self.index-1}/results.csv"
            self.save_path = f"{self.save_path}/stage{self.index}/"
            df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
            num_rows_to_keep = int(len(df) * (self.preserve_num / 100.0))
            df_subset = df.head(num_rows_to_keep)
            df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
            self.ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
            # lig_start_name = "lig"
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/ScreeningScript/Gen3D.py"
        subprocess.call(["python", script_path,  self.ligand_csv_file, self.lig_path , self.mgltool_path])
        self.lig_path = self.lig_path+ "/PDBQT/"
        script_path2 = current_path + "/app/Script/ScreeningScript/AutodockGPU.py"
        subprocess.call([f"python", script_path2, self.csv_file, self.lig_path, self.mgltool_path, self.gpu_file, self.save_path, str(self.gpu_num), str(self.n_run)])
        self.finished.emit(self.csv_file)  # Emit with parameter


class ScriptThreadledock(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, total_threads, csv_file,ligand_csv_file, lig_path, save_path, ledock_path,index, preserve_num):
        super().__init__()
        self.total_threads = total_threads
        self.csv_file = csv_file
        self.ligand_csv_file = ligand_csv_file
        self.lig_path = lig_path
        self.save_path = save_path
        self.ledock_path = ledock_path
        self.index = index
        self.preserve_num = preserve_num
        # if index > 1:  # 当在第2、3、4阶段时
        #     temp_csv_file = f"{save_path}/stage{index-1}/results.csv"
        #     self.save_path = f"{save_path}/stage{index}/"
        #     df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
        #     num_rows_to_keep = int(len(df) * (preserve_num / 100.0))
        #     df_subset = df.head(num_rows_to_keep)
        #     df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
        #     ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
    def run(self):
        self.started.emit()
        if self.index > 1:  # 当在第2、3、4阶段时
            temp_csv_file = f"{self.save_path}/stage{self.index-1}/results.csv"
            self.save_path = f"{self.save_path}/stage{self.index}/"
            df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
            num_rows_to_keep = int(len(df) * (self.preserve_num / 100.0))
            df_subset = df.head(num_rows_to_keep)
            df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
            self.ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
            # lig_start_name = "lig"
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/ScreeningScript/Gen3DMol2.py"
        subprocess.call(["python", script_path,  self.ligand_csv_file, self.lig_path])
        self.lig_path = self.lig_path+ "/MOL2/"
        script_path2 = current_path + "/app/Script/ScreeningScript/Ledock.py"
        subprocess.call([f"python", script_path2, self.total_threads, self.csv_file, self.lig_path, self.save_path, self.ledock_path])
        self.finished.emit(self.csv_file)  # Emit with parameters



class ScriptThreadvina(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, csv_file, ligand_csv_file, lig_path, tool_path, save_path, index, preserve_num):
        super().__init__()

        self.csv_file = csv_file
        self.lig_path = lig_path
        self.ligand_csv_file = ligand_csv_file
        self.tool_path = tool_path
        self.save_path = save_path
        self.index = index
        self.preserve_num = preserve_num
        # if index > 1:  # 当在第2、3、4阶段时
        #     temp_csv_file = f"{save_path}/stage{index-1}/results.csv"
        #     self.save_path = f"{save_path}/stage{index}/"
        #     df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
        #     num_rows_to_keep = int(len(df) * (preserve_num / 100.0))
        #     df_subset = df.head(num_rows_to_keep)
        #     df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
        #     ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
    def run(self):
        self.started.emit()
        if self.index > 1:  # 当在第2、3、4阶段时
            temp_csv_file = f"{self.save_path}/stage{self.index-1}/results.csv"
            self.save_path = f"{self.save_path}/stage{self.index}/"
            df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
            num_rows_to_keep = int(len(df) * (self.preserve_num / 100.0))
            df_subset = df.head(num_rows_to_keep)
            df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
            self.ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
            # lig_start_name = "lig"
        current_path = os.getcwd()
        script_path = current_path + "/app/Script/ScreeningScript/Gen3D.py"
        subprocess.call(["python", script_path,  self.ligand_csv_file, self.lig_path, self.tool_path])
        self.lig_path = self.lig_path+ "/PDBQT/"
        script_path2 = current_path + "/app/Script/ScreeningScript/AutoDock_Vina.py"
        subprocess.call([f"python", script_path2, self.csv_file, self.lig_path, self.tool_path, self.save_path])
        self.finished.emit(self.csv_file)  # Emit with parameters

# class ScriptThreadglide(QThread):
#     started = pyqtSignal()
#     finished = pyqtSignal(str)

#     def __init__(self, csv_file_path, ligand_csv_file, lig_file, schro_dir, mm_sahre_dir, save_dir, Glide_mode,threads_num, index, preserve_num):
#         super().__init__()
#         self.csv_file_path = csv_file_path
#         self.ligand_csv_file = ligand_csv_file
#         self.lig_file = lig_file
#         self.schro_dir = schro_dir
#         self.mm_sahre_dir = mm_sahre_dir
#         self.save_dir = save_dir
#         self.Glide_mode = Glide_mode
#         self.threads_num = threads_num
#         self.index = index
#         self.preserve_num = preserve_num
#         # if index > 1:  # 当在第2、3、4阶段时
#         #     temp_csv_file = f"{save_dir}/stage{index-1}/results.csv"
#         #     self.save_path = f"{save_dir}/stage{index}/"
#         #     df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
#         #     num_rows_to_keep = int(len(df) * (preserve_num / 100.0))
#         #     df_subset = df.head(num_rows_to_keep)
#         #     df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
#         #     ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径


#     def run(self):
#         self.started.emit()
#         if self.index > 1:  # 当在第2、3、4阶段时
#             temp_csv_file = f"{self.save_dir}/stage{self.index-1}/results.csv"
#             self.save_dir = f"{self.save_dir}/stage{self.index}/"
#             df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
#             num_rows_to_keep = int(len(df) * (self.preserve_num / 100.0))
#             df_subset = df.head(num_rows_to_keep)
#             df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
#             self.ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
#         current_path = os.getcwd()
#         script_path = current_path + "/app/Script/ScreeningScript/ligprepG.py"
#         subprocess.call([f"python", script_path, self.threads_num, self.ligand_csv_file, self.lig_file, self.schro_dir])
#         self.lig_file = self.lig_file + "/ligands.maegz"
#         script_path2 = current_path + "/app/Script/ScreeningScript/GlideG.py"
#         subprocess.call([f"python", script_path2, self.csv_file_path, self.lig_file, self.schro_dir, self.mm_sahre_dir, self.save_dir, self.Glide_mode, self.threads_num])
#         self.finished.emit(self.csv_file_path)  # Emit with parameters

class ScriptThreadglide(QThread):
    started = pyqtSignal()
    finished = pyqtSignal(str)

    def __init__(self, schro_dir, csv_file, file_path, ligand_csv_file, save_path,threads_num, index, preserve_num):
        super().__init__()
        self.schro_dir = schro_dir
        self.csv_file_path = csv_file
        self.file_path = file_path
        self.ligand_csv_file = ligand_csv_file
        self.save_dir = save_path
        self.threads_num = threads_num
        self.index = index
        self.preserve_num = preserve_num
        # if index > 1:  # 当在第2、3、4阶段时
        #     temp_csv_file = f"{save_dir}/stage{index-1}/results.csv"
        #     self.save_path = f"{save_dir}/stage{index}/"
        #     df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
        #     num_rows_to_keep = int(len(df) * (preserve_num / 100.0))
        #     df_subset = df.head(num_rows_to_keep)
        #     df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
        #     ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径


    def run(self):
        self.started.emit()
        if self.index > 1:  # 当在第2、3、4阶段时
            temp_csv_file = f"{self.save_dir}/stage{self.index-1}/results.csv"
            self.save_dir = f"{self.save_dir}/stage{self.index}/"
            df = pd.read_csv(temp_csv_file)  # 读取前一阶段的CSV文件或初始ligand_csv_file
            num_rows_to_keep = int(len(df) * (self.preserve_num / 100.0))
            df_subset = df.head(num_rows_to_keep)
            df_subset.to_csv(temp_csv_file, index=False)  # 保存截断后的数据到新的CSV文件
            self.ligand_csv_file = temp_csv_file  # 更新ligand_csv_file为包含截断数据的新文件路径
        current_path = os.getcwd()
        # script_path = current_path + "/app/Script/ScreeningScript/ligprepG.py"
        # subprocess.call([f"python", script_path, self.threads_num, self.ligand_csv_file, self.lig_file, self.schro_dir])
        # self.lig_file = self.lig_file + "/ligands.maegz"
        # script_path2 = current_path + "/app/Script/ScreeningScript/GlideG.py"
        # subprocess.call([f"python", script_path2, self.schro_dir, self.lig_file, self.schro_dir, self.mm_sahre_dir, self.save_dir, self.Glide_mode, self.threads_num])
        script_path = current_path + "/app/Script/ScreeningScript/Karmadock.py"
        subprocess.call ([f"{self.schro_dir}/Env/bin/python", "-u", script_path, self.schro_dir, self.csv_file_path, self.file_path, self.ligand_csv_file, self.save_dir, self.threads_num])
        self.finished.emit(self.csv_file_path)  # Emit with parameters



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = ScriptThreadglide()