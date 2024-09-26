import os
from PyQt5.QtCore import Qt, QRectF
from PyQt5.QtGui import QPixmap, QPainter, QColor, QBrush, QPainterPath, QLinearGradient
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel

from qfluentwidgets import ScrollArea, isDarkTheme, FluentIcon
from ..common.config import cfg, HELP_URL, REPO_URL, EXAMPLE_URL, FEEDBACK_URL
from ..components.link_card import LinkCardView
from ..components.sample_card import SampleCardView
from ..common.style_sheet import StyleSheet

from .SoftGUI.ligdownGUI import InputligdownDialogs
from .SoftGUI.ligprepUI import LigprepDialog
from .SoftGUI.prodownGUI import InputprodownDialogs
from .SoftGUI.deeppocket_GUI import DeeppocketDialog
from .SoftGUI.ligclusGUI import InputclusDialogs
from ..common.signal_bus import signalBus
from .SoftGUI.lig3DforGlideGUI import Input3DligGDialogs
from .SoftGUI.lig3DGUI import Input3DligDialogs
from .SoftGUI.DecoyGUI import decoydialogs
from .SoftGUI.Vina_GUI import vinadialogs
from .SoftGUI.Ledock_GUI import ledockdialogs
from .SoftGUI.Glide_GUI import glidedialogs
from .SoftGUI.AutodockGPU_GUI import autodockgpudialogs
from .SoftGUI.Evalu_GUI import TabbedGUI
from .SoftGUI.BK_GUI import BKdialogs
class BannerWidget(QWidget):
    """ Banner widget """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setFixedHeight(336)
     
        self.vBoxLayout = QVBoxLayout(self)
        self.galleryLabel = QLabel('EvaluationMaster', self)
        current_dir = os.getcwd()
        img_dir = current_dir + "/images/header1.png"
        self.banner = QPixmap(img_dir)
        if self.banner.isNull():
            print("Failed to load image:", img_dir)
        self.linkCardView = LinkCardView(self)

        self.galleryLabel.setObjectName('galleryLabel')

        self.vBoxLayout.setSpacing(0)
        self.vBoxLayout.setContentsMargins(0, 20, 0, 0)
        self.vBoxLayout.addWidget(self.galleryLabel)
        self.vBoxLayout.addWidget(self.linkCardView, 1, Qt.AlignBottom)
        self.vBoxLayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        logo_dir= os.getcwd()
        self.linkCardView.addCard(
            f'{logo_dir}/images/logo.png',
            self.tr('Getting started'),
            self.tr('An overview of app development.'),
            HELP_URL
        )

        self.linkCardView.addCard(
            FluentIcon.GITHUB,
            self.tr('GitHub repo'),
            self.tr(
                'The latest EvaluationMaster Code.'),
            REPO_URL
        )

        self.linkCardView.addCard(
            FluentIcon.CODE,
            self.tr('Code samples'),
            self.tr(
                'Find samples that demonstrate specific features.'),
            EXAMPLE_URL
        )

        self.linkCardView.addCard(
            FluentIcon.FEEDBACK,
            self.tr('Send feedback'),
            self.tr('Help us improve EvaluationMaster by providing feedback.'),
            FEEDBACK_URL
        )


    def paintEvent(self, e):
        super().paintEvent(e)
        painter = QPainter(self)
        painter.setRenderHints(
            QPainter.SmoothPixmapTransform | QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)

        path = QPainterPath()
        path.setFillRule(Qt.WindingFill)
        w, h = self.width(), self.height()
        path.addRoundedRect(QRectF(0, 0, w, h), 10, 10)
        path.addRect(QRectF(0, h-50, 50, 50))
        path.addRect(QRectF(w-50, 0, 50, 50))
        path.addRect(QRectF(w-50, h-50, 50, 50))
        path = path.simplified()

        # init linear gradient effect
        gradient = QLinearGradient(0, 0, 0, h)

        # draw background color
        if not isDarkTheme():
            gradient.setColorAt(0, QColor(207, 216, 228, 255))
            gradient.setColorAt(1, QColor(207, 216, 228, 0))
        else:
            gradient.setColorAt(0, QColor(0, 0, 0, 255))
            gradient.setColorAt(1, QColor(0, 0, 0, 0))
            
        painter.fillPath(path, QBrush(gradient))

        # draw banner image
        pixmap = self.banner.scaled(
            self.size(), transformMode=Qt.SmoothTransformation)
        painter.fillPath(path, QBrush(pixmap))


class HomeInterface(ScrollArea):
    """ Home interface """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.banner = BannerWidget(self)
        self.view = QWidget(self)
        self.vBoxLayout = QVBoxLayout(self.view)

        self.__initWidget()
        self.loadSamples()

        # Connect signal for SampleCard clicks
        signalBus.switchToSampleCard.connect(self.on_sample_card_clicked)

    def __initWidget(self):
        self.view.setObjectName('view')
        self.setObjectName('homeInterface')
        StyleSheet.HOME_INTERFACE.apply(self)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setWidget(self.view)
        self.setWidgetResizable(True)

        self.vBoxLayout.setContentsMargins(0, 0, 0, 36)
        self.vBoxLayout.setSpacing(40)
        self.vBoxLayout.addWidget(self.banner)
        self.vBoxLayout.setAlignment(Qt.AlignTop)

    def loadSamples(self):
        """ load samples """


        icoplace = os.getcwd() + "/images/controls/Button1.png"
        proicoplace = os.getcwd() + "/images/controls/Button2.png"
        LigClusplace = os.getcwd() + "/images/controls/Button3.png"
        Lig3DforGlideplace = os.getcwd() + "/images/controls/Button4.png"
        Lig3Dplace = os.getcwd() + "/images/controls/Button5.png"
        Decoyplace = os.getcwd() + "/images/controls/Button6.png"
        Vinaplace = os.getcwd() + "/images/controls/Button7.png"
        Ledockplace =os.getcwd() + "/images/controls/Button8.png"
        sitemapplace=os.getcwd() + "/images/controls/Button9.png"
        Glideplace=os.getcwd() + "/images/controls/Button10.png"
        AutodockGPUplace=os.getcwd() + "/images/controls/Button11.png"
        Tabbedplace=os.getcwd() + "/images/controls/Button11.png"
        ADCplace=os.getcwd() + "/images/controls/Button13.png"
        basicInputView = SampleCardView(
            self.tr("Protein Preparation Part"), self.view)

        basicInputView.addSampleCard(
            icon=proicoplace ,
            title="Protein_Download",
            content=self.tr("This Button streamlines the extraction and analysis of Protein structure data from RCSB_PDB."),
            routeKey="basicInputInterface",
            index=1
        )

        basicInputView.addSampleCard(
            icon=sitemapplace,
            title="Protein Binding Pocket Prediction",
            content=self.tr(
                "This Button control the Protein Binding Pocket Prediction function"),
            routeKey="basicInputInterface",
            index=2
        )


        self.vBoxLayout.addWidget(basicInputView)



        # date time samples
        dateTimeView = SampleCardView(self.tr('Molecules Handle'), self.view)

        dateTimeView.addSampleCard(

            # icon=":/gallery/images/controls/Button.png",
            icon=icoplace,
            title="Ligand_Download_and_Preparation",
            content=self.tr(
                "This button streamlines the extraction and analysis of bioactive compound data from ChEMBL."),
            routeKey="dateTimeInterface",  # Ensure this matches with the intended routeKey
            index=0  # Ensure this matches with the intended index
        )


        dateTimeView.addSampleCard(
            icon=LigClusplace,
            title="Ligand Cluster",
            content=self.tr(
                "This Button control the ligands clus function"),
            routeKey="dateTimeInterface",
            index=10
        )

        dateTimeView.addSampleCard(
            icon=Decoyplace,
            title="Decoy_Generation",
            content=self.tr("The function was used to generate decoys."),
            routeKey="dateTimeInterface",
            index=13
        )


        dateTimeView.addSampleCard(
            icon=Lig3Dplace,
            title=" Ligand 3D conformation Generation",
            content=self.tr("Ligand 3D conformation generation for sdf, mol2 and pdbqt format"),
            routeKey="dateTimeInterface",
            index=12
        )

        self.vBoxLayout.addWidget(dateTimeView)




    def on_sample_card_clicked(self, routeKey, index):
        if routeKey == "dateTimeInterface" and index == 0:  # Adjust these conditions as necessary
            self.open_inputligdown_dialog()
        if routeKey == "basicInputInterface" and index == 1:  # Adjust these conditions as necessary
            self.open_inputprodown_dialog()


        if routeKey == "basicInputInterface" and index == 2:  # Adjust these conditions as necessary
            self.DeeppocketDialog_dialog()

            
        if routeKey == "dateTimeInterface" and index == 10:  # Adjust these conditions as necessary
            self.open_inputclus_dialog()
        if routeKey == "dateTimeInterface" and index == 11:  # Adjust these conditions as necessary
            self.open_input3DligG_dialog()
        if routeKey == "dateTimeInterface" and index == 12:  # Adjust these conditions as necessary
            self.open_input3Dlig_dialog()
        if routeKey == "dateTimeInterface" and index == 13:  # Adjust these conditions as necessary
            self.open_decoy_dialog()


    def open_inputligdown_dialog(self):
        dialog = LigprepDialog(self)
        dialog.show()

    def DeeppocketDialog_dialog(self):
        dialog = DeeppocketDialog(self)
        dialog.show()
        
    def open_inputprodown_dialog(self):
        dialog = InputprodownDialogs(self)
        dialog.show()

    def open_inputclus_dialog(self):
        dialog = InputclusDialogs(self)
        dialog.exec_()

    def open_input3DligG_dialog(self):
        dialog = Input3DligGDialogs(self)
        dialog.exec_()

    def open_input3Dlig_dialog(self):
        dialog = Input3DligDialogs(self)
        dialog.exec_()
    
    def open_decoy_dialog(self):
        dialog = decoydialogs(self)
        dialog.exec_()

