import os
from PyQt5.QtCore import Qt, QRectF
from PyQt5.QtGui import QPixmap, QPainter, QColor, QBrush, QPainterPath, QLinearGradient
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel

from qfluentwidgets import ScrollArea, isDarkTheme, FluentIcon
from ..common.config import cfg, HELP_URL, REPO_URL, EXAMPLE_URL, FEEDBACK_URL
from ..components.link_card import LinkCardView
from ..components.sample_card import SampleCardView
from ..common.style_sheet import StyleSheet

from ..common.signal_bus import signalBus

from .SoftGUI.Vina_GUI import vinadialogs
from .SoftGUI.Ledock_GUI import ledockdialogs
from .SoftGUI.Glide_GUI import glidedialogs
from .SoftGUI.AutodockGPU_GUI import autodockgpudialogs
from .SoftGUI.KarmaDock_GUI import KarmaDock_Dialog
from .SoftGUI.Analysis_UI import AnalDialog
from .SoftGUI.Evalu_GUI import TabbedGUI
from .SoftGUI.Virtual_ScreeningGUI import TabbedGUI2
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


class HomeInterface2(ScrollArea):
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
        self.setObjectName('homeInterface2')
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

        Vinaplace = os.getcwd() + "/images/controls/Button7.png"
        Ledockplace =os.getcwd() + "/images/controls/Button8.png"
        Glideplace=os.getcwd() + "/images/controls/Button10.png"
        AutodockGPUplace=os.getcwd() + "/images/controls/Button11.png"
        Tabbedplace=os.getcwd() + "/images/controls/Button13.png"
        virtualscreeningplace=os.getcwd() + "/images/controls/Button14.png"
        KarmaDockplace=os.getcwd() + "/images/controls/Button15.png"
        Analplace=os.getcwd() + "/images/controls/Button16.png"

        # layout samples
        layoutView = SampleCardView(self.tr('Docking Tools For Evaluation'), self.view)
        layoutView.addSampleCard(
            icon=Vinaplace,
            title="Vina",
            content=self.tr(
                "A button to call the vina function."),
            routeKey="layoutInterface",
            index=20
        )

        layoutView.addSampleCard(
            icon=KarmaDockplace,
            title="KarmaDock",
            content=self.tr(
                "A button to call the KarmaDock function."),
            routeKey="layoutInterface",
            index=26
        )        

        layoutView.addSampleCard(
            icon=Ledockplace,
            title="LeDock",
            content=self.tr(
                "A button to call the leDock function."),
            routeKey="layoutInterface",
            index=21
        )


        # layoutView.addSampleCard(
        #     icon=Glideplace,
        #     title="Glide_Evaluation",
        #     content=self.tr(
        #         "A button to call the Glide function."),
        #     routeKey="layoutInterface",
        #     index=22
        # )



        layoutView.addSampleCard(
            icon=AutodockGPUplace,
            title="AutodockGPU",
            content=self.tr(
                "A button to call the AutoDockGPU function."),
            routeKey="layoutInterface",
            index=23
        )

        # layoutView.addSampleCard(
        #     icon=Tabbedplace,
        #     title="Evaluation_workflow",
        #     content=self.tr(
        #         "A button to call the Evaluation function."),
        #     routeKey="layoutInterface",
        #     index=24
        # )

        # layoutView.addSampleCard(
        #     icon=virtualscreeningplace,
        #     title="VirtualScreening_Evaluation",
        #     content=self.tr(
        #         "A button to call the VirtualScreening function."),
        #     routeKey="layoutInterface",
        #     index=25
        # )

        self.vBoxLayout.addWidget(layoutView)



        # material samples
        materialView = SampleCardView(self.tr('Integrated Module'), self.view)
        materialView.addSampleCard(
            icon=Tabbedplace,
            title="Evaluation_Workflow",
            content=self.tr(
                "A button to call the Evaluation function."),
            routeKey="materialInterface",
            index=24
        )

        materialView.addSampleCard(
            icon=virtualscreeningplace,
            title="Virtual_Screening",
            content=self.tr(
                "A button to call the VirtualScreening function."),
            routeKey="materialInterface",
            index=25
        )
        self.vBoxLayout.addWidget(materialView)

        materialView.addSampleCard(
            icon=Analplace,
            title="Analysis_Evaluation",
            content=self.tr(
                "A button to call the AnalysisEvaluation function."),
            routeKey="materialInterface",
            index=27
        )
        self.vBoxLayout.addWidget(materialView)


        # # menu samples
        # menuView = SampleCardView(self.tr('Menu & toolbars samples'), self.view)
        # menuView.addSampleCard(
        #     icon=":/gallery/images/controls/MenuFlyout.png",
        #     title="RoundMenu",
        #     content=self.tr(
        #         "Shows a contextual list of simple commands or options."),
        #     routeKey="menuInterface",
        #     index=0
        # )
        # menuView.addSampleCard(
        #     icon=":/gallery/images/controls/CommandBar.png",
        #     title="CommandBar",
        #     content=self.tr(
        #         "Shows a contextual list of simple commands or options."),
        #     routeKey="menuInterface",
        #     index=2
        # )
        # menuView.addSampleCard(
        #     icon=":/gallery/images/controls/CommandBarFlyout.png",
        #     title="CommandBarFlyout",
        #     content=self.tr(
        #         "A mini-toolbar displaying proactive commands, and an optional menu of commands."),
        #     routeKey="menuInterface",
        #     index=3
        # )
        # self.vBoxLayout.addWidget(menuView)

        # # navigation
        # navigationView = SampleCardView(self.tr('Navigation'), self.view)
        # navigationView.addSampleCard(
        #     icon=":/gallery/images/controls/BreadcrumbBar.png",
        #     title="BreadcrumbBar",
        #     content=self.tr(
        #         "Shows the trail of navigation taken to the current location."),
        #     routeKey="navigationViewInterface",
        #     index=0
        # )
        # navigationView.addSampleCard(
        #     icon=":/gallery/images/controls/Pivot.png",
        #     title="Pivot",
        #     content=self.tr(
        #         "Presents information from different sources in a tabbed view."),
        #     routeKey="navigationViewInterface",
        #     index=1
        # )
        # navigationView.addSampleCard(
        #     icon=":/gallery/images/controls/TabView.png",
        #     title="TabView",
        #     content=self.tr(
        #         "Presents information from different sources in a tabbed view."),
        #     routeKey="navigationViewInterface",
        #     index=3
        # )
        # self.vBoxLayout.addWidget(navigationView)

        # # scroll samples
        # scrollView = SampleCardView(self.tr('Scrolling samples'), self.view)
        # scrollView.addSampleCard(
        #     icon=":/gallery/images/controls/ScrollViewer.png",
        #     title="ScrollArea",
        #     content=self.tr(
        #         "A container control that lets the user pan and zoom its content smoothly."),
        #     routeKey="scrollInterface",
        #     index=0
        # )
        # scrollView.addSampleCard(
        #     icon=":/gallery/images/controls/PipsPager.png",
        #     title="PipsPager",
        #     content=self.tr(
        #         "A control to let the user navigate through a paginated collection when the page numbers do not need to be visually known."),
        #     routeKey="scrollInterface",
        #     index=3
        # )
        # self.vBoxLayout.addWidget(scrollView)

        # # state info samples
        # stateInfoView = SampleCardView(self.tr('Status & info samples'), self.view)
        # stateInfoView.addSampleCard(
        #     icon=":/gallery/images/controls/ProgressRing.png",
        #     title="StateToolTip",
        #     content=self.tr(
        #         "Shows the apps progress on a task, or that the app is performing ongoing work that does block user interaction."),
        #     routeKey="statusInfoInterface",
        #     index=0
        # )
        # stateInfoView.addSampleCard(
        #     icon=":/gallery/images/controls/InfoBadge.png",
        #     title="InfoBadge",
        #     content=self.tr(
        #         "An non-intrusive Ul to display notifications or bring focus to an area."),
        #     routeKey="statusInfoInterface",
        #     index=3
        # )
        # stateInfoView.addSampleCard(
        #     icon=":/gallery/images/controls/InfoBar.png",
        #     title="InfoBar",
        #     content=self.tr(
        #         "An inline message to display app-wide status change information."),
        #     routeKey="statusInfoInterface",
        #     index=4
        # )
        # stateInfoView.addSampleCard(
        #     icon=":/gallery/images/controls/ProgressBar.png",
        #     title="ProgressBar",
        #     content=self.tr(
        #         "Shows the apps progress on a task, or that the app is performing ongoing work that doesn't block user interaction."),
        #     routeKey="statusInfoInterface",
        #     index=8
        # )
        # stateInfoView.addSampleCard(
        #     icon=":/gallery/images/controls/ProgressRing.png",
        #     title="ProgressRing",
        #     content=self.tr(
        #         "Shows the apps progress on a task, or that the app is performing ongoing work that doesn't block user interaction."),
        #     routeKey="statusInfoInterface",
        #     index=10
        # )
        # stateInfoView.addSampleCard(
        #     icon=":/gallery/images/controls/ToolTip.png",
        #     title="ToolTip",
        #     content=self.tr(
        #         "Displays information for an element in a pop-up window."),
        #     routeKey="statusInfoInterface",
        #     index=1
        # )
        # self.vBoxLayout.addWidget(stateInfoView)

        # # text samples
        # textView = SampleCardView(self.tr('Text samples'), self.view)
        # textView.addSampleCard(
        #     icon=":/gallery/images/controls/TextBox.png",
        #     title="LineEdit",
        #     content=self.tr("A single-line plain text field."),
        #     routeKey="textInterface",
        #     index=0
        # )
        # textView.addSampleCard(
        #     icon=":/gallery/images/controls/PasswordBox.png",
        #     title="PasswordLineEdit",
        #     content=self.tr("A control for entering passwords."),
        #     routeKey="textInterface",
        #     index=2
        # )
        # textView.addSampleCard(
        #     icon=":/gallery/images/controls/NumberBox.png",
        #     title="SpinBox",
        #     content=self.tr(
        #         "A text control used for numeric input and evaluation of algebraic equations."),
        #     routeKey="textInterface",
        #     index=3
        # )
        # textView.addSampleCard(
        #     icon=":/gallery/images/controls/RichEditBox.png",
        #     title="TextEdit",
        #     content=self.tr(
        #         "A rich text editing control that supports formatted text, hyperlinks, and other rich content."),
        #     routeKey="textInterface",
        #     index=8
        # )
        # self.vBoxLayout.addWidget(textView)

        # # view samples
        # collectionView = SampleCardView(self.tr('View samples'), self.view)
        # collectionView.addSampleCard(
        #     icon=":/gallery/images/controls/ListView.png",
        #     title="ListView",
        #     content=self.tr(
        #         "A control that presents a collection of items in a vertical list."),
        #     routeKey="viewInterface",
        #     index=0
        # )
        # collectionView.addSampleCard(
        #     icon=":/gallery/images/controls/DataGrid.png",
        #     title="TableView",
        #     content=self.tr(
        #         "The DataGrid control provides a flexible way to display a collection of data in rows and columns."),
        #     routeKey="viewInterface",
        #     index=1
        # )
        # collectionView.addSampleCard(
        #     icon=":/gallery/images/controls/TreeView.png",
        #     title="TreeView",
        #     content=self.tr(
        #         "The TreeView control is a hierarchical list pattern with expanding and collapsing nodes that contain nested items."),
        #     routeKey="viewInterface",
        #     index=2
        # )
        # collectionView.addSampleCard(
        #     icon=":/gallery/images/controls/FlipView.png",
        #     title="FlipView",
        #     content=self.tr(
        #         "Presents a collection of items that the user can flip through,one item at a time."),
        #     routeKey="viewInterface",
        #     index=4
        # )
        # self.vBoxLayout.addWidget(collectionView)
    def on_sample_card_clicked(self, routeKey, index):
        if routeKey == "layoutInterface" and index == 20:  # Adjust these conditions as necessary
            self.open_vina_dialog()
        if routeKey == "layoutInterface" and index == 21:  # Adjust these conditions as necessary
            self.open_ledock_dialog()
        if routeKey == "layoutInterface" and index == 22:  # Adjust these conditions as necessary
            self.open_glide_dialog()
        if routeKey == "layoutInterface" and index == 23:  # Adjust these conditions as necessary
            self.open_autodockGPU_dialog()
        if routeKey == "materialInterface" and index == 24:  # Adjust these conditions as necessary
            self.open_workflow_dialog()
        if routeKey == "materialInterface" and index == 25:  # Adjust these conditions as necessary
            self.open_virtualscreening_dialog()
        if routeKey == "layoutInterface" and index == 26:  # Adjust these conditions as necessary
            self.open_KarmaDock_Dialog()
        if routeKey == "materialInterface" and index == 27:  # Adjust these conditions as necessary
            self.open_AnalDialog()

    def open_vina_dialog(self):
        dialog = vinadialogs(self)
        dialog.exec_()
    
    def open_ledock_dialog(self):
        dialog = ledockdialogs(self)
        dialog.exec_()
    
    def open_glide_dialog(self):
        dialog = glidedialogs(self)
        dialog.exec_()

    def open_autodockGPU_dialog(self):
        dialog = autodockgpudialogs(self)
        dialog.exec_()
    
    def open_workflow_dialog(self):
        self.TabbedGUI = TabbedGUI(self)
        self.TabbedGUI.show()

    def open_virtualscreening_dialog(self):
        self.TabbedGUI2 = TabbedGUI2(self)
        self.TabbedGUI2.show()

    def open_KarmaDock_Dialog(self):
        dialog = KarmaDock_Dialog(self)
        dialog.show()

    def open_AnalDialog(self):
        dialog = AnalDialog(self)
        dialog.show()

    def open_BK_dialog(self):
        dialog = BKdialogs(self)
        dialog.exec_()