# coding: utf-8
from PyQt5.QtCore import QUrl, QSize
from PyQt5.QtGui import QIcon, QDesktopServices
from PyQt5.QtWidgets import QApplication

from qfluentwidgets import (NavigationAvatarWidget, NavigationItemPosition, MessageBox, FluentWindow,
                            SplashScreen)
from qfluentwidgets import FluentIcon as FIF

from .gallery_interface import GalleryInterface
from .home_interface import HomeInterface

from .icon_interface import IconInterface

from .basic_input_interface import BasicInputInterface
from .date_time_interface import DateTimeInterface
from .dialog_interface import DialogInterface
from .layout_interface import LayoutInterface
from .home2_interface import HomeInterface2
from .material_interface import MaterialInterface
from .menu_interface import MenuInterface
from .navigation_view_interface import NavigationViewInterface
from .scroll_interface import ScrollInterface
from .status_info_interface import StatusInfoInterface
from .setting_interface import SettingInterface
from .text_interface import TextInterface
from .view_interface import ViewInterface
from ..common.config import SUPPORT_URL, cfg
from ..common.icon import Icon
from ..common.signal_bus import signalBus
from ..common.translator import Translator
from ..common import resource
import numpy as np
import os
class MainWindow(FluentWindow):

    def __init__(self):
        super().__init__()
        self.initWindow()

        # create sub interface
        # 创建原有的子界面
        self.homeInterface = HomeInterface(self)

        # 创建新的子界面
        self.homeInterface2 = HomeInterface2(self)  # 假设你已经创建了 HomeInterface2 类

        self.initNavigation()
        self.splashScreen.finish()


        self.initNavigation()
        self.splashScreen.finish()

    def initNavigation(self):
        t = Translator()
        # 添加原有的导航项
        self.addSubInterface(self.homeInterface, FIF.HOME, self.tr('Basic function'))

        # 添加新的导航项
        self.addSubInterface(self.homeInterface2, FIF.DATE_TIME, self.tr('Basic function2'))
        # 注意：如果你想使用不同的图标，可以更改 FIF.HOME 为其他图标


    def initWindow(self):
        logo_dir= os.getcwd()

        self.resize(960, 780)
        self.setMinimumWidth(760)
        self.setWindowIcon(QIcon(f'{logo_dir}/images/logo.png'))
        self.setWindowTitle('Evaluation_master')

        self.setMicaEffectEnabled(cfg.get(cfg.micaEnabled))

        # create splash screen
        self.splashScreen = SplashScreen(self.windowIcon(), self)
        self.splashScreen.setIconSize(QSize(106, 106))
        self.splashScreen.raise_()

        desktop = QApplication.desktop().availableGeometry()
        w, h = desktop.width(), desktop.height()
        self.move(w//2 - self.width()//2, h//2 - self.height()//2)
        self.show()
        QApplication.processEvents()


    def resizeEvent(self, e):
        super().resizeEvent(e)
        if hasattr(self, 'splashScreen'):
            self.splashScreen.resize(self.size())

    def switchToSample(self, routeKey, index):
        """ switch to sample """
        interfaces = self.findChildren(GalleryInterface)
        for w in interfaces:
            if w.objectName() == routeKey:
                self.stackedWidget.setCurrentWidget(w, False)
                w.scrollToCard(index)
