# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'interface.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1224, 944)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.lineEdit = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit.setGeometry(QtCore.QRect(10, 30, 113, 27))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.pushButton = QtGui.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(10, 60, 99, 27))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.verticalLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(170, 20, 1154, 1770))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.jusantepng = QtGui.QLabel(self.verticalLayoutWidget)
        self.jusantepng.setEnabled(True)
        self.jusantepng.setMinimumSize(QtCore.QSize(1152, 860))
        self.jusantepng.setMaximumSize(QtCore.QSize(1152, 860))
        self.jusantepng.setFrameShape(QtGui.QFrame.Box)
        self.jusantepng.setText(_fromUtf8(""))
        self.jusantepng.setTextFormat(QtCore.Qt.RichText)
        self.jusantepng.setPixmap(QtGui.QPixmap(_fromUtf8("interface_figs/jusante.png")))
        self.jusantepng.setScaledContents(True)
        self.jusantepng.setWordWrap(False)
        self.jusantepng.setObjectName(_fromUtf8("jusantepng"))
        self.verticalLayout.addWidget(self.jusantepng)
        self.montantepng = QtGui.QLabel(self.verticalLayoutWidget)
        self.montantepng.setText(_fromUtf8(""))
        self.montantepng.setPixmap(QtGui.QPixmap(_fromUtf8("interface_figs/montante.png")))
        self.montantepng.setObjectName(_fromUtf8("montantepng"))
        self.verticalLayout.addWidget(self.montantepng)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1224, 25))
        self.menuBar.setDefaultUp(False)
        self.menuBar.setNativeMenuBar(False)
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menuBar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.menuFile.addAction(self.actionOpen)
        self.menuBar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.lineEdit.setText(_translate("MainWindow", "grid number", None))
        self.pushButton.setText(_translate("MainWindow", "coat", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionOpen.setText(_translate("MainWindow", "Open", None))

