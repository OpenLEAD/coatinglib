from PyQt4 import QtCore, QtGui  # Import the PyQt4 module we'll need
import sys  # We need sys so that we can pass argv to QApplication

import blade_coverage

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

import qt_interface  # This file holds our MainWindow and all design related things

# it also keeps events etc that we defined in Qt Designer
import os  # For listing directory methods

dics = {
    "Faces": ['jusante','montante',blade_coverage.jusante()+blade_coverage.montante()],
    "Lados": ['borda_direita','borda_esquerda',blade_coverage.lip()],
    "Lips": ['lip','lip',blade_coverage.lip()],
    }

class ExampleApp(QtGui.QMainWindow, qt_interface.Ui_MainWindow):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
        # It sets up layout and widgets that are defined
        self.actionOpen.triggered.connect(self.browse_folder)  # When the button is presse
        self.dbs.addItems(dics.keys())
        self.dbs.currentIndexChanged[QtCore.QString].connect(self.apply_figs)
        
    def apply_figs(self,key):
        self.figesq.setPixmap(QtGui.QPixmap(_fromUtf8("interface_figs/"+dics[str(key)][0]+".png")))
        self.figdir.setPixmap(QtGui.QPixmap(_fromUtf8("interface_figs/"+dics[str(key)][1]+".png")))
        self.grid.clear()
        self.grid.addItems(map(lambda s: str(s),dics[str(key)][2]))
        
        
         
    def browse_folder(self):
        directory = QtGui.QFileDialog.getOpenFileName(self,
                                                      "Open File",
                                                      '',
                                                      'cfg (*.cfg)')
        # execute getExistingDirectory dialog and set the directory variable to be equal
        # to the user selected directory

        if directory: # if user didn't pick a directory don't continue
            for file_name in os.listdir(directory): # for all files, if any, in the directory
                self.listWidget.addItem(file_name)  # add file to the listWidget


def main():
    app = QtGui.QApplication(sys.argv)  # A new instance of QApplication
    form = ExampleApp()  # We set the form to be our ExampleApp (design)
    form.show()  # Show the form
    app.exec_()  # and execute the app


if __name__ == '__main__':  # if we're running file directly and not importing it
    main()  # run the main function
