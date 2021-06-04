import ete3
from ete3 import Tree, SeqGroup
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import re
import sys, os
from PyQt5.uic import loadUiType
from bPTP import bayesianptp
from PTP import bootstrap_ptp
from summary import partitionparser
from PyQt5.QtGui import QPixmap
from PyQt5.QtGui import *
import datetime
import re
import pandas as pd
import time
from collections import defaultdict
import tempfile
from PyQt5.QtWebEngineWidgets import QWebEngineView as QWebView,QWebEnginePage as QWebPage
from PyQt5 import QtWebEngineWidgets
from functools import partial

def factory(arg):
    result= str(arg)
    return result


def_dict = defaultdict(partial(factory, 'default value'))
real= {"None": None, "True": True, "False": False}
def_dict['None']=  None
def_dict['True']=  True
def_dict['False']=  False

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

FORM_CLASS,_=loadUiType(resource_path("ptp_new.ui"))


class Main(QDialog, FORM_CLASS):
    def __init__(self,parent=None):
        super(Main, self).__init__(parent)
        self.setupUi(self)
        self.setWindowIcon(QIcon(resource_path('PTP.ico')))
        self.box= QMessageBox()
        self.setWindowTitle("PTP")
        self.f= tempfile.TemporaryDirectory()
        self.gp.setVisible(False)
        self.reset_placement()
        self.filepath= defaultdict(lambda: None)
        self.outpath= defaultdict(lambda: None)
        self.toolButton_7.setEnabled(False)
        self.toolButton_8.setEnabled(False)
        self.toolButton_9.setEnabled(False)


        quit = QAction("Quit", self)
        quit.triggered.connect(self.closeEvent)

        self.m1.setText('10000')
        self.m2.setText('H1')
        self.m3.setText('100')
        self.m4.setText('0.1')
        self.m5.setText('123')
        self.m6.setText("None")
        self.m7.setText("False")
        self.m8.setText('0')

        self.m1.setValidator(QIntValidator())

        self.m3.setValidator(QIntValidator())

        self.m5.setValidator(QIntValidator())

        self.Handel_Buttons()


    def closeEvent(self, event):
         close = QMessageBox.question(self, "QUIT", "Are you sure you want to close the program?",QMessageBox.Yes | QMessageBox.No)
         if close == QMessageBox.Yes:

             #os.removedirs(temp_directory)
             event.accept()
             sys.exit()
         else:
             event.ignore()


    #======= SETUP slots =================================



    def Handel_Buttons(self):
        self.toolButton_6.clicked.connect(self.open_file)
        self.toolButton_7.clicked.connect(self.trigger1)
        self.toolButton_8.clicked.connect(self.save_file)
        self.toolButton_9.clicked.connect(self.save_all)
        self.toolButton_10.clicked.connect(self.clear)
        self.pushButton.clicked.connect(self.BtnHandler)
        self.checkBox.clicked.connect(self.uncheck1)
        self.checkBox_2.clicked.connect(self.uncheck2)
        self.listWidget.itemDoubleClicked.connect(self.Clicked)

    def open_file(self):
        sel = 'Select non-ultrametric phylogram tree file'
        tab = self.file_dialog(sel, ".")
        if tab:
            print(tab)
            self.filepath['path']= tab
        self.toolButton_7.setEnabled(True)

    def file_dialog(self, msg, path):
        return QFileDialog.getOpenFileName(self, msg, path)[0]

    def download2(self):
        try:
            self.toolButton_8.setEnabled(True)
            self.toolButton_9.setEnabled(True)
            m1= int(self.m1.text())
            m3= int(self.m3.text())
            m4= float(self.m4.text())
            m5= int(self.m5.text())
            m8= int(self.m8.text())
            m2= self.m2.text()
            m6= self.m6.text()
            def_dict.default_factory = partial(factory, m6)
            m7= def_dict[self.m7.text()]
            m6= def_dict[self.m6.text()]

            tmpfname=self.filepath['path']
            origfname=None

            _, filename= os.path.split(tmpfname)
            filename= filename.split(".")[0]
            output= self.f.name
            self.outpath['output']= output
            print(self.outpath['output'])
            open_file= tmpfname
            self.unique= str(int(time.time()))
            save_file = output
            with  open(open_file) as treetest:
                l1 = treetest.readline()

            if self.radioButton.isChecked() == True:
                inputformat = "nexus"
                reroot = False

            elif self.radioButton_2.isChecked() == True:
                inputformat = "raxml"
                reroot = False

            elif self.radioButton_3.isChecked() == True:
                inputformat = "raxml"
                reroot = True

            else:
                inputformat = "raxml"
                reroot = True

            bbptp = bayesianptp(filename=open_file, ftype=inputformat, reroot= reroot)
            bbptp = bayesianptp(filename=open_file,
                                ftype=inputformat,
                                reroot=reroot,
                                method=m2,
                                seed=m5,
                                thinning=m3,
                                sampling=m1,
                                burnin=m4,
                                firstktrees=m8)


            if m6!= None and len(m6) > 0:
                bbptp.remove_outgroups(m6, remove=m7, output=open_file + ".NoOutgroups")
            pars, llhs, settings = bbptp.delimit()
            pp = partitionparser(taxa_order=bbptp.taxa_order, partitions=pars, llhs=llhs, fileextension= open_file, ptp_status= "bptp")
            if bbptp.numtrees == 1:
                pp.summary(fout=os.path.join(save_file, f"bPTP_{self.unique}"),
                ML_par=bbptp.get_maxhhl_partition(),
                ml_spe_setting=bbptp.maxhhlsetting,
                sp_setting=settings)

            else:
                pp.summary(fout=os.path.join(save_file, f"bPTP_{int(time.time())}"), sp_setting=settings)
                min_no_p, max_no_p, mean_no_p = pp.hpd_numpartitions()

            onlyfiles = [self.listWidget.addItem(f) for f in os.listdir(save_file) if os.path.isfile(os.path.join(save_file, f))]

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"The species delimitation output  is not obtained because {e}")
            return
        QMessageBox.information(self, "Information", "The species delimitation results generated successfully")





    def save_file(self):
        try:

            msg = 'Please browse to output folder to save chosen files'
            QMessageBox.information(self, 'Browse output folder', msg)
            dlg = QFileDialog()
            dlg.setFileMode(QFileDialog.Directory)
            if dlg.exec_():
                filenames = dlg.selectedFiles()
                filename= QDir.toNativeSeparators(str(filenames[0]))

            items = self.listWidget.selectedItems()
            values = []
            for i in range(len(items)):
                values.append(str(self.listWidget.selectedItems()[i].text()))

            out_files= [os.path.join(filename, value) for value in values]
            in_files= [os.path.join(self.outpath['output'], value) for value in values]
            import shutil
            [shutil.copy(in_file, out_file) for in_file in in_files for out_file in out_files]

        except Exception as e:
            print(e)
            QMessageBox.warning(self, "Warning", f"The output  is not saved because {e}")


    def save_all(self):
        try:

            from pathlib import Path
            msg = 'Please browse to output folder to save all files'
            QMessageBox.information(self, 'Browse output folder', msg)

            dlg = QFileDialog()
            dlg.setFileMode(QFileDialog.Directory)
            if dlg.exec_():
                filenames = dlg.selectedFiles()
                filename= QDir.toNativeSeparators(str(filenames[0]))
            import shutil
            file_names = os.listdir(self.outpath['output'])
            for file_name in file_names:
                shutil.copy(os.path.join(self.outpath['output'], file_name), filename)

        except Exception as e:
            print(e)
            QMessageBox.warning(self, "Warning", f"The output  is not saved because {e}")


    def clear(self):
        self.toolButton_7.setEnabled(False)
        self.toolButton_8.setEnabled(False)
        self.toolButton_9.setEnabled(False)
        self.m1.setText('10000')
        self.m2.setText('H1')
        self.m3.setText('100')
        self.m4.setText('0.1')
        self.m5.setText('123')
        self.m6.setText("None")
        self.m7.setText("False")
        self.m8.setText('0')
        self.listWidget.clear()
        import os, shutil
        folder = self.f.name
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))


    def BtnHandler(self):
        if self.pushButton.isChecked():
            self.gp.setVisible(True)
            print(self.sizeHint())
        else:
            self.gp.setVisible(False)
            self.reset_placement()


    def reset_placement(self):
        g = QDesktopWidget().availableGeometry()
        self.resize(0.1 * g.width(), 0.4 * g.height())
        self.move(g.center().x() - self.width() / 2, g.center().y() - self.height() / 2)


    def uncheck1(self):
        return self.checkBox_2.setChecked(False)

    def uncheck2(self):
        return self.checkBox.setChecked(False)


    def trigger1(self):

        if self.checkBox.isChecked() == True:

            self.download2()
        elif self.checkBox_2.isChecked() == True:
            self.download1()

        else:
            QMessageBox.warning(self, "Warning", f"The selection should be made in checkbox and radiobuttons")



    def download1(self):
        try:
            self.toolButton_8.setEnabled(True)
            self.toolButton_9.setEnabled(True)
            m1= int(self.m1.text())
            m3= int(self.m3.text())
            m4= float(self.m4.text())
            m5= int(self.m5.text())
            m8= int(self.m8.text())
            m2= self.m2.text()
            m6= self.m6.text()
            def_dict.default_factory = partial(factory, m6)
            m7= def_dict[self.m7.text()]
            m6= def_dict[self.m6.text()]
            tmpfname=self.filepath['path']
            origfname=None
            _, filename= os.path.split(tmpfname)
            filename= filename.split(".")[0]
            output= self.f.name
            self.outpath['output']= output
            print(self.outpath['output'])


            open_file= tmpfname
            self.unique= str(int(time.time()))
            save_file = output

            with  open(open_file) as treetest:
                l1 = treetest.readline()

            if self.radioButton.isChecked() == True:
                inputformat = "nexus"
                reroot = False

            elif self.radioButton_2.isChecked() == True:
                inputformat = "raxml"
                reroot = False

            elif self.radioButton_3.isChecked() == True:
                inputformat = "raxml"
                reroot = True

            else:
                inputformat = "raxml"
                reroot = True

            bsptp = bootstrap_ptp(filename = open_file, ftype = inputformat, reroot = reroot, method = m2, firstktrees = m8)
            pars, settings = bsptp.delimit(spe_rate= -1.0, max_iters = 20000, min_br = 0.0001, whiten= False, strategy = m2, sprint= False, pvalue= 0.001)
            pp = partitionparser(taxa_order = bsptp.taxa_order, partitions = pars, scale = 500, fileextension= open_file, ptp_status= "ptp")
            pp.summary(fout = os.path.join(save_file, f"PTP_{self.unique}"), bnmi = False, sp_setting = settings)
            onlyfiles = [self.listWidget.addItem(f) for f in os.listdir(save_file) if os.path.isfile(os.path.join(save_file, f))]

            if bsptp.numtrees > 1:
                min_no_p, max_no_p, mean_no_p = pp.hpd_numpartitions()
                print("Estimated number of species is between " + repr(min_no_p) + " and " + repr(max_no_p))
                print("Mean: " + repr(mean_no_p))

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"The species delimitation output  is not obtained because {e}")
            return
        QMessageBox.information(self, "Information", "The species delimitation results generated successfully")


    def Clicked(self, item2):
        try:
            self.w= AnotherWindow()
            name= item2.text()
            suffix= (".svg", ".png")
            su= (".txt", ".tre", ".spart")
            if name.endswith(suffix):
                self.w.layout.addWidget(self.w.m_output)
                self.w.setLayout(self.w.layout)
                self.w.m_output.load(QUrl().fromLocalFile(os.path.join(self.outpath['output'], name)))
                self.w.setWindowIcon(QIcon(os.path.join("icon", "PTP.ico")))
                self.w.setWindowTitle("PTP")
                self.w.show()
            elif name.endswith(su):
                self.w= AnotherWindow()
                self.w.setWindowIcon(QIcon(os.path.join("icon", "PTP.ico")))
                f = open(os.path.join(self.outpath['output'], name), "rt")
                mytext1 = QGraphicsSimpleTextItem(f.read())
                self.w.scene.addItem(mytext1)
                f.close()
                self.w.layout.addWidget(self.w.graph1)
                self.w.setLayout(self.w.layout)
                self.w.setWindowIcon(QIcon(os.path.join("icon", "PTP.ico")))
                self.w.setWindowTitle("PTP")
                self.w.show()
            else:
                pass

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"The view  is not obtained because {e}")


class AnotherWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self.graph1= QGraphicsView()
        self.scene = QGraphicsScene()
        self.graph1.setScene(self.scene)
        self.m_output = QtWebEngineWidgets.QWebEngineView()


def main1():

    app=QApplication(sys.argv)
    window=Main()
    window.show()
    QApplication.processEvents()
    app.exec_()


if __name__=='__main__':
    main1()
