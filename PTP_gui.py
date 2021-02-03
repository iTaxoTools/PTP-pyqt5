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


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

FORM_CLASS,_=loadUiType("limit.ui")



class Main(QMainWindow, FORM_CLASS):
    def __init__(self,parent=None):
        QWidget.__init__(self)
        super(Main,self).__init__(parent)
        self.setWindowIcon(QIcon('PTP.ico'))
        #self.setWindowTitle('PTP')
        self.setupUi(self)
        self.Handel_Buttons()
        quit = QAction("Quit", self)
        quit.triggered.connect(self.closeEvent)



    def closeEvent(self, event):
         close = QMessageBox.question(self, "QUIT", "Are you sure you want to stop the process?",QMessageBox.Yes | QMessageBox.No)
         if close == QMessageBox.Yes:
             event.accept()

         else:
             event.ignore()



    def Handel_Buttons(self):
        self.pushButton.clicked.connect(self.browse_file1)
        self.pushButton_2.clicked.connect(self.browse_file3)
        self.pushButton_3.clicked.connect(self.trigger1)
        self.pushButton_4.clicked.connect(self.trigger2)
        self.pushButton_5.clicked.connect(self.clear)
        self.checkBox.clicked.connect(self.uncheck1)
        self.checkBox_2.clicked.connect(self.uncheck2)

        #self.pushButton_7.clicked.connect(self.takeinputs)



    def browse_file1(self):
        self.browse_file = QFileDialog.getOpenFileName(self, "browse file", directory=".",filter="All Files (*.*)")
        self.lineEdit.setText(QDir.toNativeSeparators(str(self.browse_file[0])))
        return self.browse_file[0]


    def browse_file3(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.Directory)
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            self.lineEdit_2.setText(QDir.toNativeSeparators(str(filenames[0])))

    def uncheck1(self):
        return self.checkBox_2.setChecked(False)

    def uncheck2(self):
        return self.checkBox.setChecked(False)


    def trigger1(self):

        if self.checkBox.isChecked() == True:

            self.download2()
        else:
            self.download1()



    def trigger2(self):

        if self.checkBox.isChecked() == True:
            self.view2()
        else:
            self.view1()



    def download1(self):
        try:
            self.unique= str(int(time.time()))

            open_file= self.lineEdit.text()
            save_file = self.lineEdit_2.text()
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

            strategy=  "H0"
            num_trees= 0

            bsptp = bootstrap_ptp(filename = open_file, ftype = inputformat, reroot = reroot, method = strategy, firstktrees = num_trees)


            pars, settings = bsptp.delimit(spe_rate= -1.0, max_iters = 20000, min_br = 0.0001, whiten= False, strategy = "H0", sprint= False, pvalue= 0.001)

            pp = partitionparser(taxa_order = bsptp.taxa_order, partitions = pars, scale = 500, fileextension= open_file, ptp_status= "ptp")
            pp.summary(fout = os.path.join(save_file, f"PTP_{self.unique}"), bnmi = False, sp_setting = settings)

            if bsptp.numtrees > 1:
                min_no_p, max_no_p, mean_no_p = pp.hpd_numpartitions()
                print("Estimated number of species is between " + repr(min_no_p) + " and " + repr(max_no_p))
                print("Mean: " + repr(mean_no_p))


        except Exception as e:
            print(e)
            QMessageBox.warning(self, "Warning", f"The species demlimitation output  is not obtained because {e}")
            return
        QMessageBox.information(self, "Information", "The species demlimitation results generated successfully")





    def download2(self):

        try:

            open_file= self.lineEdit.text()
            self.unique= str(int(time.time()))
            save_file = self.lineEdit_2.text()
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

        except Exception as e:

            print(e)
            QMessageBox.warning(self, "Warning", f"The species demlimitation output  is not obtained because {e}")
            return

        QMessageBox.information(self, "Information", "The species demlimitation results generated successfully")




    # def takeinputs(self):
    #     comment1, done1 = QInputDialog.getText(
    #          self, 'Input Dialog', 'Enter your first comment:')
    #
    #     comment2, done2 = QInputDialog.getText(
    #        self, 'Input Dialog', 'Enter your second comment:')
    #
    #     save_file= self.lineEdit_2.text()
    #
    #     if done1 and done2:
    #         try:
    #             if os.path.isfile(os.path.join(save_file, "myoutput.PTPMLPartition.spart")):
    #                 fin = open(os.path.join(save_file, "myoutput.PTPMLPartition.spart"), "rt")
    #             fin = open(os.path.join(save_file, "myoutput.PTPhSupportPartition.spart"), "rt")
    #             data = fin.read()
    #             data = data.replace("this is my first comment", str(comment1))
    #             data = data.replace("this is my second comment", str(comment2))
    #             fin.close()
    #             if os.path.isfile(os.path.join(save_file, "myoutput.PTPMLPartition.spart")):
    #                 fin = open(os.path.join(save_file, "myoutput.PTPhSupportPartition.spart"), "wt")
    #
    #             fin = open(os.path.join(save_file, "myoutput.PTPMLPartition.spart"), "wt")
    #             fin.write(data)
    #             fin.close()
    #
    #         except Exception:
    #             QMessageBox.warning(self, "Warning", "The spart file is not generated yet, please do the  analysis first")
    #             return
    #         QMessageBox.information(self, "Information", "spart file is updated successfully")







    def view1(self):

        try:
            save_file = self.lineEdit_2.text()
            pixmap = QPixmap(os.path.join(save_file, f"PTP_{self.unique}_PTPhSupportPartition.png"))

            if os.path.isfile(os.path.join(save_file, f"PTP_{self.unique}_PTPMLPartition.txt")):
                f = open(os.path.join(save_file, f"PTP_{self.unique}_PTPMLPartition.txt"))

            else:
                f = open(os.path.join(save_file, f"PTP_{self.unique}_PTPhSupportPartition.txt"))


            self.scene = QGraphicsScene()
            self.scene.addPixmap(QPixmap(pixmap))
            self.graphicsView_2.setScene(self.scene)
            self.scene_1 = QGraphicsScene()
            mytext1 = QGraphicsSimpleTextItem(f.read())
            self.scene_1.addItem(mytext1)
            self.graphicsView.setScene(self.scene_1)
            f.close()
        except Exception:
            QMessageBox.warning(self, "Warning", "The species demitation view is not obtained, please genetrate the analysis output first")
            return
        QMessageBox.information(self, "Information", "The species demitation result image and partition generated successfully")



    def view2(self):

        try:
            save_file = self.lineEdit_2.text()
            pixmap = QPixmap(os.path.join(save_file, f"bPTP_{self.unique}_PTPhSupportPartition.png"))

            if os.path.isfile(os.path.join(save_file, f"bPTP_{self.unique}_PTPMLPartition.txt")):
                f = open(os.path.join(save_file, f"bPTP_{self.unique}_PTPMLPartition.txt"))

            else:
                f = open(os.path.join(save_file, f"bPTP_{self.unique}_PTPhSupportPartition.txt"))


            self.scene = QGraphicsScene()
            self.scene.addPixmap(QPixmap(pixmap))
            self.graphicsView_2.setScene(self.scene)
            self.scene_1 = QGraphicsScene()
            mytext1 = QGraphicsSimpleTextItem(f.read())
            self.scene_1.addItem(mytext1)
            self.graphicsView.setScene(self.scene_1)
            f.close()
        except Exception:
            QMessageBox.warning(self, "Warning", "The species demitation view is not obtained, please genetrate the analysis output first")
            return
        QMessageBox.information(self, "Information", "The species demitation result image and partition generated successfully")





    def clear(self):
        try:

            self.lineEdit.setText("")
            self.lineEdit_2.setText("")
            self.graphicsView_2.setScene(QGraphicsScene())
            self.graphicsView.setScene(QGraphicsScene())
            self.RadioGroup = QButtonGroup()
            self.RadioGroup.addButton(self.radioButton)
            self.RadioGroup.addButton(self.radioButton_2)
            self.RadioGroup.setExclusive(False)
            self.radioButton.setChecked(False)
            self.radioButton_2.setChecked(False)
            self.RadioGroup.setExclusive(True)
        except Exception:
            QMessageBox.warning(self, "Warning", "The input data is not cleared, Please do manually")
            return
        QMessageBox.information(self, "Information", "Please start a new analysis or close the window")



def main1():

    app=QApplication(sys.argv)
    window=Main()
    window.show()
    app.exec_()


if __name__=='__main__':
    main1()
