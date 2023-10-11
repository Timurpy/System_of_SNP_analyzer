import sys
from PyQt5.QtWidgets import (QWidget, QLabel,
    QComboBox, QApplication)
from PyQt5 import QtGui, QtWidgets, QtCore

import pandas as pd
from load_data import get_aggregated_patient_data_from_excel
from parse_data import get_MAF_info_from_dbSNP
from bs4 import BeautifulSoup as BS
import requests

class ParsingThread(QtCore.QThread):
    mysignal = QtCore.pyqtSignal(list)
    
    def  __init__(self, 
                  patient_snp_list=None, 
                  snp_info_dictionary=None,
                  ):
        super().__init__()
        self.patient_snp_list = patient_snp_list
        self.snp_info_dictionary = snp_info_dictionary
        
    def run(self):
        for num, rs_id in enumerate(self.patient_snp_list):
            
            if rs_id in self.snp_info_dictionary.keys():
                allele_data = self.snp_info_dictionary[rs_id]
            else:
                allele_data = get_MAF_info_from_dbSNP(rs_id)
                self.snp_info_dictionary[rs_id] = allele_data
            self.mysignal.emit([num, rs_id, allele_data, self.snp_info_dictionary])
            
class ComboBoxWithBlock(QComboBox):
    
    def __init__(self, parent):
        QComboBox.__init__(self, parent)
        self.readonly = False
        
    def setReadonly(self, value):
        self.readonly = value

    def mousePressEvent(self, event):
        if not self.readonly:
            QComboBox.mousePressEvent(self, event)

    def keyPressEvent(self, event):
        if not self.readonly:
            QComboBox.keyPressEvent(self, event)

    def wheelEvent(self, event):
        if not self.readonly:
            QComboBox.wheelEvent(self, event)

class Example(QWidget):

    def __init__(self):
        super().__init__()
        
        df, unique_patients = get_aggregated_patient_data_from_excel(config_path="config.yaml")
        self.df = df
        self.unique_patients = sorted(unique_patients, key=lambda x: int(x.split('_')[-1]))
        self.snp_info_dictionary = {}
        
        self.initUI()


    def initUI(self):

        self.lbl = QLabel("opa", self)
        
        self.combo = ComboBoxWithBlock(self)
        self.combo.addItems(self.unique_patients)

        self.combo.move(50, 50)
        self.lbl.move(50, 150)

        self.combo.activated[str].connect(self.onActivated)

        self.setGeometry(300, 300, 700, 500)
        self.setWindowTitle('QComboBox')
        
        self.request_status = QLabel(self)
        self.request_status.setGeometry(320, 100, 500, 20)
        self.request_status.setVisible(False)
        
        self.progress = QtWidgets.QProgressBar(self)
        self.progress.setGeometry(50, 100, 250, 20)
        self.progress.setVisible(False)
        
        
        
        self.show()


    def onActivated(self, patient_id):
        
        self.lbl.setText(' ')
        self.combo.setReadonly(True)
        
        patient_info = self.df[patient_id]
        patient_snp_list = list(patient_info[patient_info>0].index)
        self.allele_data_list = []
        
        self.progress.setRange(0, len(patient_snp_list))
        self.progress.setValue(1)
        self.progress.setVisible(True)
        self.allele_data_list = []
        
        self.request_status.setVisible(True)
        self.request_status.setText(f"starting...")
        
        
        self.mythread = ParsingThread(patient_snp_list, self.snp_info_dictionary)    # Создаем экземпляр класса
        # self.mythread.started.connect(self.on_started)
        self.mythread.finished.connect(self.on_finished)
        self.mythread.mysignal.connect(self.on_change, QtCore.Qt.QueuedConnection)
        self.mythread.start() 
        
    def on_change(self, emmited_list):
        num, rs_id, allele_data, updated_info_dictionary = emmited_list
        self.progress.setValue(num+1)
        self.request_status.setText(f"parsing info on {rs_id} from dbSNP")
        self.request_status.adjustSize()
        if len(allele_data): 
            mafs = [float(entry.split('=')[-1].split()[0]) for entry in allele_data]
            mean_maf = round(sum(mafs)/len(mafs), 5)
            self.allele_data_list.append((rs_id, f"avarage MAF = {mean_maf} (from {len(mafs)} sources)"))
        self.snp_info_dictionary = updated_info_dictionary
        
    def on_finished(self):      # Вызывается при завершении потока
        
        data_to_print = [': '.join(snp_) for snp_ in self.allele_data_list]
        self.progress.setMaximum(1)
        self.progress.setValue(1)
        self.lbl.setText('\n'.join(data_to_print))
        self.lbl.adjustSize()
        self.request_status.setText(f"completed")
        self.combo.setReadonly(False)
        # self.progress.setVisible(False)
        # self.request_status.setVisible(False)

if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())