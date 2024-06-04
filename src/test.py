from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

app = QApplication([])

QCoreApplication.setOrganizationName("Dulab")
QCoreApplication.setOrganizationDomain("big.cdu.edu.cn")
QCoreApplication.setApplicationName("Test")
QSettings.setDefaultFormat(QSettings.IniFormat)

settings = QSettings()
settings.setValue('abc', [1,2,3])
print(settings.value('abc'))

print('#######')
print(settings.value('dddddd'))

win = QMainWindow()
win.show()
app.exit(app.exec())
