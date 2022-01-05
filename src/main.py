import os
import sys
import rc_icons

from PySide6.QtCore import *
from PySide6.QtWidgets import *

from config import *
from window import *

if __name__ == '__main__':
	QCoreApplication.setOrganizationName("BIG")
	QCoreApplication.setOrganizationDomain("dockey.readthedocs.io")
	QCoreApplication.setApplicationName("Dockey")

	QSettings.setDefaultFormat(QSettings.IniFormat)

	if os.name == 'nt':
		import ctypes
		myappid = "BIG.Dockey.Dockey.{}".format(DOCKEY_VERSION)
		ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

	app = QApplication(sys.argv)
	win = DockeyMainWindow()

	win.show()

	sys.exit(app.exec())