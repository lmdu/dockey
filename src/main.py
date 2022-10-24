import os
import sys
import rc_icons
import multiprocessing

from PySide6.QtCore import *
from PySide6.QtWidgets import *

from config import *
from window import *

if __name__ == '__main__':
	multiprocessing.freeze_support()

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

	#open project file from command line
	args = app.arguments()
	if len(args) > 1:
		pfile = args[1]

		if os.path.isfile(pfile):
			win.create_db_connect(pfile)
	
	sys.exit(app.exec())