import os
import sys
import pymol

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from config import *
from window import *

if __name__ == '__main__':
	if os.name == 'nt':
		import ctypes
		myappid = "BIG.Dockey.Dockey.{}".format(DOCKEY_VERSION)
		ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

	app = QApplication(sys.argv)
	win = DockeyMainWindow()

	win.show()

	sys.exit(app.exec())