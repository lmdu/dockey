import os
import sys
import rc_icons
import multiprocessing

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from config import *
from window import *

class DockeyApplication(QApplication):
	osx_open_with = Signal(str)

	def __init__(self, argv):
		super().__init__(argv)
		self.setOrganizationName("BIG")
		self.setOrganizationDomain("dockey.readthedocs.io")
		self.setApplicationName("Dockey")
		self.setApplicationVersion(DOCKEY_VERSION)

	def event(self, event):
		if sys.platform == 'darwin':
			if isinstance(event, QFileOpenEvent):
				self.osx_open_with.emit(event.file())

		return super().event(event)

if __name__ == '__main__':
	multiprocessing.freeze_support()
	QSettings.setDefaultFormat(QSettings.IniFormat)

	#fix incorrect icon display in taskbar on Windows
	if os.name == 'nt':
		import ctypes
		myappid = "BIG.Dockey.Dockey.{}".format(DOCKEY_VERSION)
		ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

	app = DockeyApplication(sys.argv)
	win = DockeyMainWindow()

	#support for macos open with when associated with dock file
	app.osx_open_with.connect(win.create_db_connect)

	#open project file from command line
	args = app.arguments()
	if len(args) > 1:
		if os.path.isfile(args[1]):
			win.create_db_connect(args[1])

	sys.exit(app.exec())