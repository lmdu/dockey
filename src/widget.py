from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

__all__ = ['BrowseInput', 'AutodockConfigDialog', 'CreateProjectDialog']

class BrowseInput(QWidget):
	def __init__(self, parent=None, is_file=True):
		super(BrowseInput, self).__init__(parent)
		self.input = QLineEdit(self)
		self.input.setReadOnly(True)
		self.browse = QPushButton(self)
		self.browse.setFlat(True)
		self.browse.setIcon(QIcon("icons/folder.svg"))

		if is_file:
			self.browse.clicked.connect(self.select_file)
		else:
			self.browse.clicked.connect(self.select_folder)

		layout = QHBoxLayout()
		layout.setSpacing(0)
		layout.setContentsMargins(0,0,0,0)
		layout.addWidget(self.input)
		layout.addWidget(self.browse)

		self.setLayout(layout)

	@Slot()
	def select_file(self):
		exe, _ = QFileDialog.getOpenFileName(self)

		if exe:
			self.input.setText(exe)

	@Slot()
	def select_folder(self):
		folder = QFileDialog.getExistingDirectory(self)

		if folder:
			self.input.setText(folder)

	def get_text(self):
		return self.input.text()

	def set_text(self, text):
		self.input.setText(text)

class CreateProjectDialog(QDialog):
	def __init__(self, parent=None):
		super(CreateProjectDialog, self).__init__(parent)
		self.setWindowTitle("Create a New Project")

		self.name_input = QLineEdit(self)
		self.location_input = BrowseInput(self, False)
		tips = (
			'<font color="gray">'
			'This will create a new folder with input name in the '
			'location you selected'
			'</font>'
		)

		self.tip_label = QLabel(tips, self)

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.check_input_name)
		action_box.rejected.connect(self.reject)

		layout = QFormLayout()
		layout.addRow("Name", self.name_input)
		layout.addRow("Location", self.location_input)
		layout.addRow(self.tip_label)
		layout.addRow(action_box)

		self.setLayout(layout)

	@Slot()
	def check_input_name(self):
		name = self.name_input.text()
		folder = self.location_input.get_text()
		if not name:
			return QMessageBox.warning(self, "Input Warning",
				"Please input a project name"
			)

		if not folder:
			return QMessageBox.warning(self, "Input Warning",
				"Please select a directory where to create project"
			)

		if ' ' in name:
			return QMessageBox.critical(self, "Input Error",
				"White-spaces are not allowed in project name"
			)
		
		if any(c in name for c in '\\/:*?"<>|'):
			return QMessageBox.critical(self, "Input Error",
				"Your input name has illegal character"
			)

		if os.path.exists(os.path.join(folder, name)):
			return QMessageBox.critical(self, "Input Error",
				"Your input project name <b>{}</b> exists in folder <b>{}</b>. Please input a new project name".format(
					name, folder
				)
			)

		self.accept()

	@classmethod
	def get_project_folder(cls, parent):
		dlg = cls(parent)
		if dlg.exec():
			name = dlg.name_input.text()
			location = dlg.location_input.get_text()
			return os.path.join(location, name)


class AutodockConfigDialog(QDialog):
	def __init__(self, parent=None):
		super(AutodockConfigDialog, self).__init__(parent)
		self.setWindowTitle("Autodock config")
		self.resize(QSize(500, -1))
		self.settings = QSettings()

		layout = QVBoxLayout()
		self.setLayout(layout)

		autodock = self.settings.value('Tools/autodock_4', '')
		autogrid = self.settings.value('Tools/autogrid_4', '')
		mgltools = self.settings.value('Tools/mgltools', '')

		if not autodock:
			self.autodock_input = BrowseInput(self)
			layout.addWidget(QLabel("Autodock executable"))
			layout.addWidget(self.autodock_input)
		else:
			self.autodock_input = None

		if not autogrid:
			self.autogrid_input = BrowseInput(self)
			layout.addWidget(QLabel("Autogrid executable"))
			layout.addWidget(self.autogrid_input)
		else:
			self.autogrid_input = None

		if not mgltools:
			self.mgltools_input = BrowseInput(self, False)
			layout.addWidget(QLabel("MGLTools directory"))
			layout.addWidget(self.mgltools_input)
		else:
			self.mgltools_input = None

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.save_config)
		action_box.rejected.connect(self.reject)

		layout.addWidget(action_box)

	@Slot()
	def save_config(self):
		configs = {}

		if self.autogrid_input:
			configs['autogrid_4'] = self.autogrid_input.get_text()

		if self.autodock_input:
			configs['autodock_4'] = self.autodock_input.get_text()

		if self.mgltools_input:
			configs['mgltools'] = self.mgltools_input.get_text()

		for k in configs:
			if not configs[k]:
				return

		for k in configs:
			self.settings.setValue('Tools/{}'.format(k), configs[k])

		self.accept()

	@classmethod
	def check_executable(cls, parent):
		settings = QSettings()
		autodock = settings.value('Tools/autodock_4', '')
		autogrid = settings.value('Tools/autogrid_4', '')
		mgltools = settings.value('Tools/mgltools', '')

		if all([autodock, autogrid, mgltools]):
			return True

		dlg = cls(parent)
		return dlg.exec()
